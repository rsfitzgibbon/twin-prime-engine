//! Twin Prime Search Engine v5 — Optimized GMP-accelerated Rust implementation
//!
//! v5 optimizations over v4:
//! - Packed bitset sieve (8× smaller working set → fits L2 cache, better locality)
//! - In-place rug::Integer operations (zero allocation in primality testing hot path)
//! - Reusable SprpCtx/TestCtx buffers per worker (no per-candidate malloc/free)
//! - BPPSW-only final confirmation (is_probably_prime(0), no extra Miller-Rabin)
//! - Non-overlapping sequential batches per round (no redundant coverage)
//! - Adaptive batch sizing per digit target
//!
//! Architecture:
//! 1. Sieve of Eratosthenes to 10^8 for prime table
//! 2. Two-tier algebraic sieve on packed bitset: base (to 10^6) + extended (to 10^8)
//!    Twin primes form (6m-1, 6m+1), sieve eliminates m where 6m±1 ≡ 0 (mod p)
//! 3. SPRP(2) on p1, then p2 (short-circuit), then SPRP(3,5,7) on survivors
//! 4. BPPSW confirmation (Lucas PRP via GMP) — deterministic for all known composites
//! 5. Rayon parallel iteration across batches and candidates

use rand::Rng;
use rayon::prelude::*;
use rug::integer::IsPrime;
use rug::Assign;
use rug::Integer;
use serde::Serialize;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Instant;

/// Sieve of Eratosthenes returning all primes up to limit.
fn prime_sieve(limit: usize) -> Vec<u64> {
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    if limit > 0 {
        is_prime[1] = false;
    }
    let sqrt_limit = (limit as f64).sqrt() as usize;
    for p in 2..=sqrt_limit {
        if is_prime[p] {
            let mut j = p * p;
            while j <= limit {
                is_prime[j] = false;
                j += p;
            }
        }
    }
    is_prime
        .iter()
        .enumerate()
        .filter(|(_, &b)| b)
        .map(|(i, _)| i as u64)
        .collect()
}

/// Compute modular inverse of 6 mod p using Fermat's little theorem.
fn inv6_mod(p: u64) -> u64 {
    let mut result: u64 = 1;
    let mut base: u64 = 6 % p;
    let mut exp = p - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % p as u128) as u64;
        }
        base = ((base as u128 * base as u128) % p as u128) as u64;
        exp >>= 1;
    }
    result
}

struct SieveData {
    primes_small: Vec<u64>,
    inv6_small: Vec<u64>,
    primes_ext: Vec<u64>,
    inv6_ext: Vec<u64>,
}

impl SieveData {
    fn build(sieve_limit: usize) -> Self {
        let all_primes = prime_sieve(sieve_limit);
        let primes_small: Vec<u64> = all_primes
            .iter()
            .copied()
            .filter(|&p| p >= 5 && p <= 1_000_000)
            .collect();
        let inv6_small: Vec<u64> = primes_small.iter().map(|&p| inv6_mod(p)).collect();
        let primes_ext: Vec<u64> = all_primes
            .iter()
            .copied()
            .filter(|&p| p > 1_000_000)
            .collect();
        let inv6_ext: Vec<u64> = primes_ext.iter().map(|&p| inv6_mod(p)).collect();
        SieveData {
            primes_small,
            inv6_small,
            primes_ext,
            inv6_ext,
        }
    }
}

/// Run base sieve (primes to 10^6) on a packed bitset.
/// Each bit in alive[] represents one candidate m value.
fn base_sieve(alive: &mut [u64], batch_size: usize, m_start: &Integer, sieve: &SieveData) {
    for (idx, &p) in sieve.primes_small.iter().enumerate() {
        let inv6 = sieve.inv6_small[idx];
        let m_mod_p = m_start.mod_u(p as u32) as u64;
        let p_us = p as usize;

        // Offset where 6(m_start + r1) - 1 ≡ 0 (mod p)
        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p + inv6 - m_mod_p) as usize
        };
        let mut j = r1;
        while j < batch_size {
            unsafe {
                *alive.get_unchecked_mut(j >> 6) &= !(1u64 << (j & 63));
            }
            j += p_us;
        }

        // Offset where 6(m_start + r2) + 1 ≡ 0 (mod p)
        let complement = p - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p + complement - m_mod_p) as usize
        };
        let mut j = r2;
        while j < batch_size {
            unsafe {
                *alive.get_unchecked_mut(j >> 6) &= !(1u64 << (j & 63));
            }
            j += p_us;
        }
    }
}

/// Run extended sieve (primes 10^6 to 10^8) on a packed bitset.
fn extended_sieve(
    alive: &mut [u64],
    batch_size: usize,
    m_start: &Integer,
    sieve: &SieveData,
) {
    for (idx, &p) in sieve.primes_ext.iter().enumerate() {
        let inv6 = sieve.inv6_ext[idx];
        let m_mod_p = m_start.mod_u(p as u32) as u64;
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p + inv6 - m_mod_p) as usize
        };
        if r1 < batch_size {
            let mut j = r1;
            while j < batch_size {
                unsafe {
                    *alive.get_unchecked_mut(j >> 6) &= !(1u64 << (j & 63));
                }
                j += p_us;
            }
        }

        let complement = p - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p + complement - m_mod_p) as usize
        };
        if r2 < batch_size {
            let mut j = r2;
            while j < batch_size {
                unsafe {
                    *alive.get_unchecked_mut(j >> 6) &= !(1u64 << (j & 63));
                }
                j += p_us;
            }
        }
    }
}

/// Collect survivor indices from packed bitset using bit tricks.
fn collect_survivors(alive: &[u64], batch_size: usize) -> Vec<u64> {
    let mut result = Vec::new();
    for (word_idx, &word) in alive.iter().enumerate() {
        if word == 0 {
            continue;
        }
        let base = (word_idx as u64) << 6;
        let mut w = word;
        while w != 0 {
            let bit = w.trailing_zeros() as u64;
            let idx = base + bit;
            if idx < batch_size as u64 {
                result.push(idx);
            }
            w &= w - 1; // clear lowest set bit
        }
    }
    result
}

/// Pre-allocated buffers for SPRP testing — eliminates per-test allocations.
struct SprpCtx {
    nm1: Integer,
    d: Integer,
    x: Integer,
    two: Integer,
}

impl SprpCtx {
    fn new() -> Self {
        SprpCtx {
            nm1: Integer::new(),
            d: Integer::new(),
            x: Integer::new(),
            two: Integer::from(2),
        }
    }

    /// Strong Probable Prime test using in-place GMP operations.
    /// Zero allocation: all computation uses pre-allocated buffers.
    fn is_sprp(&mut self, n: &Integer, base: u32) -> bool {
        // nm1 = n - 1
        self.nm1.assign(n);
        self.nm1 -= 1;

        // Write n-1 = 2^r * d
        let r = self.nm1.find_one(0).unwrap_or(0);
        self.d.assign(&self.nm1);
        self.d >>= r;

        // x = base^d mod n (in-place: assigns result to self.x)
        self.x.assign(base);
        self.x.pow_mod_mut(&self.d, n).unwrap();

        if self.x == 1 || self.x == self.nm1 {
            return true;
        }

        for _ in 0..r - 1 {
            // x = x^2 mod n (in-place squaring)
            self.x.pow_mod_mut(&self.two, n).unwrap();
            if self.x == self.nm1 {
                return true;
            }
            if self.x == 1 {
                return false;
            }
        }
        false
    }
}

/// Per-worker context with pre-allocated buffers for all operations.
/// Created once per Rayon task, reused across all candidates in the batch.
struct TestCtx {
    sprp: SprpCtx,
    p1: Integer,
    p2: Integer,
}

impl TestCtx {
    fn new() -> Self {
        TestCtx {
            sprp: SprpCtx::new(),
            p1: Integer::new(),
            p2: Integer::new(),
        }
    }

    /// Test if m yields a twin prime pair (6m-1, 6m+1).
    /// Optimized ordering: SPRP(2) both -> SPRP(3,5,7) both -> BPPSW confirmation.
    /// All arithmetic is in-place: zero heap allocation until a pair is found.
    fn test_candidate(&mut self, m_val: u64, m_start: &Integer) -> Option<(Integer, Integer)> {
        // p1 = 6*(m_start + m_val) - 1
        self.p1.assign(m_start);
        self.p1 += m_val;
        self.p1 *= 6;
        self.p1 -= 1;

        // Stage 1: SPRP(2) on p1 — cheapest filter, rejects ~65% of sieve survivors
        if !self.sprp.is_sprp(&self.p1, 2) {
            return None;
        }

        // p2 = p1 + 2
        self.p2.assign(&self.p1);
        self.p2 += 2;

        // Stage 2: SPRP(2) on p2 — short-circuit: test both with base 2 first
        if !self.sprp.is_sprp(&self.p2, 2) {
            return None;
        }

        // Stage 3: SPRP(3,5,7) on p1 — further composite filter
        for &base in &[3u32, 5, 7] {
            if !self.sprp.is_sprp(&self.p1, base) {
                return None;
            }
        }

        // Stage 4: SPRP(3,5,7) on p2
        for &base in &[3u32, 5, 7] {
            if !self.sprp.is_sprp(&self.p2, base) {
                return None;
            }
        }

        // Stage 5: BPPSW confirmation via GMP (is_probably_prime(0) = SPRP(2) + Lucas)
        // No known BPPSW pseudoprimes exist; combined with our SPRP(2,3,5,7) pre-filter,
        // the false positive probability is vanishingly small.
        if self.p1.is_probably_prime(0) != IsPrime::No
            && self.p2.is_probably_prime(0) != IsPrime::No
        {
            return Some((self.p1.clone(), self.p2.clone()));
        }
        None
    }
}

/// Adaptive parameters: batch size and sieve mode per digit target.
/// Larger batches amortize the per-prime mod_u cost in the sieve phase.
/// All batch sizes are multiples of 64 for bitset alignment.
fn get_params(target_digits: u32) -> (usize, bool) {
    match target_digits {
        0..=150 => (512_000, false),       // 500K, base sieve only (sieve-dominated)
        151..=500 => (2_048_000, true),     // 2M, full sieve
        501..=1500 => (5_120_000, true),    // 5M, full sieve
        1501..=3000 => (10_240_000, true),  // 10M, full sieve (balanced)
        _ => (20_480_000, true),            // 20M, full sieve (testing-dominated)
    }
}

#[derive(Serialize, Clone)]
struct SearchResult {
    found: bool,
    digits: u32,
    elapsed_secs: f64,
    total_raw: u64,
    total_surv: u64,
    total_tested: u64,
    batches: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    p_head: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    p_tail: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    p_digits: Option<usize>,
}

fn search_twins(
    target_digits: u32,
    max_seconds: f64,
    sieve: &SieveData,
) -> (SearchResult, Option<(Integer, Integer)>) {
    let (batch_size, use_extended) = get_params(target_digits);
    let words = batch_size / 64;

    let m_low = {
        let mut x = Integer::from(Integer::u_pow_u(10, target_digits - 1));
        x /= 6;
        x
    };
    let m_high = {
        let mut x = Integer::from(Integer::u_pow_u(10, target_digits));
        x /= 6;
        x
    };
    let m_range = Integer::from(&m_high - &m_low) - batch_size as u64;

    // Hardy-Littlewood expected twin prime trials for this digit count
    let ln_n = (target_digits as f64) * std::f64::consts::LN_10;
    let hl_raw = ln_n * ln_n / (2.0 * 0.6601618158);

    let sieve_mode = if use_extended { "full" } else { "base only" };
    println!(
        "  Target: ~{} digits | HL est: {:.0} trials | Batch: {} | Sieve: {} | Bitset: {}KB",
        target_digits,
        hl_raw,
        batch_size,
        sieve_mode,
        words * 8 / 1024
    );

    let t0 = Instant::now();
    let found_flag = Arc::new(AtomicBool::new(false));
    let mut total_raw: u64 = 0;
    let mut total_surv: u64 = 0;
    let mut total_tested: u64 = 0;
    let mut batches: u32 = 0;
    let num_workers = num_cpus::get();
    let m_range_u64 = m_range.to_u64().unwrap_or(u64::MAX);

    while t0.elapsed().as_secs_f64() < max_seconds && !found_flag.load(Ordering::Relaxed) {
        // Non-overlapping batch starts: random base, sequential stride
        let mut rng = rand::thread_rng();
        let round_base = rng.gen_range(0u64..m_range_u64);
        let stride = batch_size as u64;
        let batch_starts: Vec<Integer> = (0..num_workers)
            .map(|i| {
                let offset = (round_base + (i as u64) * stride) % m_range_u64;
                Integer::from(&m_low + offset)
            })
            .collect();

        let flag = found_flag.clone();
        let deadline = max_seconds;
        let results: Vec<_> = batch_starts
            .into_par_iter()
            .map(|m_start| {
                if flag.load(Ordering::Relaxed) {
                    return (0u64, 0u64, 0u64, None);
                }

                // Sieve phase: packed bitset (8x smaller than bool array)
                let mut alive = vec![u64::MAX; words];
                base_sieve(&mut alive, batch_size, &m_start, sieve);
                if use_extended {
                    extended_sieve(&mut alive, batch_size, &m_start, sieve);
                }

                let survivors = collect_survivors(&alive, batch_size);
                let n_surv = survivors.len() as u64;

                // Test phase: zero-allocation in-place operations
                let mut ctx = TestCtx::new();
                let mut tested = 0u64;

                for &off in &survivors {
                    if flag.load(Ordering::Relaxed) {
                        break;
                    }
                    // Check time budget every 64 candidates (amortize syscall)
                    if tested & 63 == 0 && t0.elapsed().as_secs_f64() > deadline {
                        break;
                    }
                    tested += 1;
                    if let Some((p1, p2)) = ctx.test_candidate(off, &m_start) {
                        flag.store(true, Ordering::Relaxed);
                        return (batch_size as u64, n_surv, tested, Some((p1, p2)));
                    }
                }

                (batch_size as u64, n_surv, tested, None)
            })
            .collect();

        for (raw, surv, tested, found) in results {
            total_raw += raw;
            total_surv += surv;
            total_tested += tested;
            batches += 1;

            if let Some((p1, p2)) = found {
                let elapsed = t0.elapsed().as_secs_f64();
                let p_str = p1.to_string();
                let p_len = p_str.len();
                let p_head = p_str[..40.min(p_len)].to_string();
                let p_tail = p_str[p_len.saturating_sub(40)..].to_string();

                let result = SearchResult {
                    found: true,
                    digits: target_digits,
                    elapsed_secs: elapsed,
                    total_raw,
                    total_surv,
                    total_tested,
                    batches,
                    p_head: Some(p_head),
                    p_tail: Some(p_tail),
                    p_digits: Some(p_len),
                };
                return (result, Some((p1, p2)));
            }
        }

        let el = t0.elapsed().as_secs_f64();
        let rate = if el > 0.0 {
            total_tested as f64 / el
        } else {
            0.0
        };
        let surv_pct = if total_raw > 0 {
            100.0 * total_surv as f64 / total_raw as f64
        } else {
            0.0
        };
        println!(
            "    [{:.1}s] {} raw -> {} sieved ({:.2}%) -> {} tested ({:.0}/s)",
            el, total_raw, total_surv, surv_pct, total_tested, rate
        );
    }

    let elapsed = t0.elapsed().as_secs_f64();
    let result = SearchResult {
        found: false,
        digits: target_digits,
        elapsed_secs: elapsed,
        total_raw,
        total_surv,
        total_tested,
        batches,
        p_head: None,
        p_tail: None,
        p_digits: None,
    };
    (result, None)
}

fn main() {
    let num_workers = num_cpus::get();
    println!("Twin Prime Engine v5 (Rust + GMP, optimized)");
    println!(
        "Rayon x{} | In-place SPRP(2,3,5,7) + BPPSW | Bitset sieve | Montgomery modpow",
        num_workers
    );
    println!();

    let t_pre = Instant::now();
    let sieve = SieveData::build(100_000_000);
    let setup_time = t_pre.elapsed().as_secs_f64();
    println!(
        "Sieve: {} base primes (to 10^6), {} extended (to 10^8) [{:.2}s]",
        sieve.primes_small.len(),
        sieve.primes_ext.len(),
        setup_time
    );
    println!();

    let targets: Vec<(u32, f64)> = vec![
        (100, 30.0),
        (500, 60.0),
        (1000, 300.0),
        (2000, 1800.0),
        (5000, 3600.0),
    ];

    // Baselines for comparison
    let v2_py: Vec<(u32, f64)> = vec![(100, 2.23), (500, 11.98), (1000, 75.28), (2000, 1554.0)];
    let v3_rs: Vec<(u32, f64)> = vec![(100, 0.03), (500, 2.13), (1000, 13.74), (2000, 631.0)];
    let v4_rs: Vec<(u32, f64)> = vec![(2000, 882.0)];

    let mut all_results: Vec<SearchResult> = Vec::new();
    let mut largest: Option<(Integer, Integer, SearchResult)> = None;

    for &(td, budget) in &targets {
        println!("{}", "=".repeat(70));
        let (result, primes) = search_twins(td, budget, &sieve);

        if result.found {
            println!("  FOUND in {:.2}s!", result.elapsed_secs);
            if let (Some(ref h), Some(ref t)) = (&result.p_head, &result.p_tail) {
                println!("  ({}...{}, +2)", h, t);
            }
            if let Some(d) = result.p_digits {
                println!("  {} actual digits", d);
            }
            println!(
                "  Pipeline: {} raw -> {} sieved -> {} tested",
                result.total_raw, result.total_surv, result.total_tested
            );
            if let Some(&(_, t)) = v2_py.iter().find(|&&(d, _)| d == td) {
                println!("  vs Python v2: {:.1}x faster", t / result.elapsed_secs);
            }
            if let Some(&(_, t)) = v3_rs.iter().find(|&&(d, _)| d == td) {
                println!(
                    "  vs Rust v3 (num-bigint): {:.1}x faster",
                    t / result.elapsed_secs
                );
            }
            if let Some(&(_, t)) = v4_rs.iter().find(|&&(d, _)| d == td) {
                println!(
                    "  vs v4 (GMP baseline): {:.1}x faster",
                    t / result.elapsed_secs
                );
            }
            if let Some((p1, p2)) = primes {
                let is_larger = match &largest {
                    Some((_, _, ref lr)) => result.p_digits > lr.p_digits,
                    None => true,
                };
                if is_larger {
                    largest = Some((p1, p2, result.clone()));
                }
            }
        } else {
            println!(
                "  Timeout ({:.1}s), {} tested in {} batches",
                result.elapsed_secs, result.total_tested, result.batches
            );
        }

        all_results.push(result);
    }

    if let Some((ref p1, _, ref lr)) = largest {
        println!();
        println!("{}", "=".repeat(70));
        println!("LARGEST TWIN PRIME FOUND");
        println!("{}", "=".repeat(70));
        let sp = p1.to_string();
        println!("  Digits: {}", sp.len());
        println!("  p = {}", &sp[..70.min(sp.len())]);
        if sp.len() > 140 {
            println!("      ...");
        }
        println!("      {}", &sp[sp.len().saturating_sub(70)..]);
        println!("  q = p + 2");
        println!("  Time: {:.2}s", lr.elapsed_secs);
    }

    let output = serde_json::json!({
        "version": "v5-rust-gmp-optimized",
        "workers": num_workers,
        "results": all_results,
    });
    if let Ok(json) = serde_json::to_string_pretty(&output) {
        std::fs::write("twin_prime_engine_results.json", json).ok();
    }
    println!();
    println!("Results saved to twin_prime_engine_results.json");
}
