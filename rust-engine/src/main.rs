//! Twin Prime Search Engine v5.1 — Optimized GMP-accelerated Rust implementation
//!
//! v5.1 optimizations over v5:
//! - Bitset Eratosthenes (8× less memory for prime generation → enables deeper sieve)
//! - Three-tier sieve: base (10^6) + extended (10^8) + deep (2×10^8) for high digits
//! - Compact u32 prime/inverse tables (84.5 MB cache vs ~169 MB with u64)
//! - Inline survivor testing: iterate bitset directly, no Vec<u64> allocation
//! - Hardware popcount for survivor counting
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
//! 1. Bitset Sieve of Eratosthenes to 2×10^8 for prime table
//! 2. Three-tier algebraic sieve on packed bitset:
//!    base (to 10^6) + extended (to 10^8) + deep (to 2×10^8)
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

/// Bitset Sieve of Eratosthenes — 8× less memory than bool array.
/// Uses packed u8 bitset: bit i of byte[i/8] represents number i.
fn prime_sieve(limit: usize) -> Vec<u32> {
    let num_bytes = limit / 8 + 1;
    let mut sieve = vec![0xFFu8; num_bytes];
    // Clear bits 0 and 1 (not prime)
    sieve[0] &= 0b11111100;

    let sqrt_limit = (limit as f64).sqrt() as usize;
    for p in 2..=sqrt_limit {
        if sieve[p >> 3] & (1u8 << (p & 7)) != 0 {
            let mut j = p * p;
            while j <= limit {
                sieve[j >> 3] &= !(1u8 << (j & 7));
                j += p;
            }
        }
    }

    // Collect primes using bit tricks (faster than per-number iteration)
    let est = limit / ((limit as f64).ln() as usize).max(1);
    let mut primes = Vec::with_capacity(est);
    for byte_idx in 0..num_bytes {
        let byte = sieve[byte_idx];
        if byte == 0 {
            continue;
        }
        let base = byte_idx * 8;
        let mut b = byte;
        while b != 0 {
            let bit = b.trailing_zeros() as usize;
            let num = base + bit;
            if num > limit {
                break;
            }
            if num >= 2 {
                primes.push(num as u32);
            }
            b &= b - 1;
        }
    }
    primes
}

/// Compute modular inverse of 6 mod p using Fermat's little theorem.
fn inv6_mod(p: u32) -> u32 {
    let p64 = p as u64;
    let mut result: u64 = 1;
    let mut base: u64 = 6 % p64;
    let mut exp = p64 - 2;
    while exp > 0 {
        if exp & 1 == 1 {
            result = ((result as u128 * base as u128) % p64 as u128) as u64;
        }
        base = ((base as u128 * base as u128) % p64 as u128) as u64;
        exp >>= 1;
    }
    result as u32
}

struct SieveData {
    primes_small: Vec<u32>,
    inv6_small: Vec<u32>,
    primes_ext: Vec<u32>,
    inv6_ext: Vec<u32>,
    primes_deep: Vec<u32>,
    inv6_deep: Vec<u32>,
}

impl SieveData {
    fn build(sieve_limit: usize) -> Self {
        let all_primes = prime_sieve(sieve_limit);
        let primes_small: Vec<u32> = all_primes
            .iter()
            .copied()
            .filter(|&p| p >= 5 && p <= 1_000_000)
            .collect();
        let inv6_small: Vec<u32> = primes_small.iter().map(|&p| inv6_mod(p)).collect();
        let primes_ext: Vec<u32> = all_primes
            .iter()
            .copied()
            .filter(|&p| p > 1_000_000 && p <= 100_000_000)
            .collect();
        let inv6_ext: Vec<u32> = primes_ext.iter().map(|&p| inv6_mod(p)).collect();
        let primes_deep: Vec<u32> = all_primes
            .iter()
            .copied()
            .filter(|&p| p > 100_000_000)
            .collect();
        let inv6_deep: Vec<u32> = primes_deep.iter().map(|&p| inv6_mod(p)).collect();
        SieveData {
            primes_small,
            inv6_small,
            primes_ext,
            inv6_ext,
            primes_deep,
            inv6_deep,
        }
    }

    fn cache_bytes(&self) -> usize {
        (self.primes_small.len()
            + self.inv6_small.len()
            + self.primes_ext.len()
            + self.inv6_ext.len()
            + self.primes_deep.len()
            + self.inv6_deep.len())
            * std::mem::size_of::<u32>()
    }
}

/// Run base sieve (primes to 10^6) on a packed bitset.
/// Each bit in alive[] represents one candidate m value.
fn base_sieve(alive: &mut [u64], batch_size: usize, m_start: &Integer, sieve: &SieveData) {
    for (idx, &p) in sieve.primes_small.iter().enumerate() {
        let p_u64 = p as u64;
        let inv6 = sieve.inv6_small[idx] as u64;
        let m_mod_p = m_start.mod_u(p) as u64;
        let p_us = p as usize;

        // Offset where 6(m_start + r1) - 1 ≡ 0 (mod p)
        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p_u64 + inv6 - m_mod_p) as usize
        };
        let mut j = r1;
        while j < batch_size {
            unsafe {
                *alive.get_unchecked_mut(j >> 6) &= !(1u64 << (j & 63));
            }
            j += p_us;
        }

        // Offset where 6(m_start + r2) + 1 ≡ 0 (mod p)
        let complement = p_u64 - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p_u64 + complement - m_mod_p) as usize
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
        let p_u64 = p as u64;
        let inv6 = sieve.inv6_ext[idx] as u64;
        let m_mod_p = m_start.mod_u(p) as u64;
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p_u64 + inv6 - m_mod_p) as usize
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

        let complement = p_u64 - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p_u64 + complement - m_mod_p) as usize
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

/// Run deep sieve (primes 10^8 to 2×10^8) on a packed bitset.
/// Used for high-digit targets (≥1500 digits) where SPRP testing dominates.
fn deep_sieve(
    alive: &mut [u64],
    batch_size: usize,
    m_start: &Integer,
    sieve: &SieveData,
) {
    for (idx, &p) in sieve.primes_deep.iter().enumerate() {
        let p_u64 = p as u64;
        let inv6 = sieve.inv6_deep[idx] as u64;
        let m_mod_p = m_start.mod_u(p) as u64;
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p_u64 + inv6 - m_mod_p) as usize
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

        let complement = p_u64 - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p_u64 + complement - m_mod_p) as usize
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

/// Adaptive parameters: batch size and sieve depth per digit target.
/// Larger batches amortize the per-prime mod_u cost in the sieve phase.
/// All batch sizes are multiples of 64 for bitset alignment.
/// Sieve depth: 0 = base only (to 10^6), 1 = + extended (to 10^8), 2 = + deep (to 2×10^8)
fn get_params(target_digits: u32) -> (usize, u8) {
    match target_digits {
        0..=150 => (512_000, 0),        // 500K, base sieve only (sieve-dominated)
        151..=500 => (2_048_000, 1),     // 2M, extended sieve
        501..=1500 => (5_120_000, 1),    // 5M, extended sieve
        1501..=3000 => (10_240_000, 2),  // 10M, deep sieve (testing-dominated)
        _ => (20_480_000, 2),            // 20M, deep sieve (testing-dominated)
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
    let (batch_size, sieve_depth) = get_params(target_digits);
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

    let sieve_mode = match sieve_depth {
        0 => "base only",
        1 => "extended",
        _ => "deep",
    };
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
        let depth = sieve_depth;
        let results: Vec<_> = batch_starts
            .into_par_iter()
            .map(|m_start| {
                if flag.load(Ordering::Relaxed) {
                    return (0u64, 0u64, 0u64, None);
                }

                // Sieve phase: packed bitset (8x smaller than bool array)
                let mut alive = vec![u64::MAX; words];
                base_sieve(&mut alive, batch_size, &m_start, sieve);
                if depth >= 1 {
                    extended_sieve(&mut alive, batch_size, &m_start, sieve);
                }
                if depth >= 2 {
                    deep_sieve(&mut alive, batch_size, &m_start, sieve);
                }

                // Count survivors via hardware popcount (no Vec allocation)
                let n_surv: u64 = alive.iter().map(|&w| w.count_ones() as u64).sum();

                // Test phase: iterate bitset directly, zero-allocation in-place operations
                let mut ctx = TestCtx::new();
                let mut tested = 0u64;

                'outer: for (word_idx, &word) in alive.iter().enumerate() {
                    if word == 0 {
                        continue;
                    }
                    let base = (word_idx as u64) << 6;
                    let mut w = word;
                    while w != 0 {
                        if flag.load(Ordering::Relaxed) {
                            break 'outer;
                        }
                        // Check time budget every 64 candidates (amortize syscall)
                        if tested & 63 == 0 && t0.elapsed().as_secs_f64() > deadline {
                            break 'outer;
                        }
                        let bit = w.trailing_zeros() as u64;
                        let off = base + bit;
                        if off >= batch_size as u64 {
                            break;
                        }
                        tested += 1;
                        if let Some((p1, p2)) = ctx.test_candidate(off, &m_start) {
                            flag.store(true, Ordering::Relaxed);
                            return (batch_size as u64, n_surv, tested, Some((p1, p2)));
                        }
                        w &= w - 1; // clear lowest set bit
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
    println!("Twin Prime Engine v5.1 (Rust + GMP, optimized)");
    println!(
        "Rayon x{} | In-place SPRP(2,3,5,7) + BPPSW | 3-tier sieve | Montgomery modpow",
        num_workers
    );
    println!();

    let t_pre = Instant::now();
    let sieve = SieveData::build(200_000_000);
    let setup_time = t_pre.elapsed().as_secs_f64();
    let sieve_cache_mb = sieve.cache_bytes() as f64 / (1024.0 * 1024.0);
    println!("Prime cache: {:.1} MB", sieve_cache_mb);
    println!(
        "Sieve: {} base (to 10^6), {} extended (to 10^8), {} deep (to 2×10^8) [{:.2}s]",
        sieve.primes_small.len(),
        sieve.primes_ext.len(),
        sieve.primes_deep.len(),
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
        "version": "v5.1-rust-gmp-3tier-sieve",
        "workers": num_workers,
        "results": all_results,
    });
    if let Ok(json) = serde_json::to_string_pretty(&output) {
        std::fs::write("twin_prime_engine_results.json", json).ok();
    }
    println!();
    println!("Results saved to twin_prime_engine_results.json");
}
