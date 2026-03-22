//! Twin Prime Search Engine v4 — GMP-accelerated Rust implementation
//!
//! Architecture:
//! 1. Sieve of Eratosthenes to 10^8 for prime table
//! 2. Two-tier algebraic sieve: base (to 10^6) + extended (to 10^8)
//!    Twin primes form (6m-1, 6m+1), sieve eliminates m where 6m±1 ≡ 0 (mod p)
//! 3. SPRP(2) on p1, then p2 (short-circuit), then SPRP(3,5,7) on survivors
//! 4. Lucas PRP confirmation (completing BPPSW) — no redundant SPRP(2)
//! 5. Rayon parallel iteration across batches and candidates
//!
//! Key optimizations over v3:
//! - GMP (via rug) for all big-integer arithmetic: Montgomery modpow, Toom-Cook multiply
//! - mpz_fdiv_ui for single-limb mod in sieve (no BigUint allocation)
//! - Optimized test ordering: SPRP(2) p1→p2, then SPRP(3,5,7), then Lucas
//! - GMP's native Jacobi symbol and primality testing

use rand::Rng;
use rayon::prelude::*;
use rug::integer::IsPrime;
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

/// Run base sieve (primes to 10^6) on a batch starting at m_start.
/// Uses GMP's mpz_fdiv_ui for fast single-limb mod of big integers.
fn base_sieve(alive: &mut [bool], m_start: &Integer, sieve: &SieveData) {
    let batch_size = alive.len();
    for (idx, &p) in sieve.primes_small.iter().enumerate() {
        let inv6 = sieve.inv6_small[idx];
        // GMP's fdiv_ui: single-limb mod, no allocation
        let m_mod_p = m_start.mod_u(p as u32) as u64;
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p + inv6 - m_mod_p) as usize
        };
        let mut j = r1;
        while j < batch_size {
            unsafe { *alive.get_unchecked_mut(j) = false; }
            j += p_us;
        }

        let complement = p - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p + complement - m_mod_p) as usize
        };
        let mut j = r2;
        while j < batch_size {
            unsafe { *alive.get_unchecked_mut(j) = false; }
            j += p_us;
        }
    }
}

/// Run extended sieve (primes 10^6 to 10^8) on a batch.
fn extended_sieve(alive: &mut [bool], m_start: &Integer, sieve: &SieveData) {
    let batch_size = alive.len();
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
                unsafe { *alive.get_unchecked_mut(j) = false; }
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
                unsafe { *alive.get_unchecked_mut(j) = false; }
                j += p_us;
            }
        }
    }
}

/// SPRP (Strong Probable Prime) test using GMP's modpow (Montgomery reduction).
fn is_sprp(n: &Integer, base: u32) -> bool {
    if *n < 2 {
        return false;
    }
    let nm1 = Integer::from(n - 1);

    // Write n-1 = 2^r * d
    let r = nm1.find_one(0).unwrap_or(0);
    let d = Integer::from(&nm1 >> r);

    let base_int = Integer::from(base);
    let mut x = base_int.pow_mod(&d, n).unwrap();

    if x == 1 || x == nm1 {
        return true;
    }

    for _ in 0..r - 1 {
        x = Integer::from(x.pow_mod_ref(&Integer::from(2), n).unwrap());
        if x == nm1 {
            return true;
        }
        if x == 1 {
            return false;
        }
    }
    false
}

/// Test if m yields a twin prime pair (6m-1, 6m+1).
/// Optimized test ordering: SPRP(2) both → SPRP(3,5,7) both → GMP full primality.
fn test_candidate(m_val: u64, m_start: &Integer) -> Option<(Integer, Integer)> {
    let m = Integer::from(m_start + m_val);
    let p1 = Integer::from(&m * 6) - 1;

    // Stage 1: SPRP(2) on p1 — cheapest filter, rejects ~50%
    if !is_sprp(&p1, 2) {
        return None;
    }

    let p2 = Integer::from(&p1 + 2);

    // Stage 2: SPRP(2) on p2 — short-circuit: test both with base 2 first
    if !is_sprp(&p2, 2) {
        return None;
    }

    // Stage 3: SPRP(3,5,7) on p1 — further composite filter
    for &base in &[3u32, 5, 7] {
        if !is_sprp(&p1, base) {
            return None;
        }
    }

    // Stage 4: SPRP(3,5,7) on p2
    for &base in &[3u32, 5, 7] {
        if !is_sprp(&p2, base) {
            return None;
        }
    }

    // Stage 5: GMP full primality (BPPSW internally, no redundant SPRP(2))
    // is_probably_prime with reps=25 runs GMP's internal BPPSW + extra Miller-Rabin
    if p1.is_probably_prime(25) != IsPrime::No && p2.is_probably_prime(25) != IsPrime::No {
        return Some((p1, p2));
    }
    None
}

fn get_params(target_digits: u32) -> (usize, bool) {
    match target_digits {
        0..=150 => (500_000, false),
        151..=500 => (2_000_000, true),
        501..=1500 => (5_000_000, true),
        _ => (10_000_000, true),
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

    let ln_n = (target_digits as f64) * std::f64::consts::LN_10;
    let hl_raw = ln_n * ln_n / (2.0 * 0.6601618158);

    let sieve_mode = if use_extended { "full" } else { "base only" };
    println!(
        "  Target: ~{} digits, HL trials: {:.0}, Batch: {}, Sieve: {}",
        target_digits, hl_raw, batch_size, sieve_mode
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
        let mut rng = rand::thread_rng();
        let batch_starts: Vec<Integer> = (0..num_workers)
            .map(|_| {
                let offset = rng.gen_range(0u64..m_range_u64);
                Integer::from(&m_low + offset)
            })
            .collect();

        let flag = found_flag.clone();
        let results: Vec<_> = batch_starts
            .into_par_iter()
            .map(|m_start| {
                if flag.load(Ordering::Relaxed) {
                    return (0u64, 0u64, 0u64, None);
                }

                let mut alive = vec![true; batch_size];
                base_sieve(&mut alive, &m_start, sieve);
                if use_extended {
                    extended_sieve(&mut alive, &m_start, sieve);
                }

                let survivors: Vec<u64> = alive
                    .iter()
                    .enumerate()
                    .filter(|(_, &a)| a)
                    .map(|(i, _)| i as u64)
                    .collect();

                let n_surv = survivors.len() as u64;
                let mut tested = 0u64;

                for &off in &survivors {
                    if flag.load(Ordering::Relaxed) {
                        break;
                    }
                    tested += 1;
                    if let Some((p1, p2)) = test_candidate(off, &m_start) {
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
        println!(
            "    [{:.1}s] {} raw -> {} sieved -> {} tested ({:.0}/s)",
            el, total_raw, total_surv, total_tested, rate
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
    println!("Twin Prime Engine v4 (Rust + GMP)");
    println!(
        "Rayon x{} | SPRP(2,3,5,7)+GMP BPPSW | Montgomery modpow | Toom-Cook multiply",
        num_workers
    );
    println!();

    let t_pre = Instant::now();
    let sieve = SieveData::build(100_000_000);
    let setup_time = t_pre.elapsed().as_secs_f64();
    println!(
        "Sieve: {} primes to 10^6, {} to 10^8 ({:.2}s)",
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

    // v2 Python baselines
    let v2_times: Vec<(u32, f64)> =
        vec![(100, 2.23), (500, 11.98), (1000, 75.28), (2000, 1554.0)];
    // v3 Rust (num-bigint) baselines
    let v3_times: Vec<(u32, f64)> =
        vec![(100, 0.03), (500, 2.13), (1000, 13.74), (2000, 631.0)];

    let mut all_results: Vec<SearchResult> = Vec::new();
    let mut largest: Option<(Integer, Integer, SearchResult)> = None;

    for &(td, budget) in &targets {
        println!("{}", "=".repeat(65));
        let (result, primes) = search_twins(td, budget, &sieve);

        if result.found {
            println!("  FOUND in {:.2}s!", result.elapsed_secs);
            if let (Some(ref h), Some(ref t)) = (&result.p_head, &result.p_tail) {
                println!("  ({}...{}, +2)", h, t);
            }
            if let Some(d) = result.p_digits {
                println!("  {} digits", d);
            }
            println!(
                "  Pipeline: {} raw -> {} sieved -> {} tested",
                result.total_raw, result.total_surv, result.total_tested
            );
            if let Some(&(_, v2t)) = v2_times.iter().find(|&&(d, _)| d == td) {
                println!("  vs Python v2: {:.1}x faster", v2t / result.elapsed_secs);
            }
            if let Some(&(_, v3t)) = v3_times.iter().find(|&&(d, _)| d == td) {
                println!(
                    "  vs Rust v3 (num-bigint): {:.1}x faster",
                    v3t / result.elapsed_secs
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
        println!("{}", "=".repeat(65));
        println!("LARGEST TWIN PRIME FOUND");
        println!("{}", "=".repeat(65));
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
        "version": "v4-rust-gmp",
        "workers": num_workers,
        "results": all_results,
    });
    if let Ok(json) = serde_json::to_string_pretty(&output) {
        std::fs::write("twin_prime_engine_results.json", json).ok();
    }
    println!();
    println!("Saved to twin_prime_engine_results.json");
}
