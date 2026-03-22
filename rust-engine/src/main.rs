//! Twin Prime Search Engine v3 — Rust implementation (pure Rust, no GMP)
//!
//! Architecture:
//! 1. Sieve of Eratosthenes to 10^8 for prime table
//! 2. Two-tier algebraic sieve: base (to 10^6) + extended (to 10^8)
//!    Twin primes have form (6m-1, 6m+1), sieve eliminates m where 6m±1 ≡ 0 (mod p)
//! 3. Multi-base SPRP(2,3,5,7) filter on survivors
//! 4. Full BPPSW confirmation (SPRP(2) + Lucas test)
//! 5. Rayon parallel iteration across batches and candidates

use num_bigint::BigUint;
use num_integer::Integer as NumInteger;
use num_traits::{One, Zero};
use rand::Rng;
use rayon::prelude::*;
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
fn base_sieve(alive: &mut [bool], m_start: &BigUint, sieve: &SieveData) {
    let batch_size = alive.len();
    for (idx, &p) in sieve.primes_small.iter().enumerate() {
        let inv6 = sieve.inv6_small[idx];
        let p_big = BigUint::from(p);
        let m_mod_p = (m_start % &p_big).to_u64_digits();
        let m_mod_p = if m_mod_p.is_empty() { 0u64 } else { m_mod_p[0] };
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p + inv6 - m_mod_p) as usize
        };
        let mut j = r1;
        while j < batch_size {
            alive[j] = false;
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
            alive[j] = false;
            j += p_us;
        }
    }
}

/// Run extended sieve (primes 10^6 to 10^8) on a batch.
fn extended_sieve(alive: &mut [bool], m_start: &BigUint, sieve: &SieveData) {
    let batch_size = alive.len();
    for (idx, &p) in sieve.primes_ext.iter().enumerate() {
        let inv6 = sieve.inv6_ext[idx];
        let p_big = BigUint::from(p);
        let m_mod_p = (m_start % &p_big).to_u64_digits();
        let m_mod_p = if m_mod_p.is_empty() { 0u64 } else { m_mod_p[0] };
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p + inv6 - m_mod_p) as usize
        };
        if r1 < batch_size {
            let mut j = r1;
            while j < batch_size {
                alive[j] = false;
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
                alive[j] = false;
                j += p_us;
            }
        }
    }
}

/// Modular exponentiation: base^exp mod modulus
fn mod_pow(base: &BigUint, exp: &BigUint, modulus: &BigUint) -> BigUint {
    base.modpow(exp, modulus)
}

/// SPRP (Strong Probable Prime) test for a single base.
fn is_sprp(n: &BigUint, base: u32) -> bool {
    let one = BigUint::one();
    let two = BigUint::from(2u32);

    if *n < two {
        return false;
    }

    let nm1 = n - &one;

    // Write n-1 = 2^r * d
    let mut d = nm1.clone();
    let mut r: u64 = 0;
    while d.is_even() {
        d >>= 1u32;
        r += 1;
    }

    let base_big = BigUint::from(base);
    let mut x = mod_pow(&base_big, &d, n);

    if x == one || x == nm1 {
        return true;
    }

    for _ in 0..r - 1 {
        x = mod_pow(&x, &two, n);
        if x == nm1 {
            return true;
        }
        if x == one {
            return false;
        }
    }
    false
}

/// Jacobi symbol (a/n)
fn jacobi(a: &BigUint, n: &BigUint) -> i32 {
    let zero = BigUint::zero();
    let one = BigUint::one();

    if *n == one {
        return 1;
    }

    let mut a = a % n;
    let mut n = n.clone();
    let mut result = 1i32;

    while a != zero {
        while a.is_even() {
            a >>= 1u32;
            let n_mod8 = &n % BigUint::from(8u32);
            let n_mod8_val = if n_mod8.to_u64_digits().is_empty() {
                0u64
            } else {
                n_mod8.to_u64_digits()[0]
            };
            if n_mod8_val == 3 || n_mod8_val == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        let a_mod4 = &a % BigUint::from(4u32);
        let n_mod4 = &n % BigUint::from(4u32);
        let a4 = if a_mod4.to_u64_digits().is_empty() {
            0u64
        } else {
            a_mod4.to_u64_digits()[0]
        };
        let n4 = if n_mod4.to_u64_digits().is_empty() {
            0u64
        } else {
            n_mod4.to_u64_digits()[0]
        };
        if a4 == 3 && n4 == 3 {
            result = -result;
        }
        a = &a % &n;
    }

    if n == one {
        result
    } else {
        0
    }
}

/// Lucas probable prime test (part of BPPSW).
fn is_lucas_prp(n: &BigUint) -> bool {
    let one = BigUint::one();
    let two = BigUint::from(2u32);

    // Find first D in {5, -7, 9, -11, ...} with Jacobi(D, n) = -1
    let mut d_val: i64 = 5;
    let mut sign = 1i64;
    let d_big;
    loop {
        let d_abs = BigUint::from(d_val.unsigned_abs());
        let j = if sign > 0 {
            jacobi(&d_abs, n)
        } else {
            // Jacobi(-d, n) = Jacobi(-1, n) * Jacobi(d, n)
            let n_mod4 = n % BigUint::from(4u32);
            let n4 = if n_mod4.to_u64_digits().is_empty() {
                0u64
            } else {
                n_mod4.to_u64_digits()[0]
            };
            let neg1_jac = if n4 == 1 { 1 } else { -1 };
            neg1_jac * jacobi(&d_abs, n)
        };
        if j == -1 {
            d_big = if sign > 0 {
                d_abs
            } else {
                // D is negative; we'll handle P,Q differently
                d_abs
            };
            break;
        }
        if j == 0 {
            let d_bu = BigUint::from(d_val.unsigned_abs());
            if d_bu > two && &d_bu != n {
                return false;
            }
        }
        d_val += 2;
        sign = -sign;
    }

    // P = 1, Q = (1 - D) / 4
    // For the standard strong Lucas test, we use P=1.
    // With D found, use P=1, Q=(1-D)/4
    // Since D alternates sign, handle carefully:
    // If sign > 0: D = d_val, Q = (1 - d_val) / 4
    // If sign < 0: D = -d_val, Q = (1 + d_val) / 4

    let actual_d = d_val * sign;
    let q_val = (1 - actual_d) / 4;

    // Compute U_{n+1} mod n using the standard Lucas chain
    let nm1 = n + &one; // n+1 for the Lucas test

    // Binary expansion of n+1
    let bits = nm1.bits();

    // Lucas sequence: U_k, V_k via doubling formulas
    let n_big = n;

    // We track U_k and V_k mod n
    let mut u = BigUint::one();
    let mut v = BigUint::one(); // P = 1

    let q_abs = BigUint::from(q_val.unsigned_abs());
    let q_neg = q_val < 0;

    let mut q_k = if q_neg {
        // Q mod n = n - |Q|
        if q_abs < *n_big {
            n_big - &q_abs
        } else {
            n_big - &(&q_abs % n_big)
        }
    } else {
        q_abs.clone() % n_big
    };

    for i in (0..bits - 1).rev() {
        // Double: U_{2k} = U_k * V_k mod n
        let u_new = (&u * &v) % n_big;
        // V_{2k} = V_k^2 - 2*Q^k mod n
        let v_sq = (&v * &v) % n_big;
        let two_qk = (&q_k + &q_k) % n_big;
        let v_new = if v_sq >= two_qk {
            (&v_sq - &two_qk) % n_big
        } else {
            (n_big - &((&two_qk - &v_sq) % n_big)) % n_big
        };
        q_k = (&q_k * &q_k) % n_big;

        u = u_new;
        v = v_new;

        if nm1.bit(i) {
            // Advance: U_{2k+1} = (P*U_{2k} + V_{2k}) / 2
            // V_{2k+1} = (D*U_{2k} + P*V_{2k}) / 2
            // With P=1:
            // U' = (U + V) / 2 mod n
            // V' = (D*U + V) / 2 mod n

            let u_plus_v = (&u + &v) % n_big;
            let u_next = if u_plus_v.is_even() {
                u_plus_v >> 1u32
            } else {
                (&u_plus_v + n_big) >> 1u32
            };

            // D*U mod n
            let d_abs_big = &d_big;
            let du = if sign > 0 {
                (d_abs_big * &u) % n_big
            } else {
                let tmp = (d_abs_big * &u) % n_big;
                if tmp.is_zero() {
                    BigUint::zero()
                } else {
                    n_big - &tmp
                }
            };
            let du_plus_v = (&du + &v) % n_big;
            let v_next = if du_plus_v.is_even() {
                du_plus_v >> 1u32
            } else {
                (&du_plus_v + n_big) >> 1u32
            };

            q_k = (&q_k * &q_abs) % n_big;
            if q_neg {
                q_k = if q_k.is_zero() {
                    BigUint::zero()
                } else {
                    n_big - &q_k
                };
            }

            u = u_next;
            v = v_next;
        }
    }

    // n is a Lucas PRP if U_{n+1} ≡ 0 (mod n)
    u.is_zero()
}

/// Full BPPSW test: SPRP(2) + Lucas PRP.
fn is_bppsw(n: &BigUint) -> bool {
    is_sprp(n, 2) && is_lucas_prp(n)
}

/// Test if m yields a twin prime pair (6m-1, 6m+1).
fn test_candidate(m_val: u64, m_start: &BigUint) -> Option<(BigUint, BigUint)> {
    let m = m_start + BigUint::from(m_val);
    let six = BigUint::from(6u32);
    let p1 = &m * &six - BigUint::one();

    // Multi-base SPRP filter (cheap)
    for &base in &[2u32, 3, 5, 7] {
        if !is_sprp(&p1, base) {
            return None;
        }
    }

    let p2 = &p1 + BigUint::from(2u32);
    for &base in &[2u32, 3, 5, 7] {
        if !is_sprp(&p2, base) {
            return None;
        }
    }

    // Full BPPSW confirmation
    if is_bppsw(&p1) && is_bppsw(&p2) {
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
) -> (SearchResult, Option<(BigUint, BigUint)>) {
    let (batch_size, use_extended) = get_params(target_digits);

    let ten = BigUint::from(10u32);
    let m_low = ten.pow(target_digits - 1) / BigUint::from(6u32);
    let m_high = ten.pow(target_digits) / BigUint::from(6u32);
    let m_range = &m_high - &m_low - BigUint::from(batch_size as u64);

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

    // Convert m_range to u64 if it fits, otherwise use max
    let m_range_u64 = {
        let digits = m_range.to_u64_digits();
        if digits.len() <= 1 {
            digits.first().copied().unwrap_or(1)
        } else {
            u64::MAX
        }
    };

    while t0.elapsed().as_secs_f64() < max_seconds && !found_flag.load(Ordering::Relaxed) {
        let mut rng = rand::thread_rng();
        let batch_starts: Vec<BigUint> = (0..num_workers)
            .map(|_| {
                let offset = BigUint::from(rng.gen_range(0u64..m_range_u64));
                &m_low + &offset
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
    println!("Twin Prime Engine v3 (Rust)");
    println!(
        "Rayon x{} | SPRP(2,3,5,7)+BPPSW | Adaptive sieve | Native",
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

    let v2_times: Vec<(u32, f64)> = vec![(100, 2.23), (500, 11.98), (1000, 75.28), (2000, 1554.0)];

    let mut all_results: Vec<SearchResult> = Vec::new();
    let mut largest: Option<(BigUint, BigUint, SearchResult)> = None;

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
                let speedup = v2t / result.elapsed_secs;
                println!("  vs Python v2: {:.1}x faster", speedup);
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
        "version": "v3-rust",
        "workers": num_workers,
        "results": all_results,
    });
    if let Ok(json) = serde_json::to_string_pretty(&output) {
        std::fs::write("twin_prime_engine_results.json", json).ok();
    }
    println!();
    println!("Saved to twin_prime_engine_results.json");
}
