use serde::{Deserialize, Serialize};

#[derive(Clone, Debug)]
pub struct KPlan {
    pub p: u32,
    pub inv6: u32,
    pub bad1: u32,
    pub bad2: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SurvivorBatch {
    pub n: u64,
    pub k_start: u64,
    pub k_batch_size: usize,
    pub sieve_limit: usize,
    pub batch_index: u64,
    pub batch_k0: u64,
    pub estimated_digits: u64,
    pub survival_rate: f64,
    pub total_input_k: usize,
    pub survivors: Vec<u64>,
}

pub fn align_k_start(k_start: u64) -> u64 {
    if k_start <= 3 {
        return 3;
    }
    let rem = k_start % 6;
    if rem == 3 {
        k_start
    } else {
        k_start + ((3 + 6 - rem) % 6)
    }
}

pub fn estimate_digits(k: u64, n: u64) -> u64 {
    ((k as f64).log10() + (n as f64) * std::f64::consts::LOG10_2).floor() as u64 + 1
}

pub fn batch_k0(k_start: u64, k_batch_size: usize, batch_index: u64) -> Option<u64> {
    let batch_k0_u128 = k_start as u128 + 6u128 * k_batch_size as u128 * batch_index as u128;
    if batch_k0_u128 > u64::MAX as u128 {
        None
    } else {
        Some(batch_k0_u128 as u64)
    }
}

pub fn survival_rate(plan: &[KPlan]) -> f64 {
    plan.iter()
        .fold(1.0, |acc, entry| acc * ((entry.p - 2) as f64 / entry.p as f64))
}

pub fn prime_sieve(limit: usize) -> Vec<u32> {
    let num_bytes = limit / 8 + 1;
    let mut sieve = vec![0xFFu8; num_bytes];
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

pub fn pow_mod_u64(base: u64, exp: u64, modu: u32) -> u32 {
    let modu_u64 = modu as u64;
    let mut result: u64 = 1;
    let mut base_acc = base % modu_u64;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result = ((result as u128 * base_acc as u128) % modu_u64 as u128) as u64;
        }
        base_acc = ((base_acc as u128 * base_acc as u128) % modu_u64 as u128) as u64;
        e >>= 1;
    }
    result as u32
}

pub fn mod_inv_prime(a: u32, p: u32) -> u32 {
    pow_mod_u64(a as u64, (p - 2) as u64, p)
}

pub fn build_k_sieve_plan(n: u64, sieve_limit: usize) -> Vec<KPlan> {
    prime_sieve(sieve_limit)
        .into_iter()
        .filter(|&p| p >= 5)
        .map(|p| {
            let two_n_mod_p = pow_mod_u64(2, n, p);
            let inv_two_n = mod_inv_prime(two_n_mod_p, p);
            let neg_inv = if inv_two_n == 0 { 0 } else { p - inv_two_n };
            let (bad1, bad2) = if inv_two_n <= neg_inv {
                (inv_two_n, neg_inv)
            } else {
                (neg_inv, inv_two_n)
            };
            KPlan {
                p,
                inv6: mod_inv_prime(6, p),
                bad1,
                bad2,
            }
        })
        .collect()
}

pub fn build_k_sieve_plan_range(n: u64, low_exclusive: usize, high_inclusive: usize) -> Vec<KPlan> {
    if high_inclusive <= low_exclusive {
        return Vec::new();
    }

    prime_sieve(high_inclusive)
        .into_iter()
        .filter(|&p| p >= 5 && (p as usize) > low_exclusive)
        .map(|p| {
            let two_n_mod_p = pow_mod_u64(2, n, p);
            let inv_two_n = mod_inv_prime(two_n_mod_p, p);
            let neg_inv = if inv_two_n == 0 { 0 } else { p - inv_two_n };
            let (bad1, bad2) = if inv_two_n <= neg_inv {
                (inv_two_n, neg_inv)
            } else {
                (neg_inv, inv_two_n)
            };
            KPlan {
                p,
                inv6: mod_inv_prime(6, p),
                bad1,
                bad2,
            }
        })
        .collect()
}

fn set_bit(alive: &mut [u64], idx: usize) {
    alive[idx >> 6] |= 1u64 << (idx & 63);
}

fn clear_bit(alive: &mut [u64], idx: usize) {
    alive[idx >> 6] &= !(1u64 << (idx & 63));
}

pub fn bit_is_set(alive: &[u64], idx: usize) -> bool {
    alive[idx >> 6] & (1u64 << (idx & 63)) != 0
}

pub fn self_prime_exception_indices(
    n: u64,
    p: u32,
    k_start: u64,
    k_count: usize,
) -> Vec<usize> {
    if n >= 63 {
        return Vec::new();
    }

    let two_n = 1u128 << n;
    if two_n > p as u128 + 1 {
        return Vec::new();
    }

    let mut out = Vec::new();
    for numer in [p as u128 - 1, p as u128 + 1] {
        if numer == 0 || numer % two_n != 0 {
            continue;
        }
        let k = numer / two_n;
        if k == 0 || k % 6 != 3 {
            continue;
        }
        let k_u64 = k as u64;
        if k_u64 < k_start {
            continue;
        }
        let delta = k_u64 - k_start;
        if delta % 6 != 0 {
            continue;
        }
        let idx = (delta / 6) as usize;
        if idx < k_count {
            out.push(idx);
        }
    }
    out
}

pub fn is_self_prime_exception(n: u64, p: u32, k: u64) -> bool {
    if n >= 63 {
        return false;
    }

    let two_n = 1u128 << n;
    for numer in [p as u128 - 1, p as u128 + 1] {
        if numer == 0 || numer % two_n != 0 {
            continue;
        }
        let k_exact = numer / two_n;
        if k_exact == k as u128 {
            return true;
        }
    }
    false
}

pub fn candidate_survives_plan(n: u64, k: u64, plan: &[KPlan]) -> bool {
    for entry in plan {
        let k_mod_p = (k % entry.p as u64) as u32;
        if (k_mod_p == entry.bad1 || k_mod_p == entry.bad2) && !is_self_prime_exception(n, entry.p, k) {
            return false;
        }
    }
    true
}

pub fn filter_survivors_with_plan(n: u64, survivors: &[u64], plan: &[KPlan]) -> Vec<u64> {
    if plan.is_empty() {
        return survivors.to_vec();
    }
    survivors
        .iter()
        .copied()
        .filter(|&k| candidate_survives_plan(n, k, plan))
        .collect()
}

pub fn sieve_k_batch(n: u64, k_start: u64, k_count: usize, plan: &[KPlan]) -> Vec<u64> {
    let words = (k_count + 63) / 64;
    let mut alive = vec![u64::MAX; words];
    if k_count & 63 != 0 {
        let valid = (1u64 << (k_count & 63)) - 1;
        let last = alive.len() - 1;
        alive[last] = valid;
    }

    for entry in plan {
        let p_u64 = entry.p as u64;
        let k_mod_p = k_start % p_u64;
        for residue in [entry.bad1, entry.bad2] {
            let i0 =
                (((residue as u64 + p_u64 - k_mod_p) % p_u64) * entry.inv6 as u64) % p_u64;
            let mut idx = i0 as usize;
            while idx < k_count {
                clear_bit(&mut alive, idx);
                idx += entry.p as usize;
            }
        }

        for idx in self_prime_exception_indices(n, entry.p, k_start, k_count) {
            set_bit(&mut alive, idx);
        }
    }

    alive
}

pub fn survivor_offsets(alive: &[u64], k_count: usize) -> Vec<u64> {
    let mut out = Vec::new();
    for (word_idx, &word) in alive.iter().enumerate() {
        if word == 0 {
            continue;
        }
        let base = (word_idx as u64) << 6;
        let mut w = word;
        while w != 0 {
            let bit = w.trailing_zeros() as u64;
            let off = base + bit;
            if off >= k_count as u64 {
                break;
            }
            out.push(off);
            w &= w - 1;
        }
    }
    out
}

pub fn generate_survivor_batch(
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    batch_index: u64,
) -> Option<SurvivorBatch> {
    let plan = build_k_sieve_plan(n, sieve_limit);
    generate_survivor_batch_with_plan(n, k_start, k_batch_size, sieve_limit, batch_index, &plan)
}

pub fn generate_survivor_batch_with_plan(
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    batch_index: u64,
    plan: &[KPlan],
) -> Option<SurvivorBatch> {
    let k_start = align_k_start(k_start);
    let batch_k0 = batch_k0(k_start, k_batch_size, batch_index)?;
    let alive = sieve_k_batch(n, batch_k0, k_batch_size, plan);
    let survivors = survivor_offsets(&alive, k_batch_size)
        .into_iter()
        .map(|off| batch_k0 + 6 * off)
        .collect();
    Some(SurvivorBatch {
        n,
        k_start,
        k_batch_size,
        sieve_limit,
        batch_index,
        batch_k0,
        estimated_digits: estimate_digits(batch_k0, n),
        survival_rate: survival_rate(plan),
        total_input_k: k_batch_size,
        survivors,
    })
}

pub fn formula_minus(k: u64, n: u64) -> String {
    format!("{}*2^{}-1", k, n)
}

pub fn formula_plus(k: u64, n: u64) -> String {
    format!("{}*2^{}+1", k, n)
}
