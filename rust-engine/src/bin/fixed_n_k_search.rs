use rayon::prelude::*;
use rug::integer::IsPrime;
use rug::Assign;
use rug::Integer;
use serde::{Deserialize, Serialize};
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::Instant;

const C2: f64 = 0.6601618158;
const PROTH_BASES: [u32; 8] = [2, 3, 5, 7, 11, 13, 17, 19];

fn arg_value(args: &[String], flag: &str) -> Option<String> {
    args.windows(2)
        .find(|w| w[0] == flag)
        .map(|w| w[1].clone())
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|arg| arg == flag)
}

fn align_k_start(k_start: u64) -> u64 {
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

fn prime_sieve(limit: usize) -> Vec<u32> {
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

fn pow_mod_u64(base: u64, exp: u64, modu: u32) -> u32 {
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

fn mod_inv_prime(a: u32, p: u32) -> u32 {
    pow_mod_u64(a as u64, (p - 2) as u64, p)
}

#[derive(Clone)]
struct KPlan {
    p: u32,
    inv6: u32,
    bad1: u32,
    bad2: u32,
}

fn build_k_sieve_plan(n: u64, sieve_limit: usize) -> Vec<KPlan> {
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

fn survival_rate(plan: &[KPlan]) -> f64 {
    plan.iter()
        .fold(1.0, |acc, entry| acc * ((entry.p - 2) as f64 / entry.p as f64))
}

fn set_bit(alive: &mut [u64], idx: usize) {
    alive[idx >> 6] |= 1u64 << (idx & 63);
}

fn clear_bit(alive: &mut [u64], idx: usize) {
    alive[idx >> 6] &= !(1u64 << (idx & 63));
}

fn bit_is_set(alive: &[u64], idx: usize) -> bool {
    alive[idx >> 6] & (1u64 << (idx & 63)) != 0
}

fn self_prime_exception_indices(n: u64, p: u32, k_start: u64, k_count: usize) -> Vec<usize> {
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

fn sieve_k_batch(n: u64, k_start: u64, k_count: usize, plan: &[KPlan]) -> Vec<u64> {
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
            let i0 = (((residue as u64 + p_u64 - k_mod_p) % p_u64) * entry.inv6 as u64) % p_u64;
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

fn survivor_offsets(alive: &[u64], k_count: usize) -> Vec<u64> {
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

struct SprpCtx {
    nm1: Integer,
    d: Integer,
    x: Integer,
    two: Integer,
}

impl SprpCtx {
    fn new() -> Self {
        Self {
            nm1: Integer::new(),
            d: Integer::new(),
            x: Integer::new(),
            two: Integer::from(2),
        }
    }

    fn is_sprp(&mut self, n: &Integer, base: u32) -> bool {
        self.nm1.assign(n);
        self.nm1 -= 1;

        let r = self.nm1.find_one(0).unwrap_or(0);
        self.d.assign(&self.nm1);
        self.d >>= r;

        self.x.assign(base);
        self.x.pow_mod_mut(&self.d, n).unwrap();
        if self.x == 1 || self.x == self.nm1 {
            return true;
        }

        for _ in 0..r - 1 {
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

struct FoundPair {
    p: Integer,
    q: Integer,
    plus_mode: String,
}

struct TestCtx {
    sprp: SprpCtx,
    plus: Integer,
    minus: Integer,
    exp: Integer,
}

impl TestCtx {
    fn new() -> Self {
        Self {
            sprp: SprpCtx::new(),
            plus: Integer::new(),
            minus: Integer::new(),
            exp: Integer::new(),
        }
    }

    fn test_candidate(&mut self, k: u64, n: u64) -> (bool, Option<FoundPair>) {
        self.plus.assign(k);
        self.plus <<= n as usize;
        self.plus += 1;

        if !self.sprp.is_sprp(&self.plus, 2) {
            return (false, None);
        }

        let mut plus_mode = None;
        if k & 1 == 1 && (64 - k.leading_zeros()) as u64 <= n {
            self.exp.assign(k);
            self.exp <<= (n - 1) as usize;
            for base in PROTH_BASES {
                self.sprp.x.assign(base);
                self.sprp.x.pow_mod_mut(&self.exp, &self.plus).unwrap();
                if self.sprp.x == self.sprp.nm1 {
                    plus_mode = Some(format!("proth_base_{}", base));
                    break;
                }
            }
        }

        if plus_mode.is_none() && self.plus.is_probably_prime(0) == IsPrime::No {
            return (false, None);
        }

        self.minus.assign(&self.plus);
        self.minus -= 2;
        if !self.sprp.is_sprp(&self.minus, 2) {
            return (true, None);
        }
        if self.minus.is_probably_prime(0) == IsPrime::No {
            return (true, None);
        }

        (
            true,
            Some(FoundPair {
                p: self.minus.clone(),
                q: self.plus.clone(),
                plus_mode: plus_mode.unwrap_or_else(|| "bpsw_fallback".to_string()),
            }),
        )
    }
}

#[derive(Serialize)]
struct SearchResult {
    found: bool,
    n: u64,
    estimated_digits: u64,
    elapsed_secs: f64,
    batches: u64,
    total_sieved: u64,
    tested_plus: u64,
    tested_minus: u64,
    k_start: String,
    k_batch_size: usize,
    sieve_limit: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    k: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    digits: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    p_head: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    p_tail: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    q_head: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    q_tail: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    plus_mode: Option<String>,
}

#[derive(Serialize)]
struct ValidationResult {
    n: u64,
    k_start: u64,
    k_count: usize,
    sieve_limit: usize,
    predicted_survivors: usize,
    exact_survivors: usize,
    mismatch_count: usize,
}

#[derive(Serialize)]
struct SieveOnlyResult {
    n: u64,
    batch_k0: String,
    k_batch_size: usize,
    sieve_limit: usize,
    estimated_digits: u64,
    survival_rate: f64,
    survivors: usize,
    preview_k: Vec<String>,
}

#[derive(Serialize)]
struct ExportSummary {
    n: u64,
    k_start: String,
    k_batch_size: usize,
    sieve_limit: usize,
    batch_start: u64,
    batch_count: u64,
    total_input_k: u64,
    total_survivors: u64,
    export_prefix: String,
    k_file: String,
    minus_file: String,
    plus_file: String,
    pairs_file: String,
}

#[derive(Serialize, Deserialize)]
struct SearchCheckpoint {
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    next_batch: u64,
    total_sieved: u64,
    tested_plus: u64,
    tested_minus: u64,
    total_finds: u64,
}

#[derive(Serialize)]
struct FoundPairRecord {
    batch_index: u64,
    k: String,
    digits: usize,
    plus_mode: String,
    p_head: String,
    p_tail: String,
    q_head: String,
    q_tail: String,
}

#[derive(Serialize)]
struct ContinuousSearchSummary {
    n: u64,
    k_start: String,
    k_batch_size: usize,
    sieve_limit: usize,
    elapsed_secs: f64,
    start_batch: u64,
    next_batch: u64,
    batches_completed: u64,
    total_sieved: u64,
    tested_plus: u64,
    tested_minus: u64,
    total_finds: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    results_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    checkpoint_path: Option<String>,
}

#[derive(Serialize, Deserialize)]
struct ExportCheckpoint {
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    next_batch: u64,
    exported_groups: u64,
    total_survivors: u64,
}

#[derive(Serialize)]
struct ContinuousExportSummary {
    n: u64,
    k_start: String,
    k_batch_size: usize,
    sieve_limit: usize,
    elapsed_secs: f64,
    start_batch: u64,
    next_batch: u64,
    groups_completed: u64,
    total_survivors: u64,
    output_dir: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    checkpoint_path: Option<String>,
}

fn estimate_digits(k: u64, n: u64) -> u64 {
    ((k as f64).log10() + (n as f64) * std::f64::consts::LOG10_2).floor() as u64 + 1
}

fn batch_k0(k_start: u64, k_batch_size: usize, batch_index: u64) -> Option<u64> {
    let batch_k0_u128 = k_start as u128 + 6u128 * k_batch_size as u128 * batch_index as u128;
    if batch_k0_u128 > u64::MAX as u128 {
        None
    } else {
        Some(batch_k0_u128 as u64)
    }
}

fn write_text(path: &Path, text: String) {
    fs::write(path, text).expect("failed to write export file");
}

fn path_label(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string_lossy().to_string())
}

fn append_json_line<T: Serialize>(path: &Path, value: &T) {
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .expect("failed to open append file");
    serde_json::to_writer(&mut file, value).expect("failed to serialize JSON line");
    file.write_all(b"\n").expect("failed to write newline");
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &Path) -> T {
    let text = fs::read_to_string(path).expect("failed to read JSON file");
    serde_json::from_str(&text).expect("failed to parse JSON file")
}

fn save_json<T: Serialize>(path: &Path, value: &T) {
    let text = serde_json::to_string_pretty(value).expect("failed to serialize JSON");
    fs::write(path, text).expect("failed to write JSON");
}

fn within_budget(t0: Instant, max_seconds: f64) -> bool {
    max_seconds <= 0.0 || t0.elapsed().as_secs_f64() < max_seconds
}

fn export_survivor_batches(
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    batch_start: u64,
    batch_count: u64,
    export_prefix: &Path,
) -> ExportSummary {
    let k_start = align_k_start(k_start);
    let plan = build_k_sieve_plan(n, sieve_limit);

    let k_file = export_prefix.with_extension("k.txt");
    let minus_file = export_prefix.with_extension("minus.txt");
    let plus_file = export_prefix.with_extension("plus.txt");
    let pairs_file = export_prefix.with_extension("pairs.tsv");

    let mut k_lines = String::new();
    let mut minus_lines = String::new();
    let mut plus_lines = String::new();
    let mut pairs_lines = String::from("k\tminus\tplus\n");
    let mut total_survivors = 0u64;

    for batch_index in batch_start..batch_start + batch_count {
        let batch_k0 = match batch_k0(k_start, k_batch_size, batch_index) {
            Some(value) => value,
            None => break,
        };
        let alive = sieve_k_batch(n, batch_k0, k_batch_size, &plan);
        let survivors = survivor_offsets(&alive, k_batch_size);
        total_survivors += survivors.len() as u64;

        for off in survivors {
            let k = batch_k0 + 6 * off;
            let minus_expr = format!("{}*2^{}-1", k, n);
            let plus_expr = format!("{}*2^{}+1", k, n);
            k_lines.push_str(&format!("{}\n", k));
            minus_lines.push_str(&minus_expr);
            minus_lines.push('\n');
            plus_lines.push_str(&plus_expr);
            plus_lines.push('\n');
            pairs_lines.push_str(&format!("{}\t{}\t{}\n", k, minus_expr, plus_expr));
        }
    }

    write_text(&k_file, k_lines);
    write_text(&minus_file, minus_lines);
    write_text(&plus_file, plus_lines);
    write_text(&pairs_file, pairs_lines);

    ExportSummary {
        n,
        k_start: k_start.to_string(),
        k_batch_size,
        sieve_limit,
        batch_start,
        batch_count,
        total_input_k: batch_count * k_batch_size as u64,
        total_survivors,
        export_prefix: path_label(export_prefix),
        k_file: path_label(&k_file),
        minus_file: path_label(&minus_file),
        plus_file: path_label(&plus_file),
        pairs_file: path_label(&pairs_file),
    }
}

fn validate_k_sieve(n: u64, k_start: u64, k_count: usize, sieve_limit: usize) -> ValidationResult {
    assert!(n < 63, "validation is only intended for small n");
    let k_start = align_k_start(k_start);
    let plan = build_k_sieve_plan(n, sieve_limit);
    let primes = prime_sieve(sieve_limit);
    let alive = sieve_k_batch(n, k_start, k_count, &plan);
    let two_n = 1u128 << n;
    let mut exact_survivors = 0usize;
    let mut mismatches = 0usize;

    for idx in 0..k_count {
        let k = k_start as u128 + 6u128 * idx as u128;
        let p = k * two_n - 1;
        let q = p + 2;
        let mut brute_alive = true;
        for &prime in &primes {
            if prime < 5 {
                continue;
            }
            let prime_u128 = prime as u128;
            if p % prime_u128 == 0 && p != prime_u128 {
                brute_alive = false;
                break;
            }
            if q % prime_u128 == 0 && q != prime_u128 {
                brute_alive = false;
                break;
            }
        }
        if brute_alive {
            exact_survivors += 1;
        }
        if bit_is_set(&alive, idx) != brute_alive {
            mismatches += 1;
        }
    }

    ValidationResult {
        n,
        k_start,
        k_count,
        sieve_limit,
        predicted_survivors: survivor_offsets(&alive, k_count).len(),
        exact_survivors,
        mismatch_count: mismatches,
    }
}

fn search_fixed_n(
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    max_seconds: f64,
    max_batches: Option<u64>,
) -> SearchResult {
    let k_start = align_k_start(k_start);
    let plan = build_k_sieve_plan(n, sieve_limit);
    let hl_raw = {
        let ln_n = n as f64 * std::f64::consts::LN_2;
        ln_n * ln_n / (2.0 * C2)
    };
    let survival = survival_rate(&plan);
    let estimated = estimate_digits(k_start, n);

    println!("Fixed-n k-space search");
    println!("  n = {} (~{} digits at k={})", n, estimated, k_start);
    println!("  k progression: k == 3 (mod 6)");
    println!("  sieve limit: {}", sieve_limit);
    println!("  batch size: {}", k_batch_size);
    println!("  HL raw trials/pair: {:.0}", hl_raw);
    println!("  sieve survival: {:.6}", survival);
    println!("  estimated sieved trials/pair: {:.0}", hl_raw * survival);

    let t0 = Instant::now();
    let mut batches = 0u64;
    let mut total_sieved = 0u64;
    let mut total_plus = 0u64;
    let mut total_minus = 0u64;

    while t0.elapsed().as_secs_f64() < max_seconds {
        if let Some(limit) = max_batches {
            if batches >= limit {
                break;
            }
        }

        let batch_k0 = match batch_k0(k_start, k_batch_size, batches) {
            Some(value) => value,
            None => break,
        };
        let alive = sieve_k_batch(n, batch_k0, k_batch_size, &plan);
        let survivors = survivor_offsets(&alive, k_batch_size);
        total_sieved += survivors.len() as u64;
        batches += 1;

        let stop = AtomicBool::new(false);
        let plus_counter = AtomicU64::new(0);
        let minus_counter = AtomicU64::new(0);

        let found = survivors
            .par_iter()
            .map_init(TestCtx::new, |ctx, &off| {
                if stop.load(Ordering::Relaxed) {
                    return None;
                }
                if t0.elapsed().as_secs_f64() > max_seconds {
                    stop.store(true, Ordering::Relaxed);
                    return None;
                }

                plus_counter.fetch_add(1, Ordering::Relaxed);
                let k = batch_k0 + 6 * off;
                let (minus_tested, result) = ctx.test_candidate(k, n);
                if minus_tested {
                    minus_counter.fetch_add(1, Ordering::Relaxed);
                }
                if result.is_some() {
                    stop.store(true, Ordering::Relaxed);
                    return result.map(|pair| (k, pair));
                }
                None
            })
            .find_map_any(|item| item);

        total_plus += plus_counter.load(Ordering::Relaxed);
        total_minus += minus_counter.load(Ordering::Relaxed);

        if let Some((k, pair)) = found {
            let p_str = pair.p.to_string();
            let q_str = pair.q.to_string();
            let p_len = p_str.len();
            let q_len = q_str.len();
            return SearchResult {
                found: true,
                n,
                estimated_digits: estimate_digits(k, n),
                elapsed_secs: t0.elapsed().as_secs_f64(),
                batches,
                total_sieved,
                tested_plus: total_plus,
                tested_minus: total_minus,
                k_start: k_start.to_string(),
                k_batch_size,
                sieve_limit,
                k: Some(k.to_string()),
                digits: Some(p_len.max(q_len)),
                p_head: Some(p_str[..40.min(p_len)].to_string()),
                p_tail: Some(p_str[p_len.saturating_sub(40)..].to_string()),
                q_head: Some(q_str[..40.min(q_len)].to_string()),
                q_tail: Some(q_str[q_len.saturating_sub(40)..].to_string()),
                plus_mode: Some(pair.plus_mode),
            };
        }

        if batches % 5 == 0 {
            let elapsed = t0.elapsed().as_secs_f64();
            let plus_rate = total_plus as f64 / elapsed.max(1e-9);
            println!(
                "    [{:.1}s] batches={} survivors={} plus_tests={} minus_tests={} ({:.1} +tests/s)",
                elapsed, batches, total_sieved, total_plus, total_minus, plus_rate
            );
        }
    }

    SearchResult {
        found: false,
        n,
        estimated_digits: estimated,
        elapsed_secs: t0.elapsed().as_secs_f64(),
        batches,
        total_sieved,
        tested_plus: total_plus,
        tested_minus: total_minus,
        k_start: k_start.to_string(),
        k_batch_size,
        sieve_limit,
        k: None,
        digits: None,
        p_head: None,
        p_tail: None,
        q_head: None,
        q_tail: None,
        plus_mode: None,
    }
}

fn continuous_search_fixed_n(
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    max_seconds: f64,
    max_batches: Option<u64>,
    start_batch: u64,
    stop_after_finds: Option<u64>,
    status_every_batches: u64,
    results_out: Option<&Path>,
    checkpoint_out: Option<&Path>,
) -> ContinuousSearchSummary {
    let k_start = align_k_start(k_start);
    let plan = build_k_sieve_plan(n, sieve_limit);
    let hl_raw = {
        let ln_n = n as f64 * std::f64::consts::LN_2;
        ln_n * ln_n / (2.0 * C2)
    };
    let survival = survival_rate(&plan);

    println!("Continuous fixed-n k-space search");
    println!("  n = {} (~{} digits at k={})", n, estimate_digits(k_start, n), k_start);
    println!("  start batch: {}", start_batch);
    println!("  k progression: k == 3 (mod 6)");
    println!("  sieve limit: {}", sieve_limit);
    println!("  batch size: {}", k_batch_size);
    println!("  HL raw trials/pair: {:.0}", hl_raw);
    println!("  sieve survival: {:.6}", survival);
    println!("  estimated sieved trials/pair: {:.0}", hl_raw * survival);

    let t0 = Instant::now();
    let mut batch_index = start_batch;
    let mut batches_completed = 0u64;
    let mut total_sieved = 0u64;
    let mut total_plus = 0u64;
    let mut total_minus = 0u64;
    let mut total_finds = 0u64;

    while within_budget(t0, max_seconds) {
        if let Some(limit) = max_batches {
            if batches_completed >= limit {
                break;
            }
        }
        if let Some(limit) = stop_after_finds {
            if total_finds >= limit {
                break;
            }
        }

        let batch_k0 = match batch_k0(k_start, k_batch_size, batch_index) {
            Some(value) => value,
            None => break,
        };
        let alive = sieve_k_batch(n, batch_k0, k_batch_size, &plan);
        let survivors = survivor_offsets(&alive, k_batch_size);
        total_sieved += survivors.len() as u64;

        let expired = AtomicBool::new(false);
        let plus_counter = AtomicU64::new(0);
        let minus_counter = AtomicU64::new(0);

        let found_pairs: Vec<(u64, FoundPair)> = survivors
            .par_iter()
            .map_init(TestCtx::new, |ctx, &off| {
                if expired.load(Ordering::Relaxed) {
                    return None;
                }
                if !within_budget(t0, max_seconds) {
                    expired.store(true, Ordering::Relaxed);
                    return None;
                }

                plus_counter.fetch_add(1, Ordering::Relaxed);
                let k = batch_k0 + 6 * off;
                let (minus_tested, result) = ctx.test_candidate(k, n);
                if minus_tested {
                    minus_counter.fetch_add(1, Ordering::Relaxed);
                }
                result.map(|pair| (k, pair))
            })
            .filter_map(|item| item)
            .collect();

        total_plus += plus_counter.load(Ordering::Relaxed);
        total_minus += minus_counter.load(Ordering::Relaxed);
        batches_completed += 1;
        batch_index += 1;

        if let Some(results_path) = results_out {
            for (k, pair) in &found_pairs {
                let p_str = pair.p.to_string();
                let q_str = pair.q.to_string();
                let p_len = p_str.len();
                let q_len = q_str.len();
                let record = FoundPairRecord {
                    batch_index: batch_index - 1,
                    k: k.to_string(),
                    digits: p_len.max(q_len),
                    plus_mode: pair.plus_mode.clone(),
                    p_head: p_str[..40.min(p_len)].to_string(),
                    p_tail: p_str[p_len.saturating_sub(40)..].to_string(),
                    q_head: q_str[..40.min(q_len)].to_string(),
                    q_tail: q_str[q_len.saturating_sub(40)..].to_string(),
                };
                append_json_line(results_path, &record);
            }
        }
        total_finds += found_pairs.len() as u64;

        if let Some(path) = checkpoint_out {
            let checkpoint = SearchCheckpoint {
                n,
                k_start,
                k_batch_size,
                sieve_limit,
                next_batch: batch_index,
                total_sieved,
                tested_plus: total_plus,
                tested_minus: total_minus,
                total_finds,
            };
            save_json(path, &checkpoint);
        }

        if status_every_batches > 0 && batches_completed % status_every_batches == 0 {
            let elapsed = t0.elapsed().as_secs_f64();
            let plus_rate = total_plus as f64 / elapsed.max(1e-9);
            println!(
                "    [{:.1}s] next_batch={} survivors={} plus_tests={} minus_tests={} finds={} ({:.1} +tests/s)",
                elapsed, batch_index, total_sieved, total_plus, total_minus, total_finds, plus_rate
            );
        }
    }

    ContinuousSearchSummary {
        n,
        k_start: k_start.to_string(),
        k_batch_size,
        sieve_limit,
        elapsed_secs: t0.elapsed().as_secs_f64(),
        start_batch,
        next_batch: batch_index,
        batches_completed,
        total_sieved,
        tested_plus: total_plus,
        tested_minus: total_minus,
        total_finds,
        results_path: results_out.map(path_label),
        checkpoint_path: checkpoint_out.map(path_label),
    }
}

fn continuous_export_batches(
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    start_batch: u64,
    group_batches: u64,
    max_seconds: f64,
    max_groups: Option<u64>,
    status_every_groups: u64,
    output_dir: &Path,
    checkpoint_out: Option<&Path>,
) -> ContinuousExportSummary {
    let k_start = align_k_start(k_start);
    fs::create_dir_all(output_dir).expect("failed to create export directory");

    let t0 = Instant::now();
    let mut next_batch = start_batch;
    let mut groups_completed = 0u64;
    let mut total_survivors = 0u64;

    println!("Continuous fixed-n export");
    println!("  n = {} (~{} digits at k={})", n, estimate_digits(k_start, n), k_start);
    println!("  start batch: {}", start_batch);
    println!("  group batches: {}", group_batches);
    println!("  output dir: {}", output_dir.display());

    while within_budget(t0, max_seconds) {
        if let Some(limit) = max_groups {
            if groups_completed >= limit {
                break;
            }
        }

        let prefix = output_dir.join(format!("n{}_batch{:012}", n, next_batch));
        let summary = export_survivor_batches(
            n,
            k_start,
            k_batch_size,
            sieve_limit,
            next_batch,
            group_batches,
            &prefix,
        );
        total_survivors += summary.total_survivors;
        next_batch += group_batches;
        groups_completed += 1;

        let summary_path = prefix.with_extension("summary.json");
        save_json(&summary_path, &summary);

        if let Some(path) = checkpoint_out {
            let checkpoint = ExportCheckpoint {
                n,
                k_start,
                k_batch_size,
                sieve_limit,
                next_batch,
                exported_groups: groups_completed,
                total_survivors,
            };
            save_json(path, &checkpoint);
        }

        if status_every_groups > 0 && groups_completed % status_every_groups == 0 {
            println!(
                "    [{:.1}s] next_batch={} groups={} total_survivors={}",
                t0.elapsed().as_secs_f64(),
                next_batch,
                groups_completed,
                total_survivors
            );
        }
    }

    ContinuousExportSummary {
        n,
        k_start: k_start.to_string(),
        k_batch_size,
        sieve_limit,
        elapsed_secs: t0.elapsed().as_secs_f64(),
        start_batch,
        next_batch,
        groups_completed,
        total_survivors,
        output_dir: path_label(output_dir),
        checkpoint_path: checkpoint_out.map(path_label),
    }
}

fn usage() {
    eprintln!(
        "usage: cargo run --release --bin fixed_n_k_search -- --n <n> [--k-start <u64>] [--k-batch-size <usize>] [--sieve-limit <usize>] [--max-seconds <f64>] [--max-batches <u64>] [--validate] [--sieve-only] [--preview <usize>] [--export-prefix <path>] [--export-batches <u64>] [--export-start-batch <u64>] [--continuous] [--results-out <file>] [--checkpoint-out <file>] [--resume-checkpoint <file>] [--start-batch <u64>] [--stop-after-finds <u64>] [--status-every-batches <u64>] [--continuous-export-dir <dir>] [--group-batches <u64>] [--max-groups <u64>] [--status-every-groups <u64>] [--json-out <file>]"
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let n = match arg_value(&args, "--n").and_then(|s| s.parse::<u64>().ok()) {
        Some(value) => value,
        None => {
            usage();
            return;
        }
    };

    let k_start = arg_value(&args, "--k-start")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(3);
    let k_batch_size = arg_value(&args, "--k-batch-size")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(262_144);
    let sieve_limit = arg_value(&args, "--sieve-limit")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(50_000);
    let max_seconds = arg_value(&args, "--max-seconds")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(60.0);
    let max_batches = arg_value(&args, "--max-batches").and_then(|s| s.parse::<u64>().ok());
    let preview = arg_value(&args, "--preview")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(10);
    let export_prefix = arg_value(&args, "--export-prefix").map(PathBuf::from);
    let export_batches = arg_value(&args, "--export-batches")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(1);
    let export_start_batch = arg_value(&args, "--export-start-batch")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);
    let continuous = has_flag(&args, "--continuous");
    let results_out = arg_value(&args, "--results-out").map(PathBuf::from);
    let checkpoint_out = arg_value(&args, "--checkpoint-out").map(PathBuf::from);
    let resume_checkpoint = arg_value(&args, "--resume-checkpoint").map(PathBuf::from);
    let start_batch = arg_value(&args, "--start-batch")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);
    let stop_after_finds = arg_value(&args, "--stop-after-finds").and_then(|s| s.parse::<u64>().ok());
    let status_every_batches = arg_value(&args, "--status-every-batches")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(5);
    let continuous_export_dir = arg_value(&args, "--continuous-export-dir").map(PathBuf::from);
    let group_batches = arg_value(&args, "--group-batches")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(1);
    let max_groups = arg_value(&args, "--max-groups").and_then(|s| s.parse::<u64>().ok());
    let status_every_groups = arg_value(&args, "--status-every-groups")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(1);
    let json_out = arg_value(&args, "--json-out").map(PathBuf::from);

    if has_flag(&args, "--validate") {
        let result = validate_k_sieve(n, k_start, k_batch_size, sieve_limit);
        let json = serde_json::to_string_pretty(&result).expect("failed to serialize validation result");
        println!("{}", json);
        if let Some(path) = json_out {
            fs::write(path, json).expect("failed to write JSON");
        }
        return;
    }

    if has_flag(&args, "--sieve-only") {
        let k_start = align_k_start(k_start);
        let plan = build_k_sieve_plan(n, sieve_limit);
        let alive = sieve_k_batch(n, k_start, k_batch_size, &plan);
        let survivors = survivor_offsets(&alive, k_batch_size);
        let preview_k = survivors
            .iter()
            .take(preview)
            .map(|off| (k_start + 6 * off).to_string())
            .collect();
        let result = SieveOnlyResult {
            n,
            batch_k0: k_start.to_string(),
            k_batch_size,
            sieve_limit,
            estimated_digits: estimate_digits(k_start, n),
            survival_rate: survival_rate(&plan),
            survivors: survivors.len(),
            preview_k,
        };
        let json = serde_json::to_string_pretty(&result).expect("failed to serialize sieve-only result");
        println!("{}", json);
        if let Some(path) = json_out {
            fs::write(path, json).expect("failed to write JSON");
        }
        return;
    }

    if let Some(dir) = continuous_export_dir {
        let (actual_start_batch, actual_k_start, actual_k_batch_size, actual_sieve_limit, checkpoint_path) =
            if let Some(path) = resume_checkpoint.as_deref() {
                let checkpoint: ExportCheckpoint = load_json(path);
                (
                    checkpoint.next_batch,
                    checkpoint.k_start,
                    checkpoint.k_batch_size,
                    checkpoint.sieve_limit,
                    Some(path),
                )
            } else {
                (
                    start_batch,
                    k_start,
                    k_batch_size,
                    sieve_limit,
                    checkpoint_out.as_deref(),
                )
            };

        let result = continuous_export_batches(
            n,
            actual_k_start,
            actual_k_batch_size,
            actual_sieve_limit,
            actual_start_batch,
            group_batches,
            max_seconds,
            max_groups,
            status_every_groups,
            &dir,
            checkpoint_path,
        );
        let json = serde_json::to_string_pretty(&result).expect("failed to serialize continuous export summary");
        println!("{}", json);
        if let Some(path) = json_out {
            fs::write(path, json).expect("failed to write JSON");
        }
        return;
    }

    if let Some(prefix) = export_prefix {
        let result = export_survivor_batches(
            n,
            k_start,
            k_batch_size,
            sieve_limit,
            export_start_batch,
            export_batches,
            &prefix,
        );
        let json = serde_json::to_string_pretty(&result).expect("failed to serialize export summary");
        println!("{}", json);
        if let Some(path) = json_out {
            fs::write(path, json).expect("failed to write JSON");
        }
        return;
    }

    if continuous {
        let (
            actual_start_batch,
            actual_k_start,
            actual_k_batch_size,
            actual_sieve_limit,
            checkpoint_path,
        ) = if let Some(path) = resume_checkpoint.as_deref() {
            let checkpoint: SearchCheckpoint = load_json(path);
            (
                checkpoint.next_batch,
                checkpoint.k_start,
                checkpoint.k_batch_size,
                checkpoint.sieve_limit,
                Some(path),
            )
        } else {
            (
                start_batch,
                k_start,
                k_batch_size,
                sieve_limit,
                checkpoint_out.as_deref(),
            )
        };

        let result = continuous_search_fixed_n(
            n,
            actual_k_start,
            actual_k_batch_size,
            actual_sieve_limit,
            max_seconds,
            max_batches,
            actual_start_batch,
            stop_after_finds,
            status_every_batches,
            results_out.as_deref(),
            checkpoint_path,
        );
        let json = serde_json::to_string_pretty(&result).expect("failed to serialize continuous search summary");
        println!("{}", json);
        if let Some(path) = json_out {
            fs::write(path, json).expect("failed to write JSON");
        }
        return;
    }

    let result = search_fixed_n(n, k_start, k_batch_size, sieve_limit, max_seconds, max_batches);
    let json = serde_json::to_string_pretty(&result).expect("failed to serialize search result");
    println!("{}", json);
    if let Some(path) = json_out {
        fs::write(path, json).expect("failed to write JSON");
    }
}
