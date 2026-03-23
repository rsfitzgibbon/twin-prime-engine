use rayon::prelude::*;
use rug::integer::IsPrime;
use rug::Assign;
use rug::Integer;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn arg_value(args: &[String], flag: &str) -> Option<String> {
    args.windows(2)
        .find(|w| w[0] == flag)
        .map(|w| w[1].clone())
}

fn arg_values(args: &[String], flag: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == flag {
            i += 1;
            while i < args.len() && !args[i].starts_with("--") {
                out.push(args[i].clone());
                i += 1;
            }
            continue;
        }
        i += 1;
    }
    out
}

fn sanitize_filename(s: &str) -> String {
    s.chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' || ch == '.' {
                ch
            } else {
                '_'
            }
        })
        .collect()
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

fn base_sieve(alive: &mut [u64], batch_size: usize, m_start: &Integer, sieve: &SieveData) {
    for (idx, &p) in sieve.primes_small.iter().enumerate() {
        let p_u64 = p as u64;
        let inv6 = sieve.inv6_small[idx] as u64;
        let m_mod_p = m_start.mod_u(p) as u64;
        let p_us = p as usize;

        let r1 = if inv6 >= m_mod_p {
            (inv6 - m_mod_p) as usize
        } else {
            (p_u64 + inv6 - m_mod_p) as usize
        };
        let mut j = r1;
        while j < batch_size {
            alive[j >> 6] &= !(1u64 << (j & 63));
            j += p_us;
        }

        let complement = p_u64 - inv6;
        let r2 = if complement >= m_mod_p {
            (complement - m_mod_p) as usize
        } else {
            (p_u64 + complement - m_mod_p) as usize
        };
        let mut j = r2;
        while j < batch_size {
            alive[j >> 6] &= !(1u64 << (j & 63));
            j += p_us;
        }
    }
}

fn extended_sieve(alive: &mut [u64], batch_size: usize, m_start: &Integer, sieve: &SieveData) {
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
                alive[j >> 6] &= !(1u64 << (j & 63));
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
                alive[j >> 6] &= !(1u64 << (j & 63));
                j += p_us;
            }
        }
    }
}

fn deep_sieve(alive: &mut [u64], batch_size: usize, m_start: &Integer, sieve: &SieveData) {
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
                alive[j >> 6] &= !(1u64 << (j & 63));
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
                alive[j >> 6] &= !(1u64 << (j & 63));
                j += p_us;
            }
        }
    }
}

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

        for _ in 0..r.saturating_sub(1) {
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

    fn test_candidate(&mut self, off: u64, m_start: &Integer) -> Option<(Integer, Integer)> {
        self.p1.assign(m_start);
        self.p1 += off;
        self.p1 *= 6;
        self.p1 -= 1;

        if !self.sprp.is_sprp(&self.p1, 2) {
            return None;
        }

        self.p2.assign(&self.p1);
        self.p2 += 2;
        if !self.sprp.is_sprp(&self.p2, 2) {
            return None;
        }

        for &base in &[3u32, 5, 7] {
            if !self.sprp.is_sprp(&self.p1, base) {
                return None;
            }
        }
        for &base in &[3u32, 5, 7] {
            if !self.sprp.is_sprp(&self.p2, base) {
                return None;
            }
        }

        if self.p1.is_probably_prime(0) != IsPrime::No
            && self.p2.is_probably_prime(0) != IsPrime::No
        {
            return Some((self.p1.clone(), self.p2.clone()));
        }
        None
    }
}

fn collect_survivor_offsets(alive: &[u64], batch_size: usize) -> Vec<u64> {
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
            if off >= batch_size as u64 {
                break;
            }
            out.push(off);
            w &= w - 1;
        }
    }
    out
}

#[derive(Debug, Clone, Copy)]
enum SieveDepth {
    Base,
    Extended,
    Deep,
}

impl SieveDepth {
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "base" => Ok(SieveDepth::Base),
            "extended" => Ok(SieveDepth::Extended),
            "deep" => Ok(SieveDepth::Deep),
            _ => Err(format!("unsupported sieve_depth: {}", s)),
        }
    }

    fn as_str(&self) -> &'static str {
        match self {
            SieveDepth::Base => "base",
            SieveDepth::Extended => "extended",
            SieveDepth::Deep => "deep",
        }
    }
}

#[derive(Deserialize)]
struct WorkUnitRaw {
    version: u32,
    job_id: String,
    target_digits: u32,
    m_start: Value,
    batch_size: usize,
    sieve_depth: String,
    sieve_limit: usize,
}

struct WorkUnit {
    job_id: String,
    target_digits: u32,
    m_start: Integer,
    batch_size: usize,
    sieve_depth: SieveDepth,
    sieve_limit: usize,
}

#[derive(Serialize)]
struct FoundPair {
    m: String,
    p: String,
    q: String,
    digits: usize,
}

#[derive(Serialize)]
struct WorkerResult {
    job_id: String,
    status: String,
    target_digits: u32,
    m_start: String,
    batch_size: usize,
    sieve_depth: String,
    worker_setup_seconds: f64,
    elapsed_seconds: f64,
    total_raw: usize,
    survivors: usize,
    tested: usize,
    found_count: usize,
    found_pairs: Vec<FoundPair>,
}

#[derive(Serialize)]
struct SummaryRow {
    job_id: String,
    result_path: String,
    survivors: usize,
    tested: usize,
    found_count: usize,
    elapsed_seconds: f64,
}

#[derive(Serialize)]
struct RunSummary {
    worker: String,
    work_units: usize,
    sieve_limit: usize,
    prime_cache_mb: f64,
    worker_setup_seconds: f64,
    results: Vec<SummaryRow>,
}

fn parse_integer_value(value: &Value) -> Result<Integer, String> {
    match value {
        Value::String(s) => Integer::from_str_radix(s, 10).map_err(|e| format!("invalid m_start string: {}", e)),
        Value::Number(n) => {
            if let Some(u) = n.as_u64() {
                Ok(Integer::from(u))
            } else {
                Err("m_start number must fit in u64 or be encoded as a string".to_string())
            }
        }
        _ => Err("m_start must be a decimal string or integer".to_string()),
    }
}

fn load_work_unit(path: &Path) -> Result<WorkUnit, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("failed to read {}: {}", path.display(), e))?;
    let raw: WorkUnitRaw =
        serde_json::from_str(&text).map_err(|e| format!("failed to parse {}: {}", path.display(), e))?;
    if raw.version == 0 {
        return Err(format!("{} has unsupported version 0", path.display()));
    }
    let sieve_depth = SieveDepth::from_str(&raw.sieve_depth)?;
    let m_start = parse_integer_value(&raw.m_start)?;
    Ok(WorkUnit {
        job_id: raw.job_id,
        target_digits: raw.target_digits,
        m_start,
        batch_size: raw.batch_size,
        sieve_depth,
        sieve_limit: raw.sieve_limit,
    })
}

fn run_assigned_batch(unit: &WorkUnit, sieve: &SieveData, worker_setup_seconds: f64) -> WorkerResult {
    let t0 = Instant::now();
    let words = (unit.batch_size + 63) / 64;
    let mut alive = vec![u64::MAX; words];

    base_sieve(&mut alive, unit.batch_size, &unit.m_start, sieve);
    if matches!(unit.sieve_depth, SieveDepth::Extended | SieveDepth::Deep) {
        extended_sieve(&mut alive, unit.batch_size, &unit.m_start, sieve);
    }
    if matches!(unit.sieve_depth, SieveDepth::Deep) {
        deep_sieve(&mut alive, unit.batch_size, &unit.m_start, sieve);
    }

    let survivor_offsets = collect_survivor_offsets(&alive, unit.batch_size);
    let m_start = unit.m_start.clone();
    let mut found_pairs: Vec<FoundPair> = survivor_offsets
        .par_iter()
        .map_init(TestCtx::new, |ctx, &off| {
            ctx.test_candidate(off, &m_start).map(|(p, q)| {
                let mut m_val = Integer::from(&m_start);
                m_val += off;
                let p_str = p.to_string();
                let q_str = q.to_string();
                FoundPair {
                    m: m_val.to_string(),
                    digits: p_str.len(),
                    p: p_str,
                    q: q_str,
                }
            })
        })
        .filter_map(|x| x)
        .collect();
    found_pairs.sort_by(|a, b| a.m.cmp(&b.m));

    WorkerResult {
        job_id: unit.job_id.clone(),
        status: "ok".to_string(),
        target_digits: unit.target_digits,
        m_start: unit.m_start.to_string(),
        batch_size: unit.batch_size,
        sieve_depth: unit.sieve_depth.as_str().to_string(),
        worker_setup_seconds,
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        total_raw: unit.batch_size,
        survivors: survivor_offsets.len(),
        tested: survivor_offsets.len(),
        found_count: found_pairs.len(),
        found_pairs,
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let work_units = arg_values(&args, "--work-unit");
    let result_dir = arg_value(&args, "--result-dir");
    let summary_out = arg_value(&args, "--summary-out");

    if work_units.is_empty() || result_dir.is_none() {
        eprintln!("usage: volunteer_worker --work-unit <unit1.json> [unit2.json ...] --result-dir <dir> [--summary-out <file>]");
        std::process::exit(2);
    }

    let result_dir = PathBuf::from(result_dir.unwrap());
    let loaded_units: Vec<WorkUnit> = work_units
        .iter()
        .map(|path| load_work_unit(Path::new(path)).unwrap_or_else(|e| panic!("{}", e)))
        .collect();
    let max_sieve_limit = loaded_units
        .iter()
        .map(|unit| unit.sieve_limit)
        .max()
        .unwrap_or(1_000_000);

    let setup_t0 = Instant::now();
    let sieve = SieveData::build(max_sieve_limit);
    let worker_setup_seconds = setup_t0.elapsed().as_secs_f64();
    let cache_mb = sieve.cache_bytes() as f64 / (1024.0 * 1024.0);

    println!(
        "Resident volunteer worker: {} jobs | sieve_limit={} | cache {:.1} MB | setup {:.2}s",
        loaded_units.len(),
        max_sieve_limit,
        cache_mb,
        worker_setup_seconds
    );

    fs::create_dir_all(&result_dir).expect("failed to create result dir");
    let mut summary_rows = Vec::new();

    for unit in &loaded_units {
        let result = run_assigned_batch(unit, &sieve, worker_setup_seconds);
        let file_name = format!("{}.result.json", sanitize_filename(&unit.job_id));
        let result_path = result_dir.join(file_name);
        let json = serde_json::to_string_pretty(&result).expect("failed to serialize result");
        fs::write(&result_path, json).expect("failed to write result");
        println!(
            "  {} -> survivors={} tested={} found={} [{}]",
            unit.job_id,
            result.survivors,
            result.tested,
            result.found_count,
            result_path.display()
        );
        summary_rows.push(SummaryRow {
            job_id: unit.job_id.clone(),
            result_path: format!("results/{}.result.json", sanitize_filename(&unit.job_id)),
            survivors: result.survivors,
            tested: result.tested,
            found_count: result.found_count,
            elapsed_seconds: result.elapsed_seconds,
        });
    }

    if let Some(path) = summary_out {
        let summary = RunSummary {
            worker: "resident_rust_volunteer_worker".to_string(),
            work_units: loaded_units.len(),
            sieve_limit: max_sieve_limit,
            prime_cache_mb: cache_mb,
            worker_setup_seconds,
            results: summary_rows,
        };
        let json = serde_json::to_string_pretty(&summary).expect("failed to serialize summary");
        fs::write(path, json).expect("failed to write summary");
    }
}
