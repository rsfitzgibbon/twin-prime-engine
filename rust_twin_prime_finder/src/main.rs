use std::cmp::max;
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::time::Instant;

#[derive(Clone, Debug)]
struct Candidate {
    n: usize,
    pair_right: usize,
    mode: String,
    hard_cutoff: usize,
    soft_to: Option<usize>,
    score_z: Option<usize>,
    extra_blocker: Option<usize>,
    rank_value: Option<usize>,
    sieve_score: Option<f64>,
    rank: usize,
}

#[derive(Clone, Debug)]
struct FinderResult {
    limit_n: usize,
    mode: String,
    hard_cutoff: usize,
    soft_to: Option<usize>,
    score_z: Option<usize>,
    candidate_count: usize,
    ranking: String,
    candidates: Vec<Candidate>,
}

#[derive(Clone, Debug)]
struct BenchmarkRow {
    limit_n: usize,
    mode: String,
    hard_cutoff: usize,
    soft_to: Option<usize>,
    score_z: Option<usize>,
    ranking: String,
    candidate_count: usize,
    true_twins: usize,
    total_twins: usize,
    precision: f64,
    recall: f64,
    f1: f64,
    elimination_rate: f64,
}

fn pct(num: usize, den: usize) -> f64 {
    if den == 0 {
        0.0
    } else {
        100.0 * (num as f64) / (den as f64)
    }
}

fn parse_usize(value: &str, flag: &str) -> usize {
    value.parse::<usize>().unwrap_or_else(|_| panic!("invalid integer for {}: {}", flag, value))
}

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
            break;
        }
        i += 1;
    }
    out
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|arg| arg == flag)
}

fn integer_sqrt(n: usize) -> usize {
    (n as f64).sqrt() as usize
}

fn sieve_bool(limit: usize) -> Vec<bool> {
    if limit < 2 {
        return vec![false; limit + 1];
    }
    let mut sieve = vec![true; limit + 1];
    sieve[0] = false;
    sieve[1] = false;
    let mut p = 2usize;
    while p * p <= limit {
        if sieve[p] {
            let mut k = p * p;
            while k <= limit {
                sieve[k] = false;
                k += p;
            }
        }
        p += 1;
    }
    sieve
}

fn prime_list_from_sieve(sieve: &[bool], limit: usize) -> Vec<usize> {
    (2..=limit).filter(|&p| sieve[p]).collect()
}

fn largest_prime_leq(sieve: &[bool], limit: usize) -> Option<usize> {
    (2..=limit).rev().find(|&p| sieve[p])
}

fn egcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x, y) = egcd(b, a % b);
        (g, y, x - (a / b) * y)
    }
}

fn mod_inverse(a: usize, m: usize) -> usize {
    let (g, x, _) = egcd(a as i64, m as i64);
    if g != 1 {
        panic!("no modular inverse for {} mod {}", a, m);
    }
    ((x % m as i64 + m as i64) % m as i64) as usize
}

fn resolve_cutoff(limit_n: usize, mode: &str, custom_cutoff: Option<usize>) -> (usize, String) {
    if let Some(cutoff) = custom_cutoff {
        return (cutoff, "custom".to_string());
    }
    match mode {
        "fast" => (211, "fast".to_string()),
        "high_precision" => (503, "high_precision".to_string()),
        "exact_range" => {
            let root = integer_sqrt(limit_n + 2);
            let sieve = sieve_bool(root);
            (
                largest_prime_leq(&sieve, root).expect("no prime found for exact_range"),
                "exact_range".to_string(),
            )
        }
        _ => panic!("unknown mode: {}", mode),
    }
}

fn resolve_score_z(hard_cutoff: usize, score_z_arg: Option<usize>, limit_n: usize) -> usize {
    if let Some(score_z) = score_z_arg {
        return score_z;
    }
    match hard_cutoff {
        211 => 503,
        503 => 997,
        _ => {
            let root = integer_sqrt(limit_n + 2);
            let candidate_z = hard_cutoff.saturating_mul(2);
            if candidate_z > hard_cutoff {
                candidate_z.min(root)
            } else {
                hard_cutoff + 100
            }
        }
    }
}

fn max_candidate_m(limit_n: usize) -> usize {
    if limit_n <= 1 {
        0
    } else {
        (limit_n - 1) / 6
    }
}

fn build_survivor_mask(max_m: usize, hard_cutoff: usize, primes: &[usize]) -> Vec<bool> {
    let mut alive = vec![true; max_m];
    for &p in primes {
        if p < 5 || p > hard_cutoff {
            continue;
        }
        let inv6 = mod_inverse(6, p);
        let left_exception = if p % 6 == 5 { Some((p + 1) / 6) } else { None };
        let right_exception = if p % 6 == 1 { Some((p - 1) / 6) } else { None };

        let mut m = inv6;
        while m <= max_m {
            if Some(m) != left_exception {
                alive[m - 1] = false;
            }
            m += p;
        }

        let r2 = p - inv6;
        let mut m = r2;
        while m <= max_m {
            if Some(m) != right_exception {
                alive[m - 1] = false;
            }
            m += p;
        }
    }
    alive
}

fn first_pair_blocker(n: usize, primes: &[usize], low_exclusive: usize, high_inclusive: usize) -> Option<usize> {
    let q = n + 2;
    for &p in primes {
        if p <= low_exclusive {
            continue;
        }
        if p > high_inclusive {
            return None;
        }
        if n % p == 0 && n != p {
            return Some(p);
        }
        if q % p == 0 && q != p {
            return Some(p);
        }
    }
    None
}

fn score_rec(factors: &[f64], log_r: f64, idx: usize, logd: f64, sign: f64, total: &mut f64) {
    for j in idx..factors.len() {
        let nlogd = logd + factors[j];
        if nlogd > log_r + 1e-15 {
            continue;
        }
        *total += sign * -1.0 * (1.0 - nlogd / log_r);
        score_rec(factors, log_r, j + 1, nlogd, sign * -1.0, total);
    }
}

fn sieve_score(n: usize, primes: &[usize], score_z: usize) -> f64 {
    let q = n + 2;
    let log_r = 2.0 * (score_z as f64).ln();
    let mut factors = Vec::new();
    for &p in primes {
        if p > score_z {
            break;
        }
        let x_blocked = n % p == 0 && n != p;
        let q_blocked = q % p == 0 && q != p;
        if x_blocked || q_blocked {
            factors.push((p as f64).ln());
        }
    }

    let mut total = 1.0;
    score_rec(&factors, log_r, 0, 0.0, 1.0, &mut total);
    (total * total).ln_1p()
}

fn ranking_label(soft_to: Option<usize>, hard_cutoff: usize, score_z: Option<usize>) -> String {
    if let Some(score_z) = score_z {
        return format!("sieve_score(score_z={})", score_z);
    }
    if let Some(soft_to) = soft_to {
        if soft_to > hard_cutoff {
            return format!("soft_rerank(soft_to={})", soft_to);
        }
    }
    "natural_order".to_string()
}

fn find_twin_candidates(
    limit_n: usize,
    mode: &str,
    soft_to: Option<usize>,
    top_k: Option<usize>,
    custom_cutoff: Option<usize>,
    score: bool,
    score_z_arg: Option<usize>,
) -> FinderResult {
    let (hard_cutoff, resolved_mode) = resolve_cutoff(limit_n, mode, custom_cutoff);
    let actual_score_z = if score {
        Some(resolve_score_z(hard_cutoff, score_z_arg, limit_n))
    } else {
        None
    };
    let max_m = max_candidate_m(limit_n);
    let max_prime_needed = max(
        max(hard_cutoff, soft_to.unwrap_or(hard_cutoff)),
        actual_score_z.unwrap_or(hard_cutoff),
    );
    let sieve = sieve_bool(max_prime_needed);
    let primes = prime_list_from_sieve(&sieve, max_prime_needed);
    let alive = build_survivor_mask(max_m, hard_cutoff, &primes);

    let mut candidates = Vec::new();
    for (idx, &is_alive) in alive.iter().enumerate() {
        if !is_alive {
            continue;
        }
        let m = idx + 1;
        let n = 6 * m - 1;
        let extra_blocker = if let Some(soft_limit) = soft_to {
            if soft_limit > hard_cutoff {
                first_pair_blocker(n, &primes, hard_cutoff, soft_limit)
            } else {
                None
            }
        } else {
            None
        };
        let rank_value = if let Some(soft_limit) = soft_to {
            if soft_limit > hard_cutoff {
                Some(extra_blocker.unwrap_or(soft_limit + 1))
            } else {
                None
            }
        } else {
            None
        };
        let sieve_score_value = actual_score_z.map(|score_z| sieve_score(n, &primes, score_z));
        candidates.push(Candidate {
            n,
            pair_right: n + 2,
            mode: resolved_mode.clone(),
            hard_cutoff,
            soft_to,
            score_z: actual_score_z,
            extra_blocker,
            rank_value,
            sieve_score: sieve_score_value,
            rank: 0,
        });
    }

    if actual_score_z.is_some() {
        candidates.sort_by(|a, b| {
            b.sieve_score
                .unwrap_or(f64::NEG_INFINITY)
                .partial_cmp(&a.sieve_score.unwrap_or(f64::NEG_INFINITY))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.n.cmp(&b.n))
        });
    } else if let Some(soft_limit) = soft_to {
        if soft_limit > hard_cutoff {
            candidates.sort_by(|a, b| {
                let av = a.rank_value.unwrap_or(soft_limit + 1);
                let bv = b.rank_value.unwrap_or(soft_limit + 1);
                bv.cmp(&av).then(a.n.cmp(&b.n))
            });
        }
    } else {
        candidates.sort_by_key(|c| c.n);
    }

    if let Some(k) = top_k {
        candidates.truncate(k);
    }

    for (idx, candidate) in candidates.iter_mut().enumerate() {
        candidate.rank = idx + 1;
    }

    FinderResult {
        limit_n,
        mode: resolved_mode,
        hard_cutoff,
        soft_to,
        score_z: actual_score_z,
        candidate_count: candidates.len(),
        ranking: ranking_label(soft_to, hard_cutoff, actual_score_z),
        candidates,
    }
}

fn evaluate_result(result: &FinderResult) -> BenchmarkRow {
    let limit_n = result.limit_n;
    let sieve = sieve_bool(limit_n + 2);
    let true_twins = result
        .candidates
        .iter()
        .filter(|c| sieve[c.n] && sieve[c.pair_right])
        .count();
    let total_twins = (5..=limit_n.saturating_sub(2))
        .step_by(6)
        .filter(|&n| sieve[n] && sieve[n + 2])
        .count();
    let precision = if result.candidate_count == 0 {
        0.0
    } else {
        true_twins as f64 / result.candidate_count as f64
    };
    let recall = if total_twins == 0 {
        0.0
    } else {
        true_twins as f64 / total_twins as f64
    };
    let f1 = if precision + recall == 0.0 {
        0.0
    } else {
        2.0 * precision * recall / (precision + recall)
    };
    let total_corridor = max_candidate_m(limit_n);
    let elimination_rate = if total_corridor == 0 {
        0.0
    } else {
        1.0 - (result.candidate_count as f64 / total_corridor as f64)
    };
    BenchmarkRow {
        limit_n,
        mode: result.mode.clone(),
        hard_cutoff: result.hard_cutoff,
        soft_to: result.soft_to,
        score_z: result.score_z,
        ranking: result.ranking.clone(),
        candidate_count: result.candidate_count,
        true_twins,
        total_twins,
        precision,
        recall,
        f1,
        elimination_rate,
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn result_to_json(result: &FinderResult) -> String {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("  \"limit_n\": {},\n", result.limit_n));
    out.push_str(&format!("  \"mode\": \"{}\",\n", escape_json(&result.mode)));
    out.push_str(&format!("  \"hard_cutoff\": {},\n", result.hard_cutoff));
    if let Some(soft_to) = result.soft_to {
        out.push_str(&format!("  \"soft_to\": {},\n", soft_to));
    } else {
        out.push_str("  \"soft_to\": null,\n");
    }
    if let Some(score_z) = result.score_z {
        out.push_str(&format!("  \"score_z\": {},\n", score_z));
    } else {
        out.push_str("  \"score_z\": null,\n");
    }
    out.push_str(&format!("  \"candidate_count\": {},\n", result.candidate_count));
    out.push_str("  \"engine\": \"rust\",\n");
    out.push_str(&format!("  \"ranking\": \"{}\",\n", escape_json(&result.ranking)));
    out.push_str("  \"candidates\": [\n");
    for (idx, c) in result.candidates.iter().enumerate() {
        out.push_str("    {");
        out.push_str(&format!("\"rank\": {}, ", c.rank));
        out.push_str(&format!("\"n\": {}, ", c.n));
        out.push_str(&format!("\"pair\": [{}, {}], ", c.n, c.pair_right));
        out.push_str(&format!("\"mode\": \"{}\", ", escape_json(&c.mode)));
        out.push_str(&format!("\"hard_cutoff\": {}, ", c.hard_cutoff));
        if let Some(soft_to) = c.soft_to {
            out.push_str(&format!("\"soft_to\": {}, ", soft_to));
        } else {
            out.push_str("\"soft_to\": null, ");
        }
        if let Some(score_z) = c.score_z {
            out.push_str(&format!("\"score_z\": {}, ", score_z));
        } else {
            out.push_str("\"score_z\": null, ");
        }
        if let Some(blocker) = c.extra_blocker {
            out.push_str(&format!("\"extra_blocker\": {}, ", blocker));
        } else {
            out.push_str("\"extra_blocker\": null, ");
        }
        if let Some(rank_value) = c.rank_value {
            out.push_str(&format!("\"rank_value\": {}, ", rank_value));
        } else {
            out.push_str("\"rank_value\": null, ");
        }
        if let Some(sieve_score) = c.sieve_score {
            out.push_str(&format!("\"sieve_score\": {:.12}", sieve_score));
        } else {
            out.push_str("\"sieve_score\": null");
        }
        out.push('}');
        if idx + 1 != result.candidates.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ]\n");
    out.push_str("}\n");
    out
}

fn benchmark_rows_to_json(rows: &[BenchmarkRow]) -> String {
    let mut out = String::new();
    out.push_str("[\n");
    for (idx, row) in rows.iter().enumerate() {
        out.push_str("  {\n");
        out.push_str(&format!("    \"limit_n\": {},\n", row.limit_n));
        out.push_str(&format!("    \"mode\": \"{}\",\n", escape_json(&row.mode)));
        out.push_str(&format!("    \"hard_cutoff\": {},\n", row.hard_cutoff));
        if let Some(soft_to) = row.soft_to {
            out.push_str(&format!("    \"soft_to\": {},\n", soft_to));
        } else {
            out.push_str("    \"soft_to\": null,\n");
        }
        if let Some(score_z) = row.score_z {
            out.push_str(&format!("    \"score_z\": {},\n", score_z));
        } else {
            out.push_str("    \"score_z\": null,\n");
        }
        out.push_str(&format!("    \"ranking\": \"{}\",\n", escape_json(&row.ranking)));
        out.push_str(&format!("    \"candidate_count\": {},\n", row.candidate_count));
        out.push_str(&format!("    \"true_twins\": {},\n", row.true_twins));
        out.push_str(&format!("    \"total_twins\": {},\n", row.total_twins));
        out.push_str(&format!("    \"precision\": {:.12},\n", row.precision));
        out.push_str(&format!("    \"recall\": {:.12},\n", row.recall));
        out.push_str(&format!("    \"f1\": {:.12},\n", row.f1));
        out.push_str(&format!("    \"elimination_rate\": {:.12}\n", row.elimination_rate));
        out.push_str("  }");
        if idx + 1 != rows.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("]\n");
    out
}

fn write_result_csv(result: &FinderResult, path: &str) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "rank,n,left,right,mode,hard_cutoff,ranking,score_z,sieve_score,extra_blocker,rank_value")?;
    for c in &result.candidates {
        let score_z = c.score_z.map(|v| v.to_string()).unwrap_or_default();
        let sieve_score = c
            .sieve_score
            .map(|v| format!("{:.6}", v))
            .unwrap_or_default();
        let extra = c.extra_blocker.map(|v| v.to_string()).unwrap_or_default();
        let rank_value = c.rank_value.map(|v| v.to_string()).unwrap_or_default();
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{},{},{}",
            c.rank, c.n, c.n, c.pair_right, c.mode, c.hard_cutoff, result.ranking, score_z, sieve_score, extra, rank_value
        )?;
    }
    Ok(())
}

fn write_benchmark_csv(rows: &[BenchmarkRow], path: &str) -> io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "limit_n,mode,hard_cutoff,soft_to,score_z,ranking,candidate_count,true_twins,total_twins,precision,recall,f1,elimination_rate"
    )?;
    for row in rows {
        let soft_to = row.soft_to.map(|v| v.to_string()).unwrap_or_default();
        let score_z = row.score_z.map(|v| v.to_string()).unwrap_or_default();
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{},{:.12},{:.12},{:.12},{:.12}",
            row.limit_n,
            row.mode,
            row.hard_cutoff,
            soft_to,
            score_z,
            row.ranking,
            row.candidate_count,
            row.true_twins,
            row.total_twins,
            row.precision,
            row.recall,
            row.f1,
            row.elimination_rate
        )?;
    }
    Ok(())
}

fn print_result_summary(result: &FinderResult, preview_count: usize) {
    println!("Twin Prime Finder Rust  N = {}", format_number(result.limit_n));
    println!("Mode: {}  hard_cutoff={}", result.mode, result.hard_cutoff);
    println!("Ranking: {}", result.ranking);
    if let Some(score_z) = result.score_z {
        println!("Sieve scoring: yes  score_z={}", score_z);
    }
    println!("Predicted candidate pairs: {}", format_number(result.candidate_count));
    let preview = result.candidates.iter().take(preview_count);
    println!("\nTop {} candidates:", preview_count.min(result.candidates.len()));
    for c in preview {
        let mut line = format!("  pair=({}, {})", c.n, c.pair_right);
        if let Some(sieve_score) = c.sieve_score {
            line.push_str(&format!(", score={:.4}", sieve_score));
        }
        if let Some(blocker) = c.extra_blocker {
            line.push_str(&format!(", extra_blocker={}", blocker));
        }
        println!("{}", line);
    }
}

fn print_benchmark(rows: &[BenchmarkRow]) {
    println!("Twin Prime Finder Rust Benchmark");
    println!(
        "\n  {:>10} | {:>14} | {:>6} | {:>11} | {:>8} | {:>10} | {:>8} | {:>8}",
        "N", "mode", "hard", "candidates", "twins", "precision", "recall", "F1"
    );
    println!(
        "  -----------+----------------+--------+-------------+----------+------------+----------+---------"
    );
    for row in rows {
        println!(
            "  {:>10} | {:>14} | {:>6} | {:>11} | {:>8} | {:>9.4}% | {:>7.2}% | {:>7.2}%",
            format_number(row.limit_n),
            row.mode,
            row.hard_cutoff,
            format_number(row.candidate_count),
            format_number(row.true_twins),
            pct(row.true_twins, row.candidate_count),
            100.0 * row.recall,
            100.0 * row.f1
        );
    }
}

fn format_number(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (idx, ch) in s.chars().rev().enumerate() {
        if idx > 0 && idx % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let benchmark = has_flag(&args, "--benchmark");
    let limit_n = arg_value(&args, "--limit").map(|v| parse_usize(&v, "--limit")).unwrap_or(1_000_000);
    let mode = arg_value(&args, "--mode").unwrap_or_else(|| "high_precision".to_string());
    let preview = arg_value(&args, "--preview").map(|v| parse_usize(&v, "--preview")).unwrap_or(20);
    let top_k = arg_value(&args, "--top-k").map(|v| parse_usize(&v, "--top-k"));
    let soft_to = arg_value(&args, "--soft-to").map(|v| parse_usize(&v, "--soft-to"));
    let score = has_flag(&args, "--score");
    let score_z = arg_value(&args, "--score-z").map(|v| parse_usize(&v, "--score-z"));
    let cutoff = arg_value(&args, "--cutoff").map(|v| parse_usize(&v, "--cutoff"));
    let json_out = arg_value(&args, "--json-out");
    let csv_out = arg_value(&args, "--csv-out");
    let benchmark_json_out = arg_value(&args, "--benchmark-json-out");
    let benchmark_csv_out = arg_value(&args, "--benchmark-csv-out");
    let benchmark_limits: Vec<usize> = {
        let values = arg_values(&args, "--benchmark-limits");
        if values.is_empty() {
            vec![100_000, 500_000, 1_000_000]
        } else {
            values.iter().map(|v| parse_usize(v, "--benchmark-limits")).collect()
        }
    };
    let benchmark_modes: Vec<String> = {
        let values = arg_values(&args, "--benchmark-modes");
        if values.is_empty() {
            vec!["fast".to_string(), "high_precision".to_string(), "exact_range".to_string()]
        } else {
            values
        }
    };

    let t0 = Instant::now();
    if benchmark {
        let mut rows = Vec::new();
        for &limit in &benchmark_limits {
            for mode_name in &benchmark_modes {
                let result = find_twin_candidates(limit, mode_name, soft_to, top_k, cutoff, score, score_z);
                rows.push(evaluate_result(&result));
            }
        }
        print_benchmark(&rows);
        if let Some(path) = benchmark_csv_out {
            write_benchmark_csv(&rows, &path).expect("failed to write benchmark csv");
            println!("Saved benchmark CSV to {}", path);
        }
        if let Some(path) = benchmark_json_out {
            std::fs::write(&path, benchmark_rows_to_json(&rows)).expect("failed to write benchmark json");
            println!("Saved benchmark JSON to {}", path);
        }
    } else {
        let result = find_twin_candidates(limit_n, &mode, soft_to, top_k, cutoff, score, score_z);
        print_result_summary(&result, preview);
        if let Some(path) = csv_out {
            write_result_csv(&result, &path).expect("failed to write csv");
            println!("Saved CSV to {} ({}) rows", path, format_number(result.candidate_count));
        }
        if let Some(path) = json_out {
            std::fs::write(&path, result_to_json(&result)).expect("failed to write json");
            println!("Saved JSON to {}", path);
        }
    }
    eprintln!("elapsed_ms={:.3}", t0.elapsed().as_secs_f64() * 1000.0);
}
