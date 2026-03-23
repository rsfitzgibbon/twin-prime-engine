use serde::{Deserialize, Serialize};
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use twin_prime_engine::backend::compact_export::CompactExportBackend;
use twin_prime_engine::backend::export_files::ExportFilesBackend;
use twin_prime_engine::backend::local_gmp::LocalGmpBackend;
use twin_prime_engine::backend::{BackendBatchResult, FixedNBackend};
use twin_prime_engine::fixed_n::{
    align_k_start, build_k_sieve_plan, build_k_sieve_plan_range, filter_survivors_with_plan,
    generate_survivor_batch_with_plan, survival_rate,
};

fn arg_value(args: &[String], flag: &str) -> Option<String> {
    args.windows(2)
        .find(|w| w[0] == flag)
        .map(|w| w[1].clone())
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

fn save_json<T: Serialize>(path: &Path, value: &T) {
    let text = serde_json::to_string_pretty(value).expect("failed to serialize JSON");
    fs::write(path, text).expect("failed to write JSON");
}

fn load_json<T: for<'de> Deserialize<'de>>(path: &Path) -> T {
    let text = fs::read_to_string(path).expect("failed to read JSON file");
    serde_json::from_str(&text).expect("failed to parse JSON file")
}

fn within_budget(started_at: Instant, max_seconds: f64) -> bool {
    max_seconds <= 0.0 || started_at.elapsed().as_secs_f64() < max_seconds
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CampaignConfig {
    n: u64,
    k_start: u64,
    k_batch_size: usize,
    sieve_limit: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    post_sieve_limit: Option<usize>,
    backend: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    export_dir: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CampaignCheckpoint {
    config: CampaignConfig,
    next_batch: u64,
    batches_completed: u64,
    total_base_survivors: u64,
    total_survivors: u64,
    total_hits: u64,
}

#[derive(Clone, Debug, Serialize)]
struct CampaignSummary {
    backend: String,
    n: u64,
    k_start: String,
    k_batch_size: usize,
    sieve_limit: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    post_sieve_limit: Option<usize>,
    elapsed_secs: f64,
    start_batch: u64,
    next_batch: u64,
    batches_completed: u64,
    total_base_survivors: u64,
    total_survivors: u64,
    total_hits: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    event_log: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    hit_log: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    checkpoint_path: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
struct BatchEvent {
    batch_index: u64,
    batch_k0: String,
    base_survivors: usize,
    survivors: usize,
    backend_result: BackendBatchResult,
}

fn usage() {
    eprintln!(
        "usage: fixed_n_campaign --n <n> [--backend <local_gmp|export_files|compact_export>] [--k-start <u64>] [--k-batch-size <usize>] [--sieve-limit <usize>] [--post-sieve-limit <usize>] [--max-seconds <f64>] [--max-batches <u64>] [--start-batch <u64>] [--stop-after-hits <u64>] [--status-every <u64>] [--event-log <file>] [--hit-log <file>] [--checkpoint-out <file>] [--resume-checkpoint <file>] [--json-out <file>] [--export-dir <dir>]"
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
    let backend_name = arg_value(&args, "--backend").unwrap_or_else(|| "local_gmp".to_string());
    let k_start = arg_value(&args, "--k-start")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(3);
    let k_batch_size = arg_value(&args, "--k-batch-size")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(65_536);
    let sieve_limit = arg_value(&args, "--sieve-limit")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(50_000);
    let post_sieve_limit = arg_value(&args, "--post-sieve-limit")
        .and_then(|s| s.parse::<usize>().ok());
    let max_seconds = arg_value(&args, "--max-seconds")
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(60.0);
    let max_batches = arg_value(&args, "--max-batches").and_then(|s| s.parse::<u64>().ok());
    let start_batch = arg_value(&args, "--start-batch")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);
    let stop_after_hits = arg_value(&args, "--stop-after-hits").and_then(|s| s.parse::<u64>().ok());
    let status_every = arg_value(&args, "--status-every")
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(5);
    let event_log = arg_value(&args, "--event-log").map(PathBuf::from);
    let hit_log = arg_value(&args, "--hit-log").map(PathBuf::from);
    let checkpoint_out = arg_value(&args, "--checkpoint-out").map(PathBuf::from);
    let resume_checkpoint = arg_value(&args, "--resume-checkpoint").map(PathBuf::from);
    let json_out = arg_value(&args, "--json-out").map(PathBuf::from);
    let export_dir = arg_value(&args, "--export-dir").map(PathBuf::from);

    let (
        config,
        mut next_batch,
        mut batches_completed,
        mut total_base_survivors,
        mut total_survivors,
        mut total_hits,
        checkpoint_path,
    ) =
        if let Some(path) = resume_checkpoint.as_deref() {
            let checkpoint: CampaignCheckpoint = load_json(path);
            (
                checkpoint.config,
                checkpoint.next_batch,
                checkpoint.batches_completed,
                checkpoint.total_base_survivors,
                checkpoint.total_survivors,
                checkpoint.total_hits,
                Some(path),
            )
        } else {
            (
                CampaignConfig {
                    n,
                    k_start: align_k_start(k_start),
                    k_batch_size,
                    sieve_limit,
                    post_sieve_limit,
                    backend: backend_name,
                    export_dir: export_dir
                        .as_deref()
                        .map(|path| path.to_string_lossy().to_string()),
                },
                start_batch,
                0,
                0,
                0,
                0,
                checkpoint_out.as_deref(),
            )
        };

    let export_dir = export_dir.or_else(|| config.export_dir.as_ref().map(PathBuf::from));
    let backend_box: Box<dyn FixedNBackend> = match config.backend.as_str() {
        "local_gmp" => Box::new(LocalGmpBackend),
        "export_files" => {
            let Some(dir) = export_dir else {
                eprintln!("--export-dir is required when --backend export_files");
                return;
            };
            Box::new(ExportFilesBackend::new(dir))
        }
        "compact_export" => {
            let Some(dir) = export_dir else {
                eprintln!("--export-dir is required when --backend compact_export");
                return;
            };
            Box::new(CompactExportBackend::new(dir))
        }
        other => {
            eprintln!("unknown backend: {}", other);
            return;
        }
    };

    println!("Fixed-n campaign runner");
    println!("  backend: {}", backend_box.name());
    println!("  n: {}", config.n);
    println!("  k_start: {}", config.k_start);
    println!("  k_batch_size: {}", config.k_batch_size);
    println!("  sieve_limit: {}", config.sieve_limit);
    if let Some(limit) = config.post_sieve_limit {
        println!("  post_sieve_limit: {}", limit);
    }
    println!("  start batch: {}", next_batch);

    let t0 = Instant::now();
    let start_batch_value = next_batch;
    let plan = build_k_sieve_plan(config.n, config.sieve_limit);
    let post_plan = config
        .post_sieve_limit
        .filter(|&limit| limit > config.sieve_limit)
        .map(|limit| build_k_sieve_plan_range(config.n, config.sieve_limit, limit));
    while within_budget(t0, max_seconds) {
        if let Some(limit) = max_batches {
            if batches_completed >= limit {
                break;
            }
        }
        if let Some(limit) = stop_after_hits {
            if total_hits >= limit {
                break;
            }
        }

        let Some(mut batch) = generate_survivor_batch_with_plan(
            config.n,
            config.k_start,
            config.k_batch_size,
            config.sieve_limit,
            next_batch,
            &plan,
        ) else {
            break;
        };
        let base_survivors = batch.survivors.len();
        total_base_survivors += base_survivors as u64;
        if let Some(ref extra_plan) = post_plan {
            batch.survivors = filter_survivors_with_plan(config.n, &batch.survivors, extra_plan);
            batch.sieve_limit = config.post_sieve_limit.unwrap_or(config.sieve_limit);
            batch.survival_rate *= survival_rate(extra_plan);
        }
        let result = backend_box.process_batch(
            &batch,
            t0,
            if max_seconds > 0.0 { Some(max_seconds) } else { None },
        );

        total_survivors += batch.survivors.len() as u64;
        total_hits += result.hit_count as u64;

        if let Some(path) = event_log.as_deref() {
            let event = BatchEvent {
                batch_index: batch.batch_index,
                batch_k0: batch.batch_k0.to_string(),
                base_survivors,
                survivors: batch.survivors.len(),
                backend_result: result.clone(),
            };
            append_json_line(path, &event);
        }
        if let Some(path) = hit_log.as_deref() {
            for hit in &result.hits {
                append_json_line(path, hit);
            }
        }

        next_batch += 1;
        batches_completed += 1;

        if let Some(path) = checkpoint_path {
            let checkpoint = CampaignCheckpoint {
                config: config.clone(),
                next_batch,
                batches_completed,
                total_base_survivors,
                total_survivors,
                total_hits,
            };
            save_json(path, &checkpoint);
        }

        if status_every > 0 && batches_completed % status_every == 0 {
            println!(
                "    [{:.1}s] next_batch={} batches={} survivors={} hits={}",
                t0.elapsed().as_secs_f64(),
                next_batch,
                batches_completed,
                total_survivors,
                total_hits
            );
        }

        if result.partial {
            break;
        }
    }

    let summary = CampaignSummary {
        backend: config.backend.clone(),
        n: config.n,
        k_start: config.k_start.to_string(),
        k_batch_size: config.k_batch_size,
        sieve_limit: config.sieve_limit,
        post_sieve_limit: config.post_sieve_limit,
        elapsed_secs: t0.elapsed().as_secs_f64(),
        start_batch: start_batch_value,
        next_batch,
        batches_completed,
        total_base_survivors,
        total_survivors,
        total_hits,
        event_log: event_log.as_deref().map(path_label),
        hit_log: hit_log.as_deref().map(path_label),
        checkpoint_path: checkpoint_path.map(path_label),
    };

    let json = serde_json::to_string_pretty(&summary).expect("failed to serialize summary");
    println!("{}", json);
    if let Some(path) = json_out {
        fs::write(path, json).expect("failed to write JSON");
    }
}
