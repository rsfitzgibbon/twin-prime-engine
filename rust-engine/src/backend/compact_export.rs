use serde_json::json;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::backend::{BackendArtifact, BackendBatchResult, FixedNBackend};
use crate::fixed_n::SurvivorBatch;

fn path_label(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string_lossy().to_string())
}

fn encode_varint(mut value: u64, out: &mut Vec<u8>) {
    while value >= 0x80 {
        out.push(((value as u8) & 0x7f) | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

fn write_compact_payload(path: &Path, batch: &SurvivorBatch) {
    let file = File::create(path).expect("failed to create compact export payload");
    let mut writer = BufWriter::new(file);
    let mut payload = Vec::with_capacity(batch.survivors.len() * 2 + 32);

    payload.extend_from_slice(&(batch.survivors.len() as u64).to_le_bytes());

    let mut prev_step = 0u64;
    for (i, &k) in batch.survivors.iter().enumerate() {
        let step = (k - batch.batch_k0) / 6;
        let delta = if i == 0 { step } else { step - prev_step };
        encode_varint(delta, &mut payload);
        prev_step = step;
    }

    writer
        .write_all(&payload)
        .expect("failed to write compact export payload");
    writer.flush().expect("failed to flush compact export payload");
}

pub struct CompactExportBackend {
    output_dir: PathBuf,
}

impl CompactExportBackend {
    pub fn new(output_dir: PathBuf) -> Self {
        fs::create_dir_all(&output_dir).expect("failed to create compact export directory");
        Self { output_dir }
    }
}

impl FixedNBackend for CompactExportBackend {
    fn name(&self) -> &'static str {
        "compact_export"
    }

    fn process_batch(
        &self,
        batch: &SurvivorBatch,
        _started_at: Instant,
        _max_seconds: Option<f64>,
    ) -> BackendBatchResult {
        let stem = format!("n{}_batch{:012}", batch.n, batch.batch_index);
        let prefix = self.output_dir.join(stem);
        let meta_file = prefix.with_extension("meta.json");
        let payload_file = prefix.with_extension("kdelta.bin");

        let meta = json!({
            "format": "fixed_n_compact_export_v1",
            "family": "k*2^n±1",
            "n": batch.n,
            "k_start": batch.k_start.to_string(),
            "k_batch_size": batch.k_batch_size,
            "sieve_limit": batch.sieve_limit,
            "batch_index": batch.batch_index,
            "batch_k0": batch.batch_k0.to_string(),
            "estimated_digits": batch.estimated_digits,
            "survival_rate": batch.survival_rate,
            "total_input_k": batch.total_input_k,
            "survivor_count": batch.survivors.len(),
            "encoding": {
                "kind": "delta_varint_steps",
                "step_unit": "6",
                "payload": path_label(&payload_file),
                "payload_layout": "u64_le survivor_count, then varint deltas of (k-batch_k0)/6",
            }
        });

        fs::write(
            &meta_file,
            serde_json::to_string_pretty(&meta).expect("failed to serialize compact export metadata"),
        )
        .expect("failed to write compact export metadata");
        write_compact_payload(&payload_file, batch);

        BackendBatchResult {
            backend: self.name().to_string(),
            batch_index: batch.batch_index,
            batch_k0: batch.batch_k0.to_string(),
            survivors_in: batch.survivors.len(),
            plus_tested: 0,
            minus_tested: 0,
            hit_count: 0,
            hits: Vec::new(),
            artifacts: vec![
                BackendArtifact {
                    kind: "metadata".to_string(),
                    label: path_label(&meta_file),
                },
                BackendArtifact {
                    kind: "kdelta_payload".to_string(),
                    label: path_label(&payload_file),
                },
            ],
            notes: vec![
                "Compact survivor export for external large-PRP backend".to_string(),
                "Stores k offsets only; formulas are reconstructed downstream".to_string(),
            ],
            partial: false,
            metadata: Some(json!({
                "output_dir": path_label(&self.output_dir),
                "estimated_digits": batch.estimated_digits,
                "survival_rate": batch.survival_rate,
            })),
        }
    }
}
