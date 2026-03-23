use serde_json::json;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::backend::{BackendArtifact, BackendBatchResult, FixedNBackend};
use crate::fixed_n::{formula_minus, formula_plus, SurvivorBatch};

fn path_label(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string_lossy().to_string())
}

fn write_text(path: &Path, text: String) {
    fs::write(path, text).expect("failed to write export file");
}

pub struct ExportFilesBackend {
    output_dir: PathBuf,
}

impl ExportFilesBackend {
    pub fn new(output_dir: PathBuf) -> Self {
        fs::create_dir_all(&output_dir).expect("failed to create export directory");
        Self { output_dir }
    }
}

impl FixedNBackend for ExportFilesBackend {
    fn name(&self) -> &'static str {
        "export_files"
    }

    fn process_batch(
        &self,
        batch: &SurvivorBatch,
        _started_at: Instant,
        _max_seconds: Option<f64>,
    ) -> BackendBatchResult {
        let stem = format!("n{}_batch{:012}", batch.n, batch.batch_index);
        let prefix = self.output_dir.join(stem);
        let k_file = prefix.with_extension("k.txt");
        let minus_file = prefix.with_extension("minus.txt");
        let plus_file = prefix.with_extension("plus.txt");
        let pairs_file = prefix.with_extension("pairs.tsv");

        let mut k_lines = String::new();
        let mut minus_lines = String::new();
        let mut plus_lines = String::new();
        let mut pairs_lines = String::from("k\tminus\tplus\n");

        for &k in &batch.survivors {
            let minus_expr = formula_minus(k, batch.n);
            let plus_expr = formula_plus(k, batch.n);
            k_lines.push_str(&format!("{}\n", k));
            minus_lines.push_str(&minus_expr);
            minus_lines.push('\n');
            plus_lines.push_str(&plus_expr);
            plus_lines.push('\n');
            pairs_lines.push_str(&format!("{}\t{}\t{}\n", k, minus_expr, plus_expr));
        }

        write_text(&k_file, k_lines);
        write_text(&minus_file, minus_lines);
        write_text(&plus_file, plus_lines);
        write_text(&pairs_file, pairs_lines);

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
                    kind: "k_list".to_string(),
                    label: path_label(&k_file),
                },
                BackendArtifact {
                    kind: "minus_formulas".to_string(),
                    label: path_label(&minus_file),
                },
                BackendArtifact {
                    kind: "plus_formulas".to_string(),
                    label: path_label(&plus_file),
                },
                BackendArtifact {
                    kind: "pairs_tsv".to_string(),
                    label: path_label(&pairs_file),
                },
            ],
            notes: vec!["Batch exported for external backend".to_string()],
            partial: false,
            metadata: Some(json!({
                "output_dir": path_label(&self.output_dir),
                "estimated_digits": batch.estimated_digits,
                "survival_rate": batch.survival_rate,
            })),
        }
    }
}
