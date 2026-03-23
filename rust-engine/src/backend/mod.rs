pub mod compact_export;
pub mod export_files;
pub mod local_gmp;

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::time::Instant;

use crate::fixed_n::SurvivorBatch;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProbablePrimeHit {
    pub k: String,
    pub digits: usize,
    pub plus_mode: String,
    pub p_head: String,
    pub p_tail: String,
    pub q_head: String,
    pub q_tail: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BackendArtifact {
    pub kind: String,
    pub label: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BackendBatchResult {
    pub backend: String,
    pub batch_index: u64,
    pub batch_k0: String,
    pub survivors_in: usize,
    pub plus_tested: u64,
    pub minus_tested: u64,
    pub hit_count: usize,
    pub hits: Vec<ProbablePrimeHit>,
    pub artifacts: Vec<BackendArtifact>,
    pub notes: Vec<String>,
    pub partial: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub metadata: Option<Value>,
}

pub trait FixedNBackend {
    fn name(&self) -> &'static str;
    fn process_batch(
        &self,
        batch: &SurvivorBatch,
        started_at: Instant,
        max_seconds: Option<f64>,
    ) -> BackendBatchResult;
}
