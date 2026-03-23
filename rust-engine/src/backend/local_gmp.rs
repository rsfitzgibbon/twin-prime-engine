use rayon::prelude::*;
use rug::integer::IsPrime;
use rug::Assign;
use rug::Integer;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::Instant;

use crate::backend::{BackendBatchResult, FixedNBackend, ProbablePrimeHit};
use crate::fixed_n::SurvivorBatch;

const PROTH_BASES: [u32; 8] = [2, 3, 5, 7, 11, 13, 17, 19];

fn within_budget(started_at: Instant, max_seconds: Option<f64>) -> bool {
    match max_seconds {
        Some(limit) if limit > 0.0 => started_at.elapsed().as_secs_f64() < limit,
        _ => true,
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

    fn test_candidate(&mut self, k: u64, n: u64) -> (bool, Option<ProbablePrimeHit>) {
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

        let p_str = self.minus.to_string();
        let q_str = self.plus.to_string();
        let p_len = p_str.len();
        let q_len = q_str.len();
        (
            true,
            Some(ProbablePrimeHit {
                k: k.to_string(),
                digits: p_len.max(q_len),
                plus_mode: plus_mode.unwrap_or_else(|| "bpsw_fallback".to_string()),
                p_head: p_str[..40.min(p_len)].to_string(),
                p_tail: p_str[p_len.saturating_sub(40)..].to_string(),
                q_head: q_str[..40.min(q_len)].to_string(),
                q_tail: q_str[q_len.saturating_sub(40)..].to_string(),
            }),
        )
    }
}

pub struct LocalGmpBackend;

impl FixedNBackend for LocalGmpBackend {
    fn name(&self) -> &'static str {
        "local_gmp"
    }

    fn process_batch(
        &self,
        batch: &SurvivorBatch,
        started_at: Instant,
        max_seconds: Option<f64>,
    ) -> BackendBatchResult {
        let expired = AtomicBool::new(false);
        let plus_counter = AtomicU64::new(0);
        let minus_counter = AtomicU64::new(0);

        let hits: Vec<ProbablePrimeHit> = batch
            .survivors
            .par_iter()
            .map_init(TestCtx::new, |ctx, &k| {
                if expired.load(Ordering::Relaxed) {
                    return None;
                }
                if !within_budget(started_at, max_seconds) {
                    expired.store(true, Ordering::Relaxed);
                    return None;
                }

                plus_counter.fetch_add(1, Ordering::Relaxed);
                let (minus_tested, hit) = ctx.test_candidate(k, batch.n);
                if minus_tested {
                    minus_counter.fetch_add(1, Ordering::Relaxed);
                }
                hit
            })
            .filter_map(|item| item)
            .collect();

        BackendBatchResult {
            backend: self.name().to_string(),
            batch_index: batch.batch_index,
            batch_k0: batch.batch_k0.to_string(),
            survivors_in: batch.survivors.len(),
            plus_tested: plus_counter.load(Ordering::Relaxed),
            minus_tested: minus_counter.load(Ordering::Relaxed),
            hit_count: hits.len(),
            hits,
            artifacts: Vec::new(),
            notes: vec!["Local GMP/BPSW backend".to_string()],
            partial: expired.load(Ordering::Relaxed),
            metadata: None,
        }
    }
}
