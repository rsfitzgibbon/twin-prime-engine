# Specialized Backend Plan

## Objective

Build a backend specialized for fixed-`n` twin-prime searches of the form:

- `k * 2^n - 1`
- `k * 2^n + 1`

with:

- `n` fixed for a campaign
- `k` varying over the `k == 3 (mod 6)` progression
- continuous operation
- resumable execution
- exportable and verifiable results

The immediate target is not "solve the Twin Prime Conjecture". The target is a backend that can continuously search larger and larger fixed-`n` twin-prime candidates, including million-digit scale, without the current GMP/BPSW path becoming the terminal bottleneck.

## Current State

The repo already has the right front-end pieces:

- `twin_prime_k2n_hunter.py`
  - Python proof-of-concept for sieving on `k`
- `rust-engine/src/bin/fixed_n_k_search.rs`
  - exact validation mode
  - direct `k`-space sieve
  - continuous checkpointed search
  - continuous survivor export
- `volunteer_handoff/`
  - deterministic work-unit framing for distributed use

The current limit is clear:

- the `k`-sieve works at million-digit scale
- the local GMP probable-prime backend does not

So the specialized backend must replace or sit behind the final primality stage, not the sieve stage.

## Design Principle

Do not build a full large-prime FFT engine from scratch as step one.

That is too expensive, too risky, and too slow to validate.

Instead:

1. keep the Rust `k`-sieve front-end
2. define a backend interface for large-number testing
3. first integrate an external specialized PRP engine
4. only then decide whether a native custom backend is justified

This creates a usable system quickly and keeps the long-term option open for a fully native backend.

## Backend Requirements

The specialized backend must support:

- `k * 2^n +/- 1` inputs directly
- fixed-`n` campaign mode
- batch processing of many `k` values
- checkpoint/resume
- deterministic output records
- probable-prime status plus verification metadata
- optional proof generation when supported
- continuous worker mode
- separation between "candidate generation" and "final certification"

Preferred properties:

- FFT-based arithmetic
- disk-backed checkpoints
- low per-candidate startup overhead
- clean worker CLI or library interface
- easy integration with distributed scheduling

## Proposed Architecture

### 1. Sieve Front-End

Owner:

- `rust-engine/src/bin/fixed_n_k_search.rs`

Responsibilities:

- generate `k` candidates
- apply fixed-`n` residue blocking
- checkpoint `next_batch`
- emit survivor batches

This layer is already mostly done.

### 2. Backend Adapter Layer

New module:

- `rust-engine/src/backend/`

Responsibilities:

- translate survivor batches into backend-specific jobs
- run backend workers
- parse backend outputs
- normalize results into a repo-owned JSON schema

This should be implemented behind a small trait-like interface:

```rust
trait TwinPrpBackend {
    fn name(&self) -> &'static str;
    fn submit_batch(&self, batch: SurvivorBatch) -> BackendJob;
    fn poll(&self, job: &BackendJob) -> BackendJobState;
    fn collect(&self, job: BackendJob) -> BackendBatchResult;
}
```

The repo should not depend on a single backend forever. The adapter boundary is what keeps that flexible.

### 3. Campaign Coordinator

New binary:

- `rust-engine/src/bin/fixed_n_campaign.rs`

Responsibilities:

- drive continuous search
- choose between local PRP or external backend
- own campaign metadata:
  - `n`
  - current batch
  - sieve limit
  - backend name
  - output paths
- keep append-only logs
- handle retry and resume

This becomes the "always-on" process.

### 4. Result Store

Store format:

- JSONL for append-only event logs
- periodic checkpoint JSON
- optional SQLite for indexed querying

Recommended event types:

- `batch_started`
- `batch_sieved`
- `backend_submitted`
- `backend_completed`
- `probable_prime_found`
- `proof_generated`
- `batch_failed`

### 5. Certification Layer

Separate probable-prime hits from final proof/certification.

Stages:

1. `PRP hit`
2. `retest on independent worker/backend`
3. `proof/certificate if available`
4. `publishable record`

This is important because large-scale search systems fail if they treat one raw backend result as final truth.

## Recommended Implementation Sequence

### Phase 0: Lock the interfaces

Deliverables:

- `SurvivorBatch` schema
- `BackendBatchResult` schema
- campaign checkpoint schema
- append-only result log schema

Success gate:

- one schema used by both local and external backends

### Phase 1: Add backend abstraction

Deliverables:

- `backend/mod.rs`
- `backend/local_gmp.rs`
- conversion of current `fixed_n_k_search.rs` test path to the backend interface

Success gate:

- current local GMP/BPSW path still works through the new abstraction

### Phase 2: Add first external backend adapter

Target:

- pick one real external large-PRP engine and support it first

Pragmatic choice:

- file-based adapter if the backend is CLI-oriented
- process-based adapter if it supports long-running jobs

Responsibilities:

- write backend input files
- invoke backend
- parse status and result files
- normalize output

Success gate:

- one survivor batch round-trips end-to-end through the external backend

### Phase 3: Continuous campaign runner

Deliverables:

- `fixed_n_campaign.rs`
- config file for campaign settings
- checkpointed loop
- append-only logs

Success gate:

- stop and resume without losing batch position
- continuous operation for hours without manual intervention

### Phase 4: Distributed execution

Deliverables:

- split campaign into work units
- backend-ready assignment protocol
- duplicate verification
- worker lease/retry model

This can reuse ideas already present in `volunteer_handoff/`.

Success gate:

- two workers process disjoint batches and produce consistent normalized outputs

### Phase 5: Native specialized backend

Only do this if Phases 1-4 prove the search architecture is sound.

Possible native scope:

- specialized modpow / transform path for `k * 2^n +/- 1`
- persistent FFT plans for fixed `n`
- batch-wise transform reuse
- custom checkpoint format

This is the expensive part and should come last.

## What the Native Backend Must Exploit

If a custom backend is built, it should exploit structure that generic big-int testing does not:

- fixed exponent `n`
- shared transform length across a campaign
- repeated `k * 2^n +/- 1` form
- paired testing of `-1` and `+1`
- batch locality

That means the backend should try to reuse:

- transform plans
- scratch buffers
- twiddle tables
- checkpoint format
- candidate staging memory

If it cannot reuse those, it is probably not better than an external specialized engine.

## Short-Term Deliverables

These are the next concrete steps for this repo:

1. Add a backend interface module under `rust-engine/src/backend/`.
2. Refactor the current local PRP path in `fixed_n_k_search.rs` to use that interface.
3. Add a new campaign runner binary for continuous backend-driven execution.
4. Define repo-owned JSON schemas for survivor batches and backend results.
5. Pick one external backend target and write the first adapter.

## Success Metrics

The backend project is successful if it reaches these milestones:

- can run continuously for 24h without manual intervention
- can resume from checkpoint without repeating completed batches
- can export and re-import backend jobs deterministically
- can verify hits on an independent second pass
- can process million-digit candidate batches in a routine way
- can scale from local run to distributed run without changing the core schemas

## Non-Goals

This backend project is not trying to:

- prove infinitude of twin primes
- replace all prime-search software immediately
- support every prime form at once
- build a generic arbitrary-precision algebra system

It is specifically about a reliable, continuous backend for fixed-`n` twin-prime campaigns.

## Recommendation

The best next move is:

- do not build a custom FFT backend first
- build the backend interface first
- integrate one external PRP engine second
- make continuous campaign mode stable third
- only then decide whether a native backend is worth the cost

That path gives the repo a usable large-scale engine soonest, while keeping the door open for a truly specialized backend later.
