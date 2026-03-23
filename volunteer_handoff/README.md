# Volunteer Handoff Bundle

This folder packages the twin-prime search engine into a deterministic work-unit model suitable for a BOINC or PrimeGrid-style organization.

The goal is not to send the whole repo to volunteers.
The goal is to send:

- a small worker binary or script
- a deterministic work unit
- a deterministic result format
- a validation policy

## Included

- `../rust-engine/src/bin/volunteer_worker.rs`
  - preferred resident Rust worker that reuses sieve tables across multiple jobs
- `assigned_batch_worker.py`
  - reference Python worker for one assigned batch
- `make_work_unit.py`
  - creates a JSON work unit with fixed parameters
- `WORK_UNIT_SPEC.md`
  - work-unit and result schema
- `VALIDATION_AND_OPERATIONS.md`
  - quorum, duplication, and deployment guidance
- `sample_work_unit.json`
  - small example unit

## Why this matters

PrimeGrid- or BOINC-style deployment needs work that is:

- deterministic
- resumable or at least chunked
- independently verifiable
- small enough to assign, retry, and duplicate

The existing engine already has the right arithmetic core.
What external organizations need is the assignment protocol.

## Reference usage

Create a work unit:

```powershell
python make_work_unit.py --digits 1000 --m-start 166666666666666666 --batch-size 200000 --sieve-depth extended --out sample_work_unit.json
```

Run the assigned batch:

```powershell
python assigned_batch_worker.py --work-unit sample_work_unit.json --result-out sample_result.json
```

Build and run the resident Rust worker:

```powershell
cd ..\rust-engine
cargo build --release --bin volunteer_worker
.\target\x86_64-pc-windows-gnu\release\volunteer_worker.exe --work-unit ..\volunteer_handoff\sample_work_unit.json --result-dir ..\volunteer_handoff\rust_results
```

## Scope

This is a reference handoff package, not a full BOINC server implementation.

An organization integrating it would still need:

- scheduler / feeder
- result validation
- host reputation handling
- packaging and signing
- platform-specific binaries

The preferred external deliverable is now the resident Rust worker, not the Python reference worker.
