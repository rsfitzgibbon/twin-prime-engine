# OpenPFGW Backend Bridge

## Goal

Use the exact fixed-`n` front-end in this repo to generate survivor batches, then
hand those batches to an OpenPFGW-capable worker for special-form PRP testing.

The search side stays in Rust.
The large-PRP side moves to a specialized backend.

## Pieces

### 1. Exact front-end campaign

Use the campaign runner with:

- `--sieve-limit` for the base `k`-space filter
- `--post-sieve-limit` for the exact survivor-only post-sieve
- `--backend compact_export` for compact batch export

Example:

```powershell
fixed_n_campaign.exe --n 1288907 --backend compact_export --export-dir campaign_388k_exact_post200k --k-start 3 --k-batch-size 4096 --sieve-limit 50000 --post-sieve-limit 200000 --max-batches 4 --json-out campaign_388k_exact_post200k_summary.json
```

### 2. OpenPFGW bridge

[openpfgw_batch_runner.py](C:/Users/rsfit/twin-prime-engine/openpfgw_batch_runner.py)
consumes the compact export batches and reconstructs:

- `plus.txt`
- `minus.txt`
- `plus.abcd`
- `minus.abcd`
- `manifest.json`

These files are created per batch in a clean worker directory.

## Modes

### Prepare-only mode

Use this when the machine does not have `pfgw.exe` installed yet.

```powershell
python openpfgw_batch_runner.py --compact-dir campaign_388k_exact_post200k --output-dir openpfgw_jobs_388k --prepare-only --json-out openpfgw_jobs_388k_summary.json
```

### Execute mode

If a worker has OpenPFGW installed, the same script can invoke it.

```powershell
python openpfgw_batch_runner.py --compact-dir campaign_388k_exact_post200k --output-dir openpfgw_jobs_388k --pfgw-exe C:\path\to\pfgw.exe --json-out openpfgw_jobs_388k_run_summary.json
```

Default command template:

```text
"{exe}" "{input}"
```

If the local PFGW build needs different flags, override with:

- `--plus-template`
- `--minus-template`

Placeholders:

- `{exe}`
- `{input}`
- `{output}`

### Continuous worker mode

For a long-running worker on one machine:

```powershell
python openpfgw_worker.py --compact-dir campaign_388k_exact_post200k --output-dir openpfgw_jobs_388k --pfgw-exe C:\path\to\pfgw64.exe --state-file worker_state.json --result-log worker_results.jsonl --poll-seconds 30
```

This worker:

- scans for new compact batches
- processes each batch once
- persists a checkpoint in `worker_state.json`
- appends each completed batch result to `worker_results.jsonl`

## Why this is the right architecture

This keeps the repo split cleanly:

- the engine does exact residue filtering and exact post-sieving
- the worker side does specialized PRP work

That matters because the front-end is already fast enough at `388k+` digits.
The true bottleneck is the large-number backend.

## Current status

This repo now supports:

1. exact post-sieve survivor reduction
2. compact export batches
3. OpenPFGW-oriented batch reconstruction
4. optional execution/parsing if `pfgw.exe` is available

So the missing piece is no longer the handoff design.
The missing piece is the actual OpenPFGW binary on a worker machine.
