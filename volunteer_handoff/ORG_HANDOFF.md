# Organization Handoff Note

This bundle is intended for an organization running volunteer or distributed prime-search workloads.

## What is included

- a deterministic assigned-batch worker
- a resident Rust worker that reuses sieve tables across multiple jobs
- a deterministic JSON work-unit format
- a deterministic JSON result format
- a validation and operations note

## What is not included

- a BOINC server
- account systems
- host reputation logic
- packaging/signing infrastructure
- full checkpoint/restart support

## Why this is the right handoff level

A public volunteer-computing project should not distribute the whole research repo.
It should distribute only the minimum executable workload model:

- fixed range
- fixed arithmetic pipeline
- fixed result format
- duplicate validation

That is what this bundle provides.

## Immediate evaluation path

1. Run `sample_work_unit.json`
2. Confirm deterministic result agreement across multiple machines
3. Prefer the resident Rust worker for repeated jobs
4. Generate larger work units with `make_work_unit.py`
5. Wrap the worker in the organization's scheduler or BOINC app wrapper
6. Add duplicate validation and retry handling

## External integration target

The current bundle now includes a resident Rust worker binary target.
For production deployment, the next step would be to harden it further so it:

- supports checkpoint/restart
- exposes a frozen command-line contract
- is packaged and signed per platform
