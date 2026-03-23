# Exact Replacement for the Gasket Ranking Branch

## Problem

The gasket/bridge branch is useful for coverage and structural interpretation,
but it does not improve large twin-prime search enough to justify sitting in the
hot path.

In particular:

- residue coverage from the gapless gasket is effectively universal modulo `2520`
- pure fuzzy or geometry-first ranking is weaker than arithmetic scoring
- the real bottleneck for large searches is the expensive PRP backend

So the correct replacement is not another soft geometry heuristic.
It is an **exact post-sieve** that removes more composites before PRP work.

## Replacement

The Rust fixed-`n` campaign runner now supports:

- base sieve limit: `--sieve-limit`
- exact post-sieve limit: `--post-sieve-limit`

The post-sieve uses the same exact `k*2^n +/- 1` residue arithmetic as the base
filter, but only on the already-small survivor list. That makes it a cheap,
correct second-stage composite remover.

Implemented in:

- [fixed_n.rs](C:/Users/rsfit/twin-prime-engine/rust-engine/src/fixed_n.rs)
- [fixed_n_campaign.rs](C:/Users/rsfit/twin-prime-engine/rust-engine/src/bin/fixed_n_campaign.rs)

## Benchmarks

### 388k-digit target

Test:

- `n = 1,288,907`
- base limit `50,000`
- post-sieve limit `200,000`
- `4` batches of `4096` raw `k` values

Results:

- base only: `377` survivors
- with exact post-sieve: `296` survivors

That is a reduction of `81` survivors, or about `21.5%`.

This matters because those survivors are exactly what an external large-PRP
backend would have to test.

### Live local backend check

Test:

- `n = 1000`
- `5` batches of `10,000` raw `k` values
- local GMP backend

Results:

- base only:
  - `1032` survivors
  - `1032` plus-side PRP tests
  - `30` minus-side PRP tests
  - elapsed `0.2748s`

- with exact post-sieve to `200,000`:
  - `805` survivors
  - `805` plus-side PRP tests
  - `27` minus-side PRP tests
  - elapsed `0.2997s`

Interpretation:

- PRP workload dropped by about `22%`
- wall-clock time stayed roughly flat at this small size because the extra sieve
  work is still visible
- at large sizes, where PRP dominates and one test is expensive, this trade is
  favorable

## Conclusion

This exact post-sieve is the right replacement for the gasket ranking branch in
the large-search engine:

- exact
- cheap on the survivor list
- directly reduces the expensive backend workload
- useful for the actual goal of finding the next largest twin prime

Geometry remains useful as a coverage/explanation layer, but not as the primary
search accelerator.
