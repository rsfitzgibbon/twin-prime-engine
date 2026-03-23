# 388k-Digit Fixed-n Backend Findings

## Target

We tested the fixed-`n` twin-prime pipeline past `388,000` decimal digits using

- `n = 1,288,907`
- family `k*2^n +/- 1`
- `k mod 6 = 3`

This exponent gives an estimated size of `388,001` digits for the first candidate batch.

## Main Result

The sieve still works cleanly at this size.

The local backend does not.

More precisely:

- `fixed_n_k_search --sieve-only` at `n = 1,288,907` reduced a `64`-candidate batch to `2` survivors (`k = 99`, `255`) immediately.
- `fixed_n_campaign --backend local_gmp --n 1288907 --k-start 99 --k-batch-size 1 --max-batches 1 --max-seconds 30`
  did **not finish even one survivor** before the outer process hit a `120s` wall-clock timeout.
- Therefore the current local GMP/BPSW backend is not a viable path for actually finding twin primes beyond `388k` digits on this machine.

## What Still Scales

The fixed-`n` `k`-sieve and export path still scale well.

Baseline 388k export test:

- backend: `export_files`
- `4` batches
- `4096` raw `k` values per batch
- `377` total survivors
- elapsed: `0.0284759s`

This established that the front-end is not the bottleneck.

## Speed Work Added

### 1. Campaign-level sieve-plan cache

The campaign runner was rebuilding the same small-prime `k`-sieve plan every batch.
That is now cached once per campaign and reused across all batches.

On the same 388k export workload:

- before cache reuse: `0.0284759s`
- after cache reuse: `0.0176196s`

That is a `1.62x` improvement for the repeated-batch front-end path.

### 2. Compact export backend

A new backend, `compact_export`, avoids writing duplicated giant formula text files.
Instead it writes:

- one metadata JSON file per batch
- one compact binary payload of delta-varint encoded `k` offsets

On the same 388k workload:

- `export_files`: `0.0176196s`, `31,414` bytes across `16` files
- `compact_export`: `0.0125743s`, `2,540` bytes across `8` files

Effects:

- `1.40x` faster than text export
- `12.37x` smaller on disk

This does not speed up primality itself, but it materially improves the handoff path to a specialized external backend.

## Conclusion

At `388k+` digits, the project is now in a split state:

- **Search front-end**: viable
- **Local probable-prime backend**: not viable

So the next serious work is not more sieve tuning. It is backend replacement.

## Defensible Next Moves

Ordered by value:

1. Integrate a specialized external PRP backend such as OpenPFGW/GWNUM-class tooling.
2. Keep survivor batches compact and binary; do not export repeated formula strings.
3. Move to a GPU-aware backend only if it keeps operands resident and overlaps transfers with computation.
4. Treat CPU threads as a pipeline problem, not just "use more threads":
   sieve/export, candidate staging, and PRP work should be decoupled.
5. Do **not** pursue fuzzy-logic prediction from the last 17 twin primes; that has no credible path to accelerating primality testing.
