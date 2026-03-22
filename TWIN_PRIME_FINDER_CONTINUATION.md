# Twin Prime Finder Continuation

Date: 2026-03-22

## What was added

The practical finder branch was continued in `twin_prime_finder.py`.

New capabilities:

- `--benchmark` mode for direct accuracy sweeps across multiple limits and modes
- benchmark CSV / JSON export
- richer candidate CSV export with:
  - `rank`
  - `ranking`
  - `rank_value` for soft rerank mode
  - existing sieve-score / geometry columns preserved

This continues the handoff direction instead of reopening already-tested dead ends.

## Main benchmark artifact

Saved:

- `twin_prime_finder_benchmark_results.csv`
- `twin_prime_finder_benchmark_results.json`

Results:

- `N = 100,000`
  - `fast`: `92.9331%` precision, `100%` recall
  - `high_precision`: `100%` precision, `100%` recall
  - `exact_range`: `100%` precision, `100%` recall
- `N = 500,000`
  - `fast`: `68.0077%` precision, `100%` recall
  - `high_precision`: `95.5012%` precision, `100%` recall
  - `exact_range`: `100%` precision, `100%` recall
- `N = 1,000,000`
  - `fast`: `59.1413%` precision, `100%` recall
  - `high_precision`: `85.5288%` precision, `100%` recall
  - `exact_range`: `100%` precision, `100%` recall

These match the earlier experimental picture and are now exposed through the reusable CLI.

## New ranking insight

The nontrivial continuation was to test whether sieve scoring helps once we stop looking only at the very first exact-by-cutoff segment.

Saved comparison artifacts:

- `twin_prime_finder_unscored_top4000_benchmark.csv`
- `twin_prime_finder_unscored_top4000_benchmark.json`
- `twin_prime_finder_scored_top4000_benchmark.csv`
- `twin_prime_finder_scored_top4000_benchmark.json`

At `N = 1,000,000`, `mode = high_precision`, `top_k = 4,000`:

- natural order:
  - `3,895 / 4,000` true twins
  - precision `97.3750%`
  - recall `47.69%`
- sieve-scored order:
  - `4,000 / 4,000` true twins
  - precision `100.0000%`
  - recall `48.97%`

Interpretation:

- scoring does **not** create new survivors
- but it does improve the **ordering** of the existing high-precision pool
- the gain only becomes visible once we move beyond the trivial low-`n` region where natural order is already perfect

## Current practical state

Best operational picture:

- use `fast` for a cheap wide filter
- use `high_precision` for a strong all-range heuristic
- use `high_precision --score --top-k ...` when you want the best tranche first
- use `exact_range` only when you want finite-range exactness

## Useful commands

```powershell
python Newtheory\twin_prime_finder.py --limit 1000000 --mode high_precision --score --preview 20
python Newtheory\twin_prime_finder.py --limit 1000000 --mode high_precision --score --top-k 4000 --csv-out Newtheory\top4000_scored.csv
python Newtheory\twin_prime_finder.py --benchmark --benchmark-limits 100000 500000 1000000
python Newtheory\twin_prime_finder.py --benchmark --benchmark-limits 1000000 --benchmark-modes high_precision --score --top-k 4000
```
