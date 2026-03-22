# Rust Finder Acceleration

Date: 2026-03-22

## Scope

The Rust port targets the actual hot path in the practical finder:

- corridor candidate generation
- hard small-prime masking
- optional soft rerank by next blocker

It does **not** currently replace:

- geometry-role annotation
- sieve-score ranking

Those still run through the Python engine.

## New files

- `rust_twin_prime_finder/Cargo.toml`
- `rust_twin_prime_finder/src/main.rs`
- `twin_prime_finder_rust_benchmark.py`
- `twin_prime_finder_rust_benchmark_results.json`

The existing Python CLI in `twin_prime_finder.py` now supports:

- `--engine python`
- `--engine rust`

for the unscored / non-geometry path.

## Why this target

Profiling showed the dominant cost in the Python finder was:

- `first_pair_blocker`

That function dominated `find_twin_candidates()` on the `high_precision`
workload, so the Rust port replaced the blocker scan and survivor generation
rather than touching secondary features first.

## Measured result

Benchmarked through the **same Python API path**:

- `find_twin_candidates(..., engine="python")`
- `find_twin_candidates(..., engine="rust")`

This means the Rust timings include:

- subprocess startup
- JSON serialization / parsing

So the comparison is stricter than a raw binary-only microbenchmark.

At `N = 1,000,000`:

- `fast`
  - Python: `224.35 ms`
  - Rust: `77.68 ms`
  - Speedup: `2.89x` (`188.83%`)
- `high_precision`
  - Python: `281.77 ms`
  - Rust: `59.12 ms`
  - Speedup: `4.77x` (`376.64%`)
- `exact_range`
  - Python: `464.45 ms`
  - Rust: `54.57 ms`
  - Speedup: `8.51x` (`751.17%`)

So the Rust engine clears the stated `12%` usefulness threshold by a wide margin on the practical hard-mask path.

## Honest limitation

This is **not** a benchmark against PrimeGrid.

There is no PrimeGrid executable, workload harness, or agreed comparison
protocol in the workspace, so no honest claim should be made that this
outspeeds PrimeGrid itself.

What is proven here is narrower and real:

- the Rust engine substantially outperforms the current Python engine for the
  practical twin-prime finder path used in this repo

## Practical commands

```powershell
python Newtheory\twin_prime_finder.py --limit 1000000 --mode high_precision --engine rust --preview 20
python Newtheory\twin_prime_finder.py --benchmark --benchmark-limits 100000 500000 1000000 --engine rust
python Newtheory\twin_prime_finder_rust_benchmark.py
```

## Next engineering step

If continuing this branch, the next worthwhile Rust target is:

- soft-score / ranking support in Rust, so the `high_precision --score --top-k`
  path can also move off Python

That is the remaining practical ranking path with measurable value.
