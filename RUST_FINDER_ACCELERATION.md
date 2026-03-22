# Rust Finder Acceleration

Date: 2026-03-22

## Scope

The Rust port targets the actual hot path in the practical finder:

- corridor candidate generation
- hard small-prime masking
- optional soft rerank by next blocker
- optional sieve-score ranking

It does **not** currently replace:

- geometry-role annotation

Geometry annotations still run through the Python engine.

## New files

- `rust_twin_prime_finder/Cargo.toml`
- `rust_twin_prime_finder/src/main.rs`
- `twin_prime_finder_rust_benchmark.py`
- `twin_prime_finder_rust_benchmark_results.json`

The existing Python CLI in `twin_prime_finder.py` now supports:

- `--engine python`
- `--engine rust`

for the arithmetic finder path, including `--score`.

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
  - Python: `110.60 ms`
  - Rust: `44.28 ms`
  - Speedup: `2.50x` (`149.79%`)
- `high_precision`
  - Python: `137.54 ms`
  - Rust: `35.99 ms`
  - Speedup: `3.82x` (`282.14%`)
- `exact_range`
  - Python: `175.25 ms`
  - Rust: `32.90 ms`
  - Speedup: `5.33x` (`432.73%`)
- `high_precision --score`
  - Python: `277.48 ms`
  - Rust: `46.87 ms`
  - Speedup: `5.92x` (`492.06%`)
- `high_precision --score --top-k 4000`
  - Python: `272.19 ms`
  - Rust: `32.53 ms`
  - Speedup: `8.37x` (`736.86%`)

So the Rust engine clears the stated `12%` usefulness threshold by a wide margin on both the hard-mask and scored-ranking paths.

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
python twin_prime_finder.py --limit 1000000 --mode high_precision --engine rust --preview 20
python twin_prime_finder.py --limit 1000000 --mode high_precision --engine rust --score --top-k 4000 --preview 20
python twin_prime_finder.py --benchmark --benchmark-limits 100000 500000 1000000 --engine rust
python twin_prime_finder_rust_benchmark.py
```

## Next engineering step

If continuing this branch, the next worthwhile Rust target is:

- geometry-role annotation in the Rust result path
- reducing subprocess / JSON overhead by exposing the Rust finder as an importable extension or library call

The scored path is already off Python; the remaining gap is feature parity and interface overhead.
