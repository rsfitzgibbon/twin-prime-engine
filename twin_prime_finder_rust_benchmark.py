"""
Benchmark Python vs Rust execution engines for twin_prime_finder.py.

This measures the real top-level API path:
    find_twin_candidates(..., engine="python" | "rust")

So the Rust timings include subprocess invocation and JSON parsing overhead,
not just the raw standalone executable runtime.
"""

import json
import time

import twin_prime_finder as finder


WORKLOADS = [
    {"limit_n": 1_000_000, "mode": "fast"},
    {"limit_n": 1_000_000, "mode": "high_precision"},
    {"limit_n": 1_000_000, "mode": "exact_range"},
]

REPEATS = 7


def time_engine(limit_n, mode, engine):
    best = None
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        result = finder.find_twin_candidates(limit_n=limit_n, mode=mode, engine=engine)
        elapsed = time.perf_counter() - t0
        best = elapsed if best is None else min(best, elapsed)
    return best, result["candidate_count"]


def main():
    print("Twin Prime Finder Rust Benchmark")
    print(f"Repeats per engine: {REPEATS}")
    print()

    rows = []
    for workload in WORKLOADS:
        limit_n = workload["limit_n"]
        mode = workload["mode"]

        py_time, py_count = time_engine(limit_n, mode, "python")
        rust_time, rust_count = time_engine(limit_n, mode, "rust")
        if py_count != rust_count:
            raise RuntimeError(
                f"Candidate count mismatch for {mode}: python={py_count}, rust={rust_count}"
            )

        speedup_factor = py_time / rust_time if rust_time else float("inf")
        speedup_percent = 100.0 * (speedup_factor - 1.0)

        row = {
            "limit_n": limit_n,
            "mode": mode,
            "python_best_seconds": py_time,
            "rust_best_seconds": rust_time,
            "candidate_count": py_count,
            "speedup_factor": speedup_factor,
            "speedup_percent": speedup_percent,
            "clears_12_percent_goal": speedup_percent >= 12.0,
        }
        rows.append(row)

        print(
            f"{mode:>14} @ N={limit_n:,}: "
            f"python={py_time*1000:.2f} ms, rust={rust_time*1000:.2f} ms, "
            f"speedup={speedup_factor:.2f}x ({speedup_percent:.2f}%)"
        )

    out = {
        "repeats": REPEATS,
        "workloads": rows,
        "notes": [
            "Rust timings include subprocess invocation and JSON parsing through the Python API",
            "The goal threshold is a 12% speedup over the Python engine",
        ],
    }
    out_path = "twin_prime_finder_rust_benchmark_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print()
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
