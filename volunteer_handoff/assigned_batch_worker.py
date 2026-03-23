import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import sys
import time
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import twin_prime_engine as eng


def load_work_unit(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        unit = json.load(f)

    required = [
        "version",
        "job_id",
        "target_digits",
        "m_start",
        "batch_size",
        "sieve_depth",
        "sieve_limit",
    ]
    missing = [key for key in required if key not in unit]
    if missing:
        raise ValueError(f"work unit missing fields: {missing}")
    if unit["sieve_depth"] not in {"base", "extended"}:
        raise ValueError("sieve_depth must be base or extended")
    return unit


def build_sieve_tables(sieve_depth: str, sieve_limit: int):
    all_primes = eng.prime_sieve(sieve_limit)
    primes_small = [int(p) for p in all_primes if 5 <= p <= 1_000_000]
    inv6_small = {p: pow(6, -1, p) for p in primes_small}

    if sieve_depth == "base":
        primes_ext = np.array([], dtype=np.int64)
        inv6_ext = np.array([], dtype=np.int64)
    else:
        primes_ext_list = [int(p) for p in all_primes if p > 1_000_000]
        primes_ext = np.array(primes_ext_list, dtype=np.int64)
        inv6_ext = np.array([pow(6, -1, p) for p in primes_ext_list], dtype=np.int64)

    return primes_small, inv6_small, primes_ext, inv6_ext


def run_assigned_batch(unit: dict) -> dict:
    m_start = int(unit["m_start"])
    batch_size = int(unit["batch_size"])
    sieve_depth = unit["sieve_depth"]
    sieve_limit = int(unit["sieve_limit"])

    t0 = time.perf_counter()
    primes_small, inv6_small, primes_ext, inv6_ext = build_sieve_tables(sieve_depth, sieve_limit)
    setup_seconds = time.perf_counter() - t0

    m_start_mpz = eng.gmpy2.mpz(m_start)

    alive = eng.base_sieve_fast(m_start_mpz, m_start, batch_size, primes_small, inv6_small)
    if sieve_depth == "extended":
        alive = eng.extended_sieve_fast(alive, m_start_mpz, batch_size, primes_ext, inv6_ext)

    offsets = np.where(alive)[0]
    m_values = [m_start + int(off) for off in offsets]

    worker_t0 = time.perf_counter()
    found_pairs = []
    if m_values:
        with ThreadPoolExecutor(max_workers=eng.NUM_WORKERS) as pool:
            for result in pool.map(eng.test_candidate, m_values, chunksize=1):
                if result is None:
                    continue
                m_val, p_val, q_val = result
                found_pairs.append(
                    {
                        "m": str(m_val),
                        "p": str(p_val),
                        "q": str(q_val),
                        "digits": len(str(p_val)),
                    }
                )

    found_pairs.sort(key=lambda row: row["m"])
    elapsed = time.perf_counter() - worker_t0

    return {
        "job_id": unit["job_id"],
        "status": "ok",
        "target_digits": int(unit["target_digits"]),
        "m_start": str(m_start),
        "batch_size": batch_size,
        "sieve_depth": sieve_depth,
        "setup_seconds": setup_seconds,
        "elapsed_seconds": elapsed,
        "total_raw": batch_size,
        "survivors": len(m_values),
        "tested": len(m_values),
        "found_count": len(found_pairs),
        "found_pairs": found_pairs,
    }


def main():
    parser = argparse.ArgumentParser(description="Run one deterministic assigned batch for twin-prime search.")
    parser.add_argument("--work-unit", type=Path, required=True, help="Input work unit JSON.")
    parser.add_argument("--result-out", type=Path, required=True, help="Output result JSON.")
    args = parser.parse_args()

    unit = load_work_unit(args.work_unit)
    result = run_assigned_batch(unit)

    args.result_out.parent.mkdir(parents=True, exist_ok=True)
    with args.result_out.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"Saved result to {args.result_out}")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
