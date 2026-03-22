"""
Reusable twin-prime finder for the current corridor + arithmetic framework.

Modes:
- fast: hard mask at y = 211
- high_precision: hard mask at y = 503
- exact_range: hard mask at the largest prime <= sqrt(N)

Optional reranking:
- soft_to = R > y ranks survivors by the first additional pair-blocker in (y, R]
- candidates with no blocker up to R rank first

Optional sieve scoring (from PSM step 7 insight):
- --score enables Selberg-style divisor-sum weights for continuous ranking
- --score-z sets the upper bound for the sieve weight (default: auto)
- Sieve weight uses inclusion-exclusion over ALL prime factors, not just first blocker
- Higher weight = fewer/smaller prime factors = higher confidence

This is useful as a finite-range predictor/ranker.
It is not a proof mechanism for infinitely many twin primes.
"""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

from spectral_twin_prime_v2 import SEED_A, SEED_B, build_bridge, generate_gasket, sieve_bool


MODE_CUTOFFS = {
    "fast": 211,
    "high_precision": 503,
}

# When --score is enabled without --score-z, auto-pick a z above the hard cutoff.
SCORE_Z_AUTO = {
    211: 503,
    503: 997,
}

RUST_FINDER_EXE = (
    Path(__file__).resolve().parent
    / "rust_twin_prime_finder"
    / "target"
    / "release"
    / ("rust_twin_prime_finder.exe" if os.name == "nt" else "rust_twin_prime_finder")
)


def prime_list_from_sieve(is_prime, limit_n):
    return [p for p in range(2, limit_n + 1) if is_prime[p]]


def largest_prime_leq(is_prime, limit_n):
    for p in range(limit_n, 1, -1):
        if is_prime[p]:
            return p
    return None


def resolve_cutoff(limit_n, mode, custom_cutoff=None):
    if custom_cutoff is not None:
        return int(custom_cutoff), "custom"
    if mode in MODE_CUTOFFS:
        return MODE_CUTOFFS[mode], mode
    if mode == "exact_range":
        root = int(math.isqrt(limit_n + 2))
        is_prime = sieve_bool(root)
        return largest_prime_leq(is_prime, root), mode
    raise ValueError(f"Unknown mode: {mode}")


def corridor_candidates(limit_n):
    return range(5, limit_n - 1, 6)


def pct(num, den):
    return 100.0 * num / den if den else 0.0


def first_pair_blocker(n, primes, low_exclusive, high_inclusive):
    q = n + 2
    for p in primes:
        if p <= low_exclusive:
            continue
        if p > high_inclusive:
            return None
        if n % p == 0 and n != p:
            return p
        if q % p == 0 and q != p:
            return p
    return None


def hard_mask_survivor(n, primes, cutoff):
    return first_pair_blocker(n, primes, 0, cutoff) is None


def soft_rank_key(n, primes, hard_cutoff, soft_to):
    if soft_to is None or soft_to <= hard_cutoff:
        return None, None
    blocker = first_pair_blocker(n, primes, hard_cutoff, soft_to)
    # No blocker is best; encode that with a key larger than any actual blocker.
    rank_value = soft_to + 1 if blocker is None else blocker
    return blocker, rank_value


def sieve_score(n, primes, score_z):
    """
    Selberg-style divisor-sum weight for pair (n, n+2).

    Adapted from psm_step7_pairwise_phi_experiment.py.
    Uses inclusion-exclusion over all prime factors of both endpoints up to score_z.
    Respects the self-prime exception: if an endpoint equals p, that is not a blocker.
    Higher weight = cleaner (fewer/smaller prime factors) = more likely twin prime.

    Returns log1p(w^2) where w is the truncated Möbius sum.
    """
    x, q = n, n + 2
    log_R = 2.0 * math.log(score_z)

    factors = []
    for p in primes:
        if p > score_z:
            break
        x_blocked = (x % p == 0 and x != p)
        q_blocked = (q % p == 0 and q != p)
        if x_blocked or q_blocked:
            factors.append(math.log(p))

    total = 1.0

    def rec(idx, logd, sign):
        nonlocal total
        for j in range(idx, len(factors)):
            nlogd = logd + factors[j]
            if nlogd > log_R + 1e-15:
                continue
            total += sign * (-1.0) * (1.0 - nlogd / log_R)
            rec(j + 1, nlogd, sign * (-1.0))

    rec(0, 0.0, 1.0)
    return math.log1p(total * total)


def resolve_score_z(hard_cutoff, score_z_arg, limit_n):
    """Pick the sieve-score upper bound z."""
    if score_z_arg is not None:
        return int(score_z_arg)
    if hard_cutoff in SCORE_Z_AUTO:
        return SCORE_Z_AUTO[hard_cutoff]
    # For exact_range or custom cutoffs, go to next round number above
    root = int(math.isqrt(limit_n + 2))
    candidate_z = hard_cutoff * 2
    return min(candidate_z, root) if candidate_z > hard_cutoff else hard_cutoff + 100


def geometry_role_builder(limit_n):
    A = generate_gasket(SEED_A, limit_n)
    B = generate_gasket(SEED_B, limit_n)
    bridge_ab, _ = build_bridge(A, B, limit_n)
    bridge_ba, _ = build_bridge(B, A, limit_n)

    def role_pair(n):
        if n in A:
            left = "A"
        elif n in bridge_ab:
            left = "AB"
        elif n in bridge_ba:
            left = "BA"
        else:
            left = "other"

        q = n + 2
        if q in B:
            right = "B"
        elif q in bridge_ab:
            right = "AB"
        elif q in bridge_ba:
            right = "BA"
        else:
            right = "other"
        return f"{left}+{right}"

    return role_pair


def ranking_label(result):
    if result.get("score_z"):
        return f"sieve_score(score_z={result['score_z']})"
    if result["soft_to"] is not None and result["soft_to"] > result["hard_cutoff"]:
        return f"soft_rerank(soft_to={result['soft_to']})"
    return "natural_order"


def rust_finder_available():
    return RUST_FINDER_EXE.exists()


def find_twin_candidates_rust(limit_n, mode="high_precision", soft_to=None, top_k=None,
                              custom_cutoff=None, score=False, score_z=None):
    if not rust_finder_available():
        raise FileNotFoundError(f"Rust finder binary not found: {RUST_FINDER_EXE}")

    with tempfile.NamedTemporaryFile("w+", suffix=".json", delete=False) as tmp:
        json_path = tmp.name

    try:
        cmd = [
            str(RUST_FINDER_EXE),
            "--limit", str(limit_n),
            "--mode", mode,
            "--preview", "0",
            "--json-out", json_path,
        ]
        if soft_to is not None:
            cmd += ["--soft-to", str(soft_to)]
        if top_k is not None:
            cmd += ["--top-k", str(top_k)]
        if custom_cutoff is not None:
            cmd += ["--cutoff", str(custom_cutoff)]
        if score:
            cmd += ["--score"]
        if score_z is not None:
            cmd += ["--score-z", str(score_z)]

        subprocess.run(cmd, check=True, capture_output=True, text=True)

        with open(json_path, "r", encoding="utf-8") as f:
            result = json.load(f)
        return result
    finally:
        try:
            os.unlink(json_path)
        except OSError:
            pass


def find_twin_candidates(limit_n, mode="high_precision", soft_to=None, top_k=None,
                         include_geometry=False, custom_cutoff=None,
                         score=False, score_z=None, engine="python"):
    if engine == "rust":
        if include_geometry:
            raise ValueError("Rust engine currently supports arithmetic candidate generation only (no geometry annotations).")
        return find_twin_candidates_rust(
            limit_n=limit_n,
            mode=mode,
            soft_to=soft_to,
            top_k=top_k,
            custom_cutoff=custom_cutoff,
            score=score,
            score_z=score_z,
        )
    if engine != "python":
        raise ValueError(f"Unknown engine: {engine}")

    hard_cutoff, resolved_mode = resolve_cutoff(limit_n, mode, custom_cutoff)

    # Determine the sieve-score upper bound if scoring is enabled.
    actual_score_z = None
    if score:
        actual_score_z = resolve_score_z(hard_cutoff, score_z, limit_n)

    max_prime_needed = max(hard_cutoff, soft_to or hard_cutoff, actual_score_z or hard_cutoff)
    is_prime = sieve_bool(max_prime_needed)
    primes = prime_list_from_sieve(is_prime, max_prime_needed)

    role_fn = geometry_role_builder(limit_n) if include_geometry else None

    candidates = []
    for n in corridor_candidates(limit_n):
        if not hard_mask_survivor(n, primes, hard_cutoff):
            continue
        extra_blocker, rank_value = soft_rank_key(n, primes, hard_cutoff, soft_to)
        record = {
            "n": n,
            "pair": [n, n + 2],
            "mode": resolved_mode,
            "hard_cutoff": hard_cutoff,
            "soft_to": soft_to,
            "extra_blocker": extra_blocker,
            "rank_value": rank_value,
        }
        if actual_score_z is not None:
            record["sieve_score"] = sieve_score(n, primes, actual_score_z)
            record["score_z"] = actual_score_z
        if role_fn is not None:
            record["geometry_role"] = role_fn(n)
        candidates.append(record)

    # Sort: sieve score (best first) > soft rerank > natural order.
    if actual_score_z is not None:
        candidates.sort(key=lambda x: x["sieve_score"], reverse=True)
    elif soft_to is not None and soft_to > hard_cutoff:
        candidates.sort(key=lambda x: (x["rank_value"], -x["n"]), reverse=True)
    else:
        candidates.sort(key=lambda x: x["n"])

    if top_k is not None:
        candidates = candidates[: int(top_k)]

    for idx, item in enumerate(candidates, start=1):
        item["rank"] = idx

    return {
        "limit_n": limit_n,
        "mode": resolved_mode,
        "hard_cutoff": hard_cutoff,
        "soft_to": soft_to,
        "score_z": actual_score_z,
        "candidate_count": len(candidates),
        "engine": "python",
        "ranking": ranking_label(
            {
                "hard_cutoff": hard_cutoff,
                "soft_to": soft_to,
                "score_z": actual_score_z,
            }
        ),
        "candidates": candidates,
    }


def write_csv(result, csv_path):
    candidates = result["candidates"]
    if not candidates:
        print("No candidates to export.")
        return
    fieldnames = ["rank", "n", "left", "right", "mode", "hard_cutoff", "ranking"]
    if result.get("score_z"):
        fieldnames += ["sieve_score", "score_z"]
    if result["soft_to"] is not None and result["soft_to"] > result["hard_cutoff"]:
        fieldnames += ["extra_blocker", "rank_value"]
    if "geometry_role" in candidates[0]:
        fieldnames.append("geometry_role")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for item in candidates:
            row = {
                "rank": item["rank"],
                "n": item["n"],
                "left": item["pair"][0],
                "right": item["pair"][1],
                "mode": item["mode"],
                "hard_cutoff": item["hard_cutoff"],
                "ranking": result["ranking"],
            }
            if "sieve_score" in item:
                row["sieve_score"] = f"{item['sieve_score']:.6f}"
                row["score_z"] = item["score_z"]
            if result["soft_to"] is not None and result["soft_to"] > result["hard_cutoff"]:
                row["extra_blocker"] = item["extra_blocker"]
                row["rank_value"] = item["rank_value"]
            if "geometry_role" in item:
                row["geometry_role"] = item["geometry_role"]
            writer.writerow(row)
    print(f"Saved CSV to {csv_path} ({len(candidates):,} rows)")


def print_summary(result, preview_count):
    print(f"Twin Prime Finder  N = {result['limit_n']:,}")
    print(f"Mode: {result['mode']}  hard_cutoff={result['hard_cutoff']}")
    print(f"Engine: {result.get('engine', 'python')}")
    print(f"Ranking: {result['ranking']}")
    if result.get("score_z"):
        print(f"Sieve scoring: yes  score_z={result['score_z']}")
    if result["soft_to"] is not None and result["soft_to"] > result["hard_cutoff"]:
        print(f"Soft rerank: yes  soft_to={result['soft_to']}")
    print(f"Predicted candidate pairs: {result['candidate_count']:,}")

    preview = result["candidates"][:preview_count]
    if not preview:
        print("No candidates.")
        return

    print(f"\nTop {len(preview)} candidates:")
    for item in preview:
        bits = [f"pair=({item['pair'][0]}, {item['pair'][1]})"]
        if "sieve_score" in item:
            bits.append(f"score={item['sieve_score']:.4f}")
        if item.get("soft_to") and item["soft_to"] > item["hard_cutoff"]:
            bits.append(f"extra_blocker={item['extra_blocker']}")
        if "geometry_role" in item:
            bits.append(f"role={item['geometry_role']}")
        print("  " + ", ".join(bits))


def evaluate_result(result):
    limit_n = result["limit_n"]
    is_prime = sieve_bool(limit_n + 2)
    candidates = result["candidates"]
    twin_count = 0
    for item in candidates:
        left, right = item["pair"]
        if is_prime[left] and is_prime[right]:
            twin_count += 1

    total_twins = sum(
        1 for n in corridor_candidates(limit_n)
        if is_prime[n] and is_prime[n + 2]
    )
    candidate_count = len(candidates)
    precision = twin_count / candidate_count if candidate_count else 0.0
    recall = twin_count / total_twins if total_twins else 0.0
    f1 = 0.0 if precision + recall == 0 else 2.0 * precision * recall / (precision + recall)
    total_corridor = max(0, (limit_n - 5) // 6 + 1)
    return {
        "limit_n": limit_n,
        "mode": result["mode"],
        "hard_cutoff": result["hard_cutoff"],
        "soft_to": result["soft_to"],
        "score_z": result["score_z"],
        "ranking": result["ranking"],
        "candidate_count": candidate_count,
        "true_twins": twin_count,
        "total_twins": total_twins,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "elimination_rate": 1.0 - candidate_count / total_corridor if total_corridor else 0.0,
    }


def benchmark_finder(limits, modes, soft_to=None, top_k=None,
                     include_geometry=False, custom_cutoff=None,
                     score=False, score_z=None, engine="python"):
    rows = []
    for limit_n in limits:
        for mode in modes:
            result = find_twin_candidates(
                limit_n=limit_n,
                mode=mode,
                soft_to=soft_to,
                top_k=top_k,
                include_geometry=include_geometry,
                custom_cutoff=custom_cutoff,
                score=score,
                score_z=score_z,
                engine=engine,
            )
            rows.append(evaluate_result(result))
    return rows


def write_benchmark_csv(rows, csv_path):
    fieldnames = [
        "limit_n",
        "mode",
        "hard_cutoff",
        "soft_to",
        "score_z",
        "ranking",
        "candidate_count",
        "true_twins",
        "total_twins",
        "precision",
        "recall",
        "f1",
        "elimination_rate",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["precision"] = f"{row['precision']:.6f}"
            out["recall"] = f"{row['recall']:.6f}"
            out["f1"] = f"{row['f1']:.6f}"
            out["elimination_rate"] = f"{row['elimination_rate']:.6f}"
            writer.writerow(out)
    print(f"Saved benchmark CSV to {csv_path} ({len(rows):,} rows)")


def print_benchmark(rows):
    print("Twin Prime Finder Benchmark")
    print(
        f"\n  {'N':>10} | {'mode':>12} | {'hard':>6} | {'candidates':>11} | "
        f"{'twins':>8} | {'precision':>10} | {'recall':>8} | {'F1':>8}"
    )
    print(
        f"  {'-'*10}-+-{'-'*12}-+-{'-'*6}-+-{'-'*11}-+-{'-'*8}-+-{'-'*10}-+-{'-'*8}-+-{'-'*8}"
    )
    for row in rows:
        print(
            f"  {row['limit_n']:>10,} | {row['mode']:>12} | {row['hard_cutoff']:>6} | "
            f"{row['candidate_count']:>11,} | {row['true_twins']:>8,} | "
            f"{pct(row['precision'], 1):>9.4f}% | {pct(row['recall'], 1):>7.2f}% | "
            f"{pct(row['f1'], 1):>7.2f}%"
        )


def main():
    parser = argparse.ArgumentParser(description="Twin-prime finder for the corridor + arithmetic framework.")
    parser.add_argument("--limit", type=int, default=1_000_000, help="Upper bound N.")
    parser.add_argument("--mode", choices=["fast", "high_precision", "exact_range"], default="high_precision")
    parser.add_argument("--cutoff", type=int, default=None, help="Custom hard-mask cutoff y.")
    parser.add_argument("--soft-to", type=int, default=None, help="Optional rerank cutoff R > y.")
    parser.add_argument("--top-k", type=int, default=None, help="Return only the top k candidates.")
    parser.add_argument("--score", action="store_true", help="Enable Selberg-style sieve scoring for ranking.")
    parser.add_argument("--score-z", type=int, default=None, help="Upper bound z for sieve scoring (default: auto).")
    parser.add_argument("--include-geometry", action="store_true", help="Annotate candidates with geometry roles.")
    parser.add_argument("--engine", choices=["python", "rust"], default="python",
                        help="Execution engine for unscored candidate generation.")
    parser.add_argument("--preview", type=int, default=20, help="How many candidates to print.")
    parser.add_argument("--json-out", type=str, default=None, help="Optional output JSON path.")
    parser.add_argument("--csv-out", type=str, default=None, help="Optional output CSV path.")
    parser.add_argument("--benchmark", action="store_true", help="Run accuracy benchmarks instead of listing candidates.")
    parser.add_argument("--benchmark-limits", type=int, nargs="+", default=[100_000, 500_000, 1_000_000],
                        help="Benchmark limits to evaluate.")
    parser.add_argument("--benchmark-modes", nargs="+",
                        choices=["fast", "high_precision", "exact_range"],
                        default=["fast", "high_precision", "exact_range"],
                        help="Finder modes to benchmark.")
    parser.add_argument("--benchmark-json-out", type=str, default=None, help="Optional benchmark JSON path.")
    parser.add_argument("--benchmark-csv-out", type=str, default=None, help="Optional benchmark CSV path.")
    args = parser.parse_args()

    if args.benchmark:
        rows = benchmark_finder(
            limits=args.benchmark_limits,
            modes=args.benchmark_modes,
            soft_to=args.soft_to,
            top_k=args.top_k,
            include_geometry=args.include_geometry,
            custom_cutoff=args.cutoff,
            score=args.score,
            score_z=args.score_z,
            engine=args.engine,
        )
        print_benchmark(rows)

        if args.benchmark_csv_out:
            write_benchmark_csv(rows, args.benchmark_csv_out)

        if args.benchmark_json_out:
            with open(args.benchmark_json_out, "w", encoding="utf-8") as f:
                json.dump(rows, f, indent=2)
            print(f"Saved benchmark JSON to {args.benchmark_json_out}")
        return

    result = find_twin_candidates(
        limit_n=args.limit,
        mode=args.mode,
        soft_to=args.soft_to,
        top_k=args.top_k,
        include_geometry=args.include_geometry,
        custom_cutoff=args.cutoff,
        score=args.score,
        score_z=args.score_z,
        engine=args.engine,
    )
    print_summary(result, args.preview)

    if args.csv_out:
        write_csv(result, args.csv_out)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print(f"Saved JSON to {args.json_out}")


if __name__ == "__main__":
    main()
