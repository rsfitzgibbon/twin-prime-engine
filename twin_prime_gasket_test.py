"""
Head-to-head: Baseline vs Selberg-scored vs Gasket-aligned twin prime search.

Strategy A (Baseline): Random batch, sieve to 5×10^7, SPRP in arrival order.
Strategy B (Selberg):  Same sieve, then quick trial divide survivors against
                       200 extra primes. Candidates passing extra check tested
                       first (2-tier priority).
Strategy C (Gasket):   Use Apollonian gasket curvature residues mod 2520 as a
                       secondary admissibility filter during sieve. Only test
                       m-values whose 6m-1 residue mod 2520 matches a gasket-
                       or bridge-covered class.

Metric: SPRP tests to first twin (fewer = more efficient targeting).
"""

import math
import random
import sys
import time
import json
from multiprocessing.pool import ThreadPool

import gmpy2
import numpy as np

sys.set_int_max_str_digits(1_000_000)

C2 = 0.6601618158
NUM_WORKERS = 12

# ---- Apollonian gasket generation (from spectral_twin_prime_v2) ----

SEED_A = (-1, 2, 2, 3)
SEED_B = (-2, 3, 6, 7)

def generate_gasket(seed, limit_n, max_depth=50):
    curvatures = set()
    for s in seed:
        if 0 < s <= limit_n:
            curvatures.add(s)
    queue = [tuple(seed)]
    seen = {tuple(sorted(seed))}
    depth = 0
    while queue and depth < max_depth:
        next_queue = []
        for quad in queue:
            vals = list(quad)
            quad_sum = sum(vals)
            for i in range(4):
                new_val = 2 * (quad_sum - vals[i]) - vals[i]
                if 0 < new_val <= limit_n:
                    new_quad = vals.copy()
                    new_quad[i] = new_val
                    key = tuple(sorted(new_quad))
                    if key not in seen:
                        seen.add(key)
                        next_queue.append(tuple(new_quad))
                        curvatures.add(new_val)
        queue = next_queue
        depth += 1
    return curvatures

def build_bridge(A, B, limit_n):
    import bisect
    a_sorted = sorted(A)
    bridge = set()
    for h in range(1, limit_n + 1):
        if h in B:
            continue
        idx = bisect.bisect_left(a_sorted, h)
        left = [a_sorted[j] for j in range(max(0, idx - 3), idx) if a_sorted[j] < h]
        if not left:
            continue
        a = max(left)
        c = 2 * h - a
        if 1 <= c <= limit_n and c not in (A | B) and c not in bridge:
            bridge.add(c)
    return bridge


# ---- Build gasket residue tables ----

def build_gasket_residue_set(modulus=2520):
    """Generate gasket + bridge curvatures, extract residue classes mod modulus.
    Returns set of residues r such that r mod modulus is 'gasket-covered'."""
    limit = 100_000  # enough curvatures to saturate residue classes
    A = generate_gasket(SEED_A, limit)
    B = generate_gasket(SEED_B, limit)
    bridge_ab = build_bridge(A, B, limit)
    bridge_ba = build_bridge(B, A, limit)
    covered = A | B | bridge_ab | bridge_ba

    # Which residue classes mod modulus appear in the covered set?
    residues = set()
    for c in covered:
        residues.add(c % modulus)

    return residues, len(A), len(B), len(bridge_ab), len(bridge_ba), len(covered)


# ---- Sieve ----

def prime_sieve(limit):
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
    return np.where(sieve)[0]


def base_sieve(m_start, batch_size, primes_small, inv6_small):
    alive = np.ones(batch_size, dtype=bool)
    for p in primes_small:
        inv6 = inv6_small[p]
        r1 = (inv6 - m_start % p) % p
        alive[r1::p] = False
        r2 = ((p - inv6) - m_start % p) % p
        alive[r2::p] = False
    return alive


def extended_sieve(alive, m_start, batch_size, primes_ext, inv6_ext):
    m_mod = np.array([m_start % int(p) for p in primes_ext], dtype=np.int64)
    r1_all = (inv6_ext - m_mod) % primes_ext
    r2_all = (primes_ext - inv6_ext - m_mod) % primes_ext
    K = max(1, batch_size // int(primes_ext.min()) + 1)
    K = min(K, 50)
    for k in range(K):
        pos1 = r1_all + k * primes_ext
        valid1 = pos1 < batch_size
        if np.any(valid1):
            alive[pos1[valid1]] = False
        pos2 = r2_all + k * primes_ext
        valid2 = pos2 < batch_size
        if np.any(valid2):
            alive[pos2[valid2]] = False
    return alive


# ---- Candidate tester ----

def test_candidate(m):
    p1 = gmpy2.mpz(6) * gmpy2.mpz(m) - 1
    if not gmpy2.is_strong_prp(p1, 2):
        return None
    p2 = p1 + 2
    if not gmpy2.is_strong_prp(p2, 2):
        return None
    if gmpy2.is_bpsw_prp(p1) and gmpy2.is_bpsw_prp(p2):
        return (int(m), int(p1), int(p2))
    return None


# ---- Strategy A: Baseline (random order) ----

def strategy_baseline(m_values, pool):
    chunk = max(1, len(m_values) // (NUM_WORKERS * 4))
    tested = 0
    for result in pool.imap_unordered(test_candidate, m_values, chunksize=chunk):
        tested += 1
        if result is not None:
            return result, tested
    return None, tested


# ---- Strategy B: Selberg-scored (extra trial division priority) ----

def selberg_score_batch(m_values, extra_primes):
    """Score each m by how many extra primes DON'T divide 6m-1 or 6m+1.
    Higher score = fewer small factors beyond sieve = more likely prime."""
    scores = []
    for m in m_values:
        p1 = 6 * m - 1
        p2 = 6 * m + 1
        clean = 0
        for ep in extra_primes:
            if p1 % ep != 0 and p2 % ep != 0:
                clean += 1
        scores.append(clean)
    return scores


def strategy_selberg(m_values, pool, extra_primes):
    # Score candidates
    scores = selberg_score_batch(m_values, extra_primes)
    # Sort by score descending (most "clean" first)
    indexed = sorted(zip(scores, m_values), reverse=True)
    sorted_m = [m for _, m in indexed]

    chunk = max(1, len(sorted_m) // (NUM_WORKERS * 4))
    tested = 0
    for result in pool.imap_unordered(test_candidate, sorted_m, chunksize=chunk):
        tested += 1
        if result is not None:
            return result, tested
    return None, tested


# ---- Strategy C: Gasket-aligned (residue filter) ----

def strategy_gasket(m_values, pool, gasket_residues, modulus):
    """Filter m_values to only those where (6m-1) mod modulus is gasket-covered.
    Then test filtered set. If no twin found, fall back to the rest."""
    aligned = []
    rest = []
    for m in m_values:
        r = (6 * m - 1) % modulus
        if r in gasket_residues:
            aligned.append(m)
        else:
            rest.append(m)

    # Test aligned first
    chunk = max(1, len(aligned) // (NUM_WORKERS * 4)) if aligned else 1
    tested = 0
    if aligned:
        for result in pool.imap_unordered(test_candidate, aligned, chunksize=chunk):
            tested += 1
            if result is not None:
                return result, tested, len(aligned), len(rest)

    # Fallback to rest
    chunk = max(1, len(rest) // (NUM_WORKERS * 4)) if rest else 1
    if rest:
        for result in pool.imap_unordered(test_candidate, rest, chunksize=chunk):
            tested += 1
            if result is not None:
                return result, tested, len(aligned), len(rest)

    return None, tested, len(aligned), len(rest)


# ---- Main benchmark ----

def run_benchmark(target_digits, max_seconds, primes_small, inv6_small,
                  primes_ext, inv6_ext, extra_primes, gasket_residues,
                  gasket_modulus, pool):

    m_low = 10 ** (target_digits - 1) // 6
    m_high = 10 ** target_digits // 6
    batch_size = min(5_000_000, (m_high - m_low) // 20)
    batch_size = max(batch_size, 500_000)

    lnN = target_digits * math.log(10)
    hl_raw = lnN * lnN / (2 * C2)
    print(f"  Target: ~{target_digits} digits, HL trials: {int(hl_raw):,}, batch: {batch_size:,}")

    # Run each strategy with equal time budget
    budget_each = max_seconds / 3
    results = {}

    for name, run_fn in [("baseline", "A"), ("selberg", "B"), ("gasket", "C")]:
        t0 = time.time()
        total_raw = 0
        total_surv = 0
        total_tested = 0
        found_result = None
        batches = 0
        extra_info = {}

        while time.time() - t0 < budget_each:
            m_start = random.randint(m_low, m_high - batch_size)
            alive = base_sieve(m_start, batch_size, primes_small, inv6_small)
            alive = extended_sieve(alive, m_start, batch_size, primes_ext, inv6_ext)
            ext_count = int(np.count_nonzero(alive))
            total_surv += ext_count
            total_raw += batch_size

            offsets = np.where(alive)[0]
            if len(offsets) == 0:
                batches += 1
                continue

            m_values = [m_start + int(off) for off in offsets]

            if run_fn == "A":
                result, tested = strategy_baseline(m_values, pool)
            elif run_fn == "B":
                result, tested = strategy_selberg(m_values, pool, extra_primes)
            elif run_fn == "C":
                out = strategy_gasket(m_values, pool, gasket_residues, gasket_modulus)
                result, tested = out[0], out[1]
                if len(out) > 2:
                    extra_info["last_aligned"] = out[2]
                    extra_info["last_rest"] = out[3]

            total_tested += tested
            batches += 1

            if result is not None:
                found_result = result
                break

        elapsed = time.time() - t0

        r = {
            "strategy": name,
            "elapsed": elapsed,
            "total_raw": total_raw,
            "total_surv": total_surv,
            "total_tested": total_tested,
            "batches": batches,
            "found": found_result is not None,
        }
        if found_result:
            m_val, p_val, q_val = found_result
            sp = str(p_val)
            r["digits"] = len(sp)
            r["p_head"] = sp[:40]
            r["p_tail"] = sp[-40:]
            r["m"] = m_val
        r.update(extra_info)

        results[name] = r

        status = f"FOUND in {elapsed:.2f}s, {total_tested:,} SPRP tests" if found_result else f"TIMEOUT {elapsed:.1f}s, {total_tested:,} tests"
        extra_str = ""
        if name == "gasket" and "last_aligned" in extra_info:
            extra_str = f" (aligned:{extra_info['last_aligned']}, rest:{extra_info['last_rest']})"
        print(f"    [{name:>8}] {status}{extra_str}")

    return results


def main():
    print("Twin Prime Gasket Benchmark")
    print("=" * 65)
    print("A=Baseline  B=Selberg-scored  C=Gasket-aligned")
    print()

    # Build sieve
    t0 = time.time()
    SIEVE_LIMIT = 50_000_000
    all_primes = prime_sieve(SIEVE_LIMIT)
    primes_small = [int(p) for p in all_primes if 5 <= p <= 1_000_000]
    inv6_small = {p: pow(6, -1, p) for p in primes_small}
    primes_ext_list = [int(p) for p in all_primes if p > 1_000_000]
    primes_ext = np.array(primes_ext_list, dtype=np.int64)
    inv6_ext = np.array([pow(6, -1, p) for p in primes_ext_list], dtype=np.int64)
    print(f"Sieve: {len(primes_small):,} + {len(primes_ext):,} primes ({time.time()-t0:.2f}s)")

    # Extra primes for Selberg scoring (200 primes just above sieve limit)
    extra_sieve = prime_sieve(SIEVE_LIMIT + 100_000)
    extra_primes = [int(p) for p in extra_sieve if p > SIEVE_LIMIT][:200]
    print(f"Selberg extra primes: {len(extra_primes)} ({extra_primes[0]:,} to {extra_primes[-1]:,})")

    # Gasket residue set
    t_g = time.time()
    GASKET_MOD = 2520  # LCM(24, 2520 = 2^3 * 3^2 * 5 * 7)
    gasket_residues, nA, nB, nAB, nBA, nCov = build_gasket_residue_set(GASKET_MOD)
    gasket_coverage = len(gasket_residues) / GASKET_MOD * 100
    print(f"Gasket: |A|={nA:,}, |B|={nB:,}, |bridge_AB|={nAB:,}, |bridge_BA|={nBA:,}, |covered|={nCov:,}")
    print(f"Gasket mod {GASKET_MOD}: {len(gasket_residues)}/{GASKET_MOD} residue classes covered ({gasket_coverage:.1f}%)")
    print(f"Gasket setup: {time.time()-t_g:.2f}s")
    print()

    pool = ThreadPool(NUM_WORKERS)

    targets = [
        (100, 90),     # 30s each strategy
        (500, 180),    # 60s each
        (1000, 360),   # 120s each
        (2000, 900),   # 300s each
    ]

    all_results = []

    for td, budget in targets:
        print(f"{'='*65}")
        r = run_benchmark(td, budget,
                          primes_small=primes_small, inv6_small=inv6_small,
                          primes_ext=primes_ext, inv6_ext=inv6_ext,
                          extra_primes=extra_primes,
                          gasket_residues=gasket_residues,
                          gasket_modulus=GASKET_MOD,
                          pool=pool)
        all_results.append({"target_digits": td, "strategies": r})
        print()

    pool.close()
    pool.join()

    # Summary table
    print("=" * 65)
    print("SUMMARY: SPRP tests to first twin")
    print("=" * 65)
    print(f"{'Digits':>8} | {'Baseline':>14} | {'Selberg':>14} | {'Gasket':>14}")
    print(f"{'-'*8}-+-{'-'*14}-+-{'-'*14}-+-{'-'*14}")
    for entry in all_results:
        td = entry["target_digits"]
        strats = entry["strategies"]
        cells = []
        for name in ["baseline", "selberg", "gasket"]:
            s = strats[name]
            if s["found"]:
                cells.append(f"{s['total_tested']:>8,} ({s['elapsed']:.1f}s)")
            else:
                cells.append(f"{'TIMEOUT':>8} ({s['elapsed']:.0f}s)")
        print(f"{td:>8} | {cells[0]:>14} | {cells[1]:>14} | {cells[2]:>14}")

    print()
    print("Efficiency = fewer SPRP tests per twin found = better targeting")

    # Winner analysis
    print()
    for entry in all_results:
        td = entry["target_digits"]
        strats = entry["strategies"]
        found = {n: s for n, s in strats.items() if s["found"]}
        if found:
            best = min(found.items(), key=lambda x: x[1]["total_tested"])
            print(f"  {td}-digit winner: {best[0]} ({best[1]['total_tested']:,} tests, {best[1]['elapsed']:.2f}s)")
        else:
            print(f"  {td}-digit: all timed out")

    with open("twin_prime_gasket_test_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print()
    print("Saved to twin_prime_gasket_test_results.json")


if __name__ == "__main__":
    main()
