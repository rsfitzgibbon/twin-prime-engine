"""
Twin prime search engine v2 -- improved from gasket benchmark findings.

Improvements over v1:
1. Sieve to 10^8 (from Selberg insight: extra trial division integrated into sieve)
2. gmpy2.mpz for big-int mod in sieve (2-4x faster m_start%p at large digits)
3. Two-tier sieve: base loop to 10^6, vectorized extension to 10^8
4. Adaptive batch sizing per digit target
5. ThreadPool x12 (gmpy2 releases GIL)
6. SPRP(2) on p1 first, p2 only for survivors, BPPSW confirmation

Benchmark showed: Selberg scoring uses 15x fewer SPRP tests but scoring overhead
kills wall-clock. Fix: bake the extra trial division INTO the sieve (zero per-candidate
overhead). Gasket filter covers 100% of mod-2520 classes (universal covering, no-op).
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


def prime_sieve(limit):
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p*p::p] = False
    return np.where(sieve)[0]


# ---- Fast sieve with gmpy2 acceleration ----

def base_sieve_fast(m_start_mpz, m_start_int, batch_size, primes_small, inv6_small):
    """Base sieve using gmpy2 for fast big-int mod."""
    alive = np.ones(batch_size, dtype=bool)
    for p in primes_small:
        inv6 = inv6_small[p]
        m_mod_p = int(gmpy2.f_mod(m_start_mpz, p))
        r1 = (inv6 - m_mod_p) % p
        alive[r1::p] = False
        r2 = ((p - inv6) - m_mod_p) % p
        alive[r2::p] = False
    return alive


def extended_sieve_fast(alive, m_start_mpz, batch_size, primes_ext, inv6_ext):
    """Vectorized sieve for large primes using gmpy2 for the mod computation."""
    # Compute m_start % p for all extended primes using gmpy2 (faster than Python int)
    m_mod = np.empty(len(primes_ext), dtype=np.int64)
    for i, p in enumerate(primes_ext):
        m_mod[i] = int(gmpy2.f_mod(m_start_mpz, int(p)))

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
    """SPRP(2) on p1 first, then p2. BPPSW confirmation for survivors."""
    p1 = gmpy2.mpz(6) * gmpy2.mpz(m) - 1
    if not gmpy2.is_strong_prp(p1, 2):
        return None
    p2 = p1 + 2
    if not gmpy2.is_strong_prp(p2, 2):
        return None
    if gmpy2.is_bpsw_prp(p1) and gmpy2.is_bpsw_prp(p2):
        return (int(m), int(p1), int(p2))
    return None


# ---- Main search ----

def search_twins(target_digits, max_seconds, primes_small, inv6_small,
                 primes_ext, inv6_ext, pool):
    m_low = 10 ** (target_digits - 1) // 6
    m_high = 10 ** target_digits // 6

    lnN = target_digits * math.log(10)
    hl_raw = lnN * lnN / (2 * C2)

    batch_size = min(5_000_000, (m_high - m_low) // 20)
    batch_size = max(batch_size, 500_000)

    print(f"  Target: ~{target_digits} digits, HL trials: {int(hl_raw):,}")
    print(f"  Batch: {batch_size:,}, Workers: {NUM_WORKERS}")

    t0 = time.time()
    total_raw = 0
    surv = 0
    tested = 0
    batches = 0

    while time.time() - t0 < max_seconds:
        m_start = random.randint(m_low, m_high - batch_size)
        m_start_mpz = gmpy2.mpz(m_start)

        t_sieve = time.time()
        alive = base_sieve_fast(m_start_mpz, m_start, batch_size, primes_small, inv6_small)
        alive = extended_sieve_fast(alive, m_start_mpz, batch_size, primes_ext, inv6_ext)
        sieve_time = time.time() - t_sieve

        ext_count = int(np.count_nonzero(alive))
        surv += ext_count
        total_raw += batch_size

        offsets = np.where(alive)[0]
        if len(offsets) == 0:
            batches += 1
            continue

        m_values = [m_start + int(off) for off in offsets]
        tested += len(m_values)

        chunk = max(1, len(m_values) // (NUM_WORKERS * 4))
        for result in pool.imap_unordered(test_candidate, m_values, chunksize=chunk):
            if result is not None:
                m_val, p_val, q_val = result
                elapsed = time.time() - t0
                sp = str(p_val)
                return {
                    "found": True, "m": m_val, "p": p_val, "q": q_val,
                    "digits": len(sp), "p_head": sp[:40], "p_tail": sp[-40:],
                    "elapsed": elapsed, "total_raw": total_raw,
                    "surv": surv, "tested": tested, "batches": batches + 1,
                    "sieve_time": sieve_time,
                }

        batches += 1
        if batches % 3 == 0 or time.time() - t0 > max_seconds * 0.9:
            el = time.time() - t0
            rate = tested / el if el > 0 else 0
            print(f"    [{el:.1f}s] {total_raw:,} raw -> {surv:,} sieved -> {tested:,} tested ({rate:.0f}/s) [sieve {sieve_time:.1f}s]")

    elapsed = time.time() - t0
    return {
        "found": False, "digits": target_digits, "elapsed": elapsed,
        "total_raw": total_raw, "surv": surv, "tested": tested,
        "batches": batches,
    }


def main():
    print("Twin Prime Engine v2")
    print(f"ThreadPool x{NUM_WORKERS} | SPRP(2)+BPPSW | Sieve 10^8 | gmpy2-accelerated")
    print()

    # Build sieve data -- extended to 10^8 (Selberg insight)
    t_pre = time.time()
    SIEVE_LIMIT = 100_000_000
    all_primes = prime_sieve(SIEVE_LIMIT)
    primes_small = [int(p) for p in all_primes if 5 <= p <= 1_000_000]
    inv6_small = {p: pow(6, -1, p) for p in primes_small}
    primes_ext_list = [int(p) for p in all_primes if p > 1_000_000]
    primes_ext = np.array(primes_ext_list, dtype=np.int64)
    inv6_ext = np.array([pow(6, -1, p) for p in primes_ext_list], dtype=np.int64)
    setup_time = time.time() - t_pre
    print(f"Sieve: {len(primes_small):,} primes to 10^6, {len(primes_ext):,} to 10^8 ({setup_time:.2f}s)")

    # Survival rate estimate
    surv = 1.0
    for p in primes_small:
        surv *= (p - 2) / p
    for p in primes_ext_list[:20000]:
        surv *= (p - 2) / p
    print(f"Estimated sieve survival: ~{surv*100:.4f}%")
    print()

    pool = ThreadPool(NUM_WORKERS)

    targets = [
        (100, 30),
        (500, 60),
        (1000, 300),
        (2000, 900),
        (5000, 3600),
    ]

    results = []
    largest = None

    for td, budget in targets:
        print("=" * 65)
        r = search_twins(td, max_seconds=budget,
                         primes_small=primes_small, inv6_small=inv6_small,
                         primes_ext=primes_ext, inv6_ext=inv6_ext,
                         pool=pool)

        if r["found"]:
            print(f"  FOUND in {r['elapsed']:.2f}s!")
            print(f"  ({r['p_head']}...{r['p_tail']}, +2)")
            print(f"  {r['digits']:,} digits")
            print(f"  Pipeline: {r['total_raw']:,} raw -> {r['surv']:,} sieved -> {r['tested']:,} tested")
            if largest is None or r["digits"] > largest["digits"]:
                largest = r
        else:
            print(f"  Timeout ({r['elapsed']:.1f}s), {r['tested']:,} tested in {r['batches']} batches")

        results.append({k: v for k, v in r.items() if k not in ("p", "q")})

    pool.close()
    pool.join()

    if largest:
        print()
        print("=" * 65)
        print("LARGEST TWIN PRIME FOUND")
        print("=" * 65)
        sp = str(largest["p"])
        print(f"  Digits: {largest['digits']:,}")
        print(f"  p = {sp[:70]}")
        if len(sp) > 140:
            print("      ...")
        print(f"      {sp[-70:]}")
        print(f"  q = p + 2")
        print(f"  Time: {largest['elapsed']:.2f}s")
        print()
        # Baselines
        prev_fast = {100: 0.03, 500: 1.48, 1000: 57.32, 2000: 303.80}
        prev_alg = {100: 0.59, 500: 16.40, 1000: 38.50, 2000: 280.34}
        prev_bench = {100: 1.05, 500: 6.98, 1000: 170.21, 2000: 912.15}
        for dk in sorted(prev_fast.keys()):
            for res in results:
                if res.get("found") and res.get("digits", 0) >= dk:
                    t = res["elapsed"]
                    parts = [f"{dk}-digit: {t:.2f}s"]
                    if dk in prev_fast:
                        parts.append(f"vs fast_hunter {prev_fast[dk]/t:.1f}x")
                    if dk in prev_alg:
                        parts.append(f"vs algebraic {prev_alg[dk]/t:.1f}x")
                    if dk in prev_bench:
                        parts.append(f"vs bench_baseline {prev_bench[dk]/t:.1f}x")
                    print(f"  {', '.join(parts)}")
                    break

    out = {"results": results}
    if largest:
        out["largest_p"] = str(largest["p"])
        out["largest_q"] = str(largest["q"])
        out["largest_digits"] = largest["digits"]

    with open("twin_prime_engine_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print()
    print("Saved to twin_prime_engine_results.json")


if __name__ == "__main__":
    main()
