"""
Twin-prime search in k-space for pairs of the form k * 2^n +/- 1, with n fixed.

This shifts the sieve from literal integer ranges to the coefficient k.

Core arithmetic:
- For fixed n and prime p >= 5, k * 2^n +/- 1 is divisible by p iff
  k * 2^n == +/-1 (mod p), i.e.
  k == +/- 2^(-n) (mod p).
- Restricting to odd multiples of 3 means k = 6i + 3, so the sieve can strike
  positions i in arithmetic progressions, exactly like the corridor sieve acts
  on m for pairs (6m-1, 6m+1).

This file provides:
- a precomputed k-sieve plan for fixed n
- a validation mode against brute-force small-prime blocking
- a search loop using Proth-friendly +1 testing first, then -1 confirmation
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import gmpy2


sys.set_int_max_str_digits(1_000_000)

C2 = 0.6601618158
PROTH_BASES = (2, 3, 5, 7, 11, 13, 17, 19)


def small_prime_sieve(limit: int) -> list[int]:
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [p for p in range(2, limit + 1) if is_prime[p]]


def align_k_start(k_start: int) -> int:
    if k_start <= 3:
        return 3
    return k_start + ((3 - k_start) % 6)


def build_k_sieve_plan(n: int, sieve_primes: list[int]) -> list[tuple[int, int, tuple[int, ...]]]:
    """
    Precompute the fixed-n sieve plan.

    Each entry is:
      (p, inv6_mod_p, bad_k_residues_mod_p)
    """
    plan = []
    for p in sieve_primes:
        if p < 5:
            continue
        two_n_mod_p = pow(2, n, p)
        inv_two_n = pow(two_n_mod_p, -1, p)
        residues = tuple(sorted({inv_two_n % p, (-inv_two_n) % p}))
        inv6 = pow(6, -1, p)
        plan.append((p, inv6, residues))
    return plan


def survival_rate_from_plan(plan: list[tuple[int, int, tuple[int, ...]]]) -> float:
    rate = 1.0
    for p, _, _ in plan:
        rate *= (p - 2) / p
    return rate


def self_prime_exception_indices(n: int, p: int, k_start: int, k_count: int) -> list[int]:
    """
    Rare exact exceptions where k * 2^n +/- 1 equals the sieving prime itself.

    These matter only for very small n; for realistic record-hunting n they
    disappear because 2^n exceeds the entire sieve limit.
    """
    if n >= 63:
        return []
    two_n = 1 << n
    if two_n > p + 1:
        return []

    exceptions = []
    for numer in (p - 1, p + 1):
        if numer <= 0 or numer % two_n != 0:
            continue
        k_val = numer // two_n
        if k_val <= 0 or k_val % 6 != 3:
            continue
        if k_val < k_start:
            continue
        delta = k_val - k_start
        if delta % 6 != 0:
            continue
        idx = delta // 6
        if 0 <= idx < k_count:
            exceptions.append(idx)
    return exceptions


def sieve_k_batch(n: int, k_start: int, k_count: int, plan: list[tuple[int, int, tuple[int, ...]]]) -> bytearray:
    """
    Sieve the k values:
      k = k_start + 6*i, i = 0..k_count-1
    """
    k_start = align_k_start(k_start)
    alive = bytearray(b"\x01") * k_count

    for p, inv6, bad_residues in plan:
        k_mod_p = k_start % p
        for residue in bad_residues:
            i0 = ((residue - k_mod_p) * inv6) % p
            if i0 < k_count:
                alive[i0::p] = b"\x00" * len(alive[i0::p])

        for idx in self_prime_exception_indices(n, p, k_start, k_count):
            alive[idx] = 1

    return alive


def candidate_k_values(k_start: int, alive: bytearray) -> list[int]:
    k_start = align_k_start(k_start)
    return [k_start + 6 * idx for idx, flag in enumerate(alive) if flag]


def plus_candidate(k: int, two_n: gmpy2.mpz) -> gmpy2.mpz:
    return gmpy2.mpz(k) * two_n + 1


def minus_candidate(k: int, two_n: gmpy2.mpz) -> gmpy2.mpz:
    return gmpy2.mpz(k) * two_n - 1


def proth_or_bpsw_plus(k: int, n: int, two_n: gmpy2.mpz) -> tuple[bool, str]:
    p = plus_candidate(k, two_n)

    # Fast prefilter.
    if not gmpy2.is_strong_prp(p, 2):
        return False, "sprp2_reject"

    # Proth path when k is odd and k < 2^n. This is the usual record-search regime.
    if k % 2 == 1 and k.bit_length() <= n:
        half_pm1 = p >> 1
        for base in PROTH_BASES:
            if gmpy2.powmod(base, half_pm1, p) == p - 1:
                return True, f"proth_base_{base}"

    # Fallback preserves recall if the fixed Proth bases do not hit a witness.
    if gmpy2.is_bpsw_prp(p):
        return True, "bpsw_fallback"
    return False, "bpsw_reject"


def minus_is_probable_prime(k: int, two_n: gmpy2.mpz) -> bool:
    q = minus_candidate(k, two_n)
    if not gmpy2.is_strong_prp(q, 2):
        return False
    return bool(gmpy2.is_bpsw_prp(q))


def estimate_trials(n: int, plan: list[tuple[int, int, tuple[int, ...]]]) -> dict:
    ln_n = n * math.log(2)
    raw_trials = ln_n * ln_n / (2 * C2)
    survival = survival_rate_from_plan(plan)
    sieved_trials = raw_trials * survival
    return {
        "raw_trials": raw_trials,
        "sieve_survival": survival,
        "sieved_trials": sieved_trials,
    }


def validate_k_sieve(n: int, k_start: int, k_count: int, sieve_limit: int) -> dict:
    primes = small_prime_sieve(sieve_limit)
    plan = build_k_sieve_plan(n, primes)
    k_start = align_k_start(k_start)
    alive = sieve_k_batch(n, k_start, k_count, plan)

    if n >= 63:
        raise ValueError("validation mode is only intended for small n where brute force is tractable")

    two_n = 1 << n
    mismatches = []
    exact_survivors = 0

    for idx in range(k_count):
        k = k_start + 6 * idx
        p = k * two_n - 1
        q = p + 2
        brute_alive = True
        for prime in primes:
            if prime < 5:
                continue
            if p % prime == 0 and p != prime:
                brute_alive = False
                break
            if q % prime == 0 and q != prime:
                brute_alive = False
                break
        if brute_alive:
            exact_survivors += 1
        if bool(alive[idx]) != brute_alive:
            mismatches.append(
                {
                    "index": idx,
                    "k": k,
                    "sieve_alive": bool(alive[idx]),
                    "brute_alive": brute_alive,
                }
            )

    return {
        "n": n,
        "k_start": k_start,
        "k_count": k_count,
        "sieve_limit": sieve_limit,
        "predicted_survivors": int(sum(alive)),
        "exact_survivors": exact_survivors,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches[:20],
    }


def search_fixed_n(
    n: int,
    k_start: int,
    k_batch_size: int,
    sieve_limit: int,
    max_seconds: float,
    max_batches: int | None = None,
) -> dict:
    primes = small_prime_sieve(sieve_limit)
    plan = build_k_sieve_plan(n, primes)
    estimate = estimate_trials(n, plan)
    k_start = align_k_start(k_start)
    two_n = gmpy2.mpz(1) << n

    digits = int(n * math.log10(2)) + 1
    print(f"Fixed-n search: n={n:,} (~{digits:,} digits)")
    print(f"  k-space sieve limit: {sieve_limit:,}")
    print(f"  k batch size: {k_batch_size:,} values in the k == 3 (mod 6) progression")
    print(f"  HL raw trials/pair: {estimate['raw_trials']:,.0f}")
    print(f"  Sieve survival: {estimate['sieve_survival']:.6f}")
    print(f"  Estimated sieved trials/pair: {estimate['sieved_trials']:,.0f}")

    t0 = time.perf_counter()
    total_sieved = 0
    tested_plus = 0
    tested_minus = 0
    batches = 0

    while time.perf_counter() - t0 < max_seconds:
        if max_batches is not None and batches >= max_batches:
            break

        batch_start = k_start + batches * k_batch_size * 6
        alive = sieve_k_batch(n, batch_start, k_batch_size, plan)
        survivors = candidate_k_values(batch_start, alive)
        total_sieved += len(survivors)
        batches += 1

        for k in survivors:
            if time.perf_counter() - t0 > max_seconds:
                break

            tested_plus += 1
            plus_ok, plus_mode = proth_or_bpsw_plus(k, n, two_n)
            if not plus_ok:
                continue

            tested_minus += 1
            if minus_is_probable_prime(k, two_n):
                elapsed = time.perf_counter() - t0
                p = minus_candidate(k, two_n)
                q = p + 2
                return {
                    "found": True,
                    "k": str(k),
                    "n": n,
                    "digits": len(str(p)),
                    "p_head": str(p)[:40],
                    "p_tail": str(p)[-40:],
                    "q_head": str(q)[:40],
                    "q_tail": str(q)[-40:],
                    "elapsed": elapsed,
                    "batches": batches,
                    "total_sieved": total_sieved,
                    "tested_plus": tested_plus,
                    "tested_minus": tested_minus,
                    "plus_mode": plus_mode,
                    "k_start": str(k_start),
                    "k_batch_size": k_batch_size,
                    "sieve_limit": sieve_limit,
                }

        if batches % 5 == 0:
            elapsed = time.perf_counter() - t0
            rate = tested_plus / elapsed if elapsed else 0.0
            print(
                f"    [{elapsed:.1f}s] batches={batches:,} survivors={total_sieved:,} "
                f"plus_tests={tested_plus:,} minus_tests={tested_minus:,} ({rate:.1f} +tests/s)"
            )

    elapsed = time.perf_counter() - t0
    return {
        "found": False,
        "n": n,
        "digits": digits,
        "elapsed": elapsed,
        "batches": batches,
        "total_sieved": total_sieved,
        "tested_plus": tested_plus,
        "tested_minus": tested_minus,
        "k_start": str(k_start),
        "k_batch_size": k_batch_size,
        "sieve_limit": sieve_limit,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Search twin primes of the form k * 2^n +/- 1 with a sieve on k.")
    parser.add_argument("--n", type=int, required=True, help="Fixed exponent n.")
    parser.add_argument("--k-start", type=int, default=3, help="Starting k; will be aligned to 3 mod 6.")
    parser.add_argument("--k-batch-size", type=int, default=50_000, help="Number of k values per batch in the k == 3 mod 6 progression.")
    parser.add_argument("--sieve-limit", type=int, default=50_000, help="Small-prime sieve limit.")
    parser.add_argument("--max-seconds", type=float, default=60.0, help="Search budget.")
    parser.add_argument("--max-batches", type=int, default=None, help="Optional max batch count.")
    parser.add_argument("--validate", action="store_true", help="Validate the k-sieve against brute-force small-prime blocking.")
    parser.add_argument("--json-out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.validate:
        result = validate_k_sieve(
            n=args.n,
            k_start=args.k_start,
            k_count=args.k_batch_size,
            sieve_limit=args.sieve_limit,
        )
    else:
        result = search_fixed_n(
            n=args.n,
            k_start=args.k_start,
            k_batch_size=args.k_batch_size,
            sieve_limit=args.sieve_limit,
            max_seconds=args.max_seconds,
            max_batches=args.max_batches,
        )

    print(json.dumps(result, indent=2))
    if args.json_out is not None:
        with args.json_out.open("w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print(f"Saved JSON to {args.json_out}")


if __name__ == "__main__":
    main()
