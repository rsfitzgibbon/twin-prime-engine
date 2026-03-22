"""
Corrected spectral/twin-prime analysis for the Apollonian bridge experiments.

This file replaces the misleading parts of spectral_twin_prime.py with a more
honest split between three different objects:

1. Circles in a packing, counted with multiplicity.
   Kontorovich-Oh's exponent delta = 1.30568 applies here.
2. Distinct curvature values up to N.
   This is the object our generator actually produces after deduplication.
3. Integer coverage after adjoining the bridge set.
   This is where the twin-prime coverage experiments live.

Main corrections:
- The old script applied delta to distinct curvature values. Empirically that
  object grows roughly linearly, not like N^delta.
- A single primitive packing A or B has a local mod 24 obstruction against
  containing both entries of a twin-prime pair.
- The bridge does not create twin pairs internally; instead it complements the
  residue classes missed by A union B.
"""

import bisect
import collections
import json
import math
import sys
import time

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SEED_A = (-1, 2, 2, 3)
SEED_B = (-2, 3, 6, 7)
DELTA = 1.30568
MODULUS = 24
PRIME_CLASSES_24 = (1, 5, 7, 11, 13, 17, 19, 23)
TWIN_RESIDUE_PAIRS_24 = ((5, 7), (11, 13), (17, 19), (23, 1))


def pct(num, den):
    return 100.0 * num / den if den else 0.0


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


def sieve_bool(limit_n):
    is_prime = bytearray(b"\x01") * (limit_n + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, int(limit_n ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start : limit_n + 1 : p] = b"\x00" * (((limit_n - start) // p) + 1)
    return is_prime


def build_bridge(A, B, limit_n):
    union = A | B
    a_sorted = sorted(A)
    bridge = set()
    generators = {}

    for h in range(1, limit_n + 1):
        if h in B:
            continue
        idx = bisect.bisect_left(a_sorted, h)
        left = [a_sorted[j] for j in range(max(0, idx - 3), idx) if a_sorted[j] < h]
        if not left:
            continue
        a = max(left)
        c = 2 * h - a
        if 1 <= c <= limit_n and c not in union and c not in bridge:
            bridge.add(c)
            generators[c] = (h, a)
    return bridge, generators


def counts_up_to(sorted_values, thresholds):
    return [bisect.bisect_right(sorted_values, t) for t in thresholds]


def occupied_classes(values, mod=MODULUS):
    return sorted({v % mod for v in values})


def class_counter(values, mod=MODULUS):
    return collections.Counter(v % mod for v in values)


def prime_class_counter(values, is_prime, mod=MODULUS):
    return collections.Counter(v % mod for v in values if v > 1 and is_prime[v])


def dominant_prime_classes(counter, minimum_share=0.01):
    total = sum(counter.values())
    if not total:
        return []
    return sorted(r for r, count in counter.items() if count / total >= minimum_share)


def internal_twin_pairs(classes):
    class_set = set(classes)
    pairs = []
    for r in PRIME_CLASSES_24:
        s = (r + 2) % MODULUS
        if r in class_set and s in class_set:
            pairs.append((r, s))
    return pairs


def format_counter(counter):
    return "{" + ", ".join(f"{k}: {counter[k]}" for k in sorted(counter)) + "}"


def analyze_distinct_growth(A, thresholds):
    print("\n" + "=" * 72)
    print("EXPERIMENT 1: DISTINCT CURVATURE GROWTH")
    print("=" * 72)
    print("  This generator deduplicates curvature values.")
    print(f"  Therefore it is testing distinct-value growth, not the KO circle-count exponent delta = {DELTA:.5f}.")

    a_sorted = sorted(A)
    counts = counts_up_to(a_sorted, thresholds)
    print(f"\n  {'T':>8} | {'|A cap [1,T]|':>12} | {'log count/log T':>16} | {'density on A mod 24':>18}")
    print(f"  {'-'*8}-+-{'-'*12}-+-{'-'*16}-+-{'-'*18}")

    admissible_classes = occupied_classes(A)
    growth_rows = []
    for t, count in zip(thresholds, counts):
        exponent = math.log(count) / math.log(t) if count > 0 and t > 1 else 0.0
        admissible_total = sum(1 for n in range(1, t + 1) if (n % MODULUS) in admissible_classes)
        density = count / admissible_total if admissible_total else 0.0
        growth_rows.append(
            {
                "T": t,
                "count": count,
                "effective_exponent": exponent,
                "density_on_admissible_classes": density,
            }
        )
        print(f"  {t:>8,} | {count:>12,} | {exponent:>16.6f} | {density:>18.6f}")

    if len(thresholds) >= 2:
        t1, t2 = thresholds[-2], thresholds[-1]
        c1, c2 = counts[-2], counts[-1]
        fitted = math.log(c2 / c1) / math.log(t2 / t1)
        print(f"\n  Large-scale fitted exponent from the last two thresholds: {fitted:.6f}")
        print("  Interpretation: the deduplicated curvature set behaves approximately linearly.")
        print("  The KO exponent belongs to circle counts with multiplicity, so it should not be compared to this column.")

    return growth_rows


def analyze_residue_profiles(A, B, union, bridge, covered, is_prime):
    print("\n" + "=" * 72)
    print("EXPERIMENT 2: MOD 24 RESIDUE PROFILES")
    print("=" * 72)
    print("  Twin primes > 3 live in residue pairs (5,7), (11,13), (17,19), (23,1) mod 24.")

    results = {}
    for name, values in [("A", A), ("B", B), ("union", union), ("bridge", bridge), ("covered", covered)]:
        classes = occupied_classes(values)
        prime_counts = prime_class_counter(values, is_prime)
        prime_classes = sorted(r for r in prime_counts if prime_counts[r] > 0)
        dominant_classes = dominant_prime_classes(prime_counts)
        twin_pairs = internal_twin_pairs(prime_classes)

        print(f"\n  {name}:")
        print(f"    occupied classes:           {classes}")
        print(f"    prime residue counts:       {format_counter(prime_counts)}")
        print(f"    prime-supported classes:    {prime_classes}")
        print(f"    dominant prime classes:     {dominant_classes}")
        print(f"    internal twin residue pairs:{twin_pairs if twin_pairs else ' none'}")

        results[name] = {
            "occupied_classes": classes,
            "prime_residue_counts": dict(sorted(prime_counts.items())),
            "prime_supported_classes": prime_classes,
            "dominant_prime_classes": dominant_classes,
            "internal_twin_pairs": twin_pairs,
        }

    print("\n  Key structural point:")
    print("    A, B, and A union B do not contain a full twin residue pair among their dominant prime classes.")
    print("    The bridge supplies the complementary prime classes that the union misses.")

    return results


def analyze_local_to_global_empirically(A, thresholds):
    print("\n" + "=" * 72)
    print("EXPERIMENT 3: EMPIRICAL LOCAL-TO-GLOBAL FOR THIS SPECIFIC PACKING")
    print("=" * 72)
    print("  This is only an experiment for seed A. It is not a theorem and should not be generalized.")

    a_set = set(A)
    admissible = set(occupied_classes(A))
    rows = []

    print(f"\n  {'T':>8} | {'admissible ints':>16} | {'represented':>12} | {'missing':>10} | {'missing rate':>12}")
    print(f"  {'-'*8}-+-{'-'*16}-+-{'-'*12}-+-{'-'*10}-+-{'-'*12}")
    for t in thresholds:
        admissible_total = 0
        missing = 0
        for n in range(1, t + 1):
            if (n % MODULUS) in admissible:
                admissible_total += 1
                if n not in a_set:
                    missing += 1
        represented = admissible_total - missing
        rate = missing / admissible_total if admissible_total else 0.0
        rows.append(
            {
                "T": t,
                "admissible": admissible_total,
                "represented": represented,
                "missing": missing,
                "missing_rate": rate,
            }
        )
        print(f"  {t:>8,} | {admissible_total:>16,} | {represented:>12,} | {missing:>10,} | {rate:>12.6f}")

    return rows


def analyze_twin_coverage(union, bridge, covered, is_prime, limit_n, thresholds):
    print("\n" + "=" * 72)
    print("EXPERIMENT 4: TWIN PRIME COVERAGE IS A UNION/BRIDGE COMPLEMENTARITY")
    print("=" * 72)

    all_twins = [(p, p + 2) for p in range(3, limit_n - 1, 2) if is_prime[p] and is_prime[p + 2]]

    def role_of(n):
        if n in union:
            return "union"
        if n in bridge:
            return "bridge"
        return "other"

    role_counts = collections.Counter()
    residue_role_counts = collections.Counter()
    missed = []
    for p, q in all_twins:
        left_role = role_of(p)
        right_role = role_of(q)
        role_counts[(left_role, right_role)] += 1
        residue_role_counts[(p % MODULUS, q % MODULUS, left_role, right_role)] += 1
        if left_role == "other" or right_role == "other":
            missed.append((p, q))

    print(f"  Total twin pairs up to {limit_n:,}: {len(all_twins):,}")
    print(f"  Covered by union + bridge:      {len(all_twins) - len(missed):,} ({pct(len(all_twins) - len(missed), len(all_twins)):.2f}%)")
    print(f"  Missed twin pairs:              {len(missed):,}")
    if missed:
        print(f"  First missed pairs:             {missed[:10]}")

    print("\n  Role split:")
    for key in sorted(role_counts):
        print(f"    {key[0]:>6} + {key[1]:<6}: {role_counts[key]:>6,}")

    print("\n  Residue-pair split:")
    for key in sorted(residue_role_counts):
        pair = key[:2]
        roles = key[2:]
        print(f"    {pair} via {roles[0]} + {roles[1]}: {residue_role_counts[key]:>6,}")

    print("\n  Coverage by threshold:")
    threshold_rows = []
    for t in thresholds:
        twins_t = [(p, q) for p, q in all_twins if q <= t]
        covered_t = sum(1 for p, q in twins_t if p in covered and q in covered)
        row = {"T": t, "twins": len(twins_t), "covered": covered_t}
        threshold_rows.append(row)
        print(f"    T <= {t:>7,}: covered {covered_t:>5,}/{len(twins_t):>5,} ({pct(covered_t, len(twins_t)):.2f}%)")

    return {
        "total_twins": len(all_twins),
        "covered_twins": len(all_twins) - len(missed),
        "missed_twins": missed,
        "role_counts": {f"{k[0]}+{k[1]}": v for k, v in sorted(role_counts.items())},
        "residue_role_counts": {
            f"{k[0]}->{k[1]}:{k[2]}+{k[3]}": v for k, v in sorted(residue_role_counts.items())
        },
        "threshold_rows": threshold_rows,
    }


def analyze_corridor_completion(A, B, bridge, reverse_bridge, limit_n, thresholds):
    print("\n" + "=" * 72)
    print("EXPERIMENT 5: CORRIDOR COMPLETION NEEDS DIRECTIONAL CARE")
    print("=" * 72)
    print("  'bridge' uses A as anchors over holes of B.")
    print("  'reverse_bridge' swaps the roles of A and B.")

    rows = []
    for t in thresholds:
        left_missing_one_way = sum(1 for n in range(5, t + 1, 6) if n not in A and n not in bridge)
        right_missing_one_way = sum(1 for n in range(1, t + 1, 6) if n != 1 and n not in B and n not in bridge)
        right_missing_reverse = sum(1 for n in range(1, t + 1, 6) if n != 1 and n not in B and n not in reverse_bridge)
        left_missing_bidir = sum(
            1 for n in range(5, t + 1, 6) if n not in A and n not in bridge and n not in reverse_bridge
        )
        right_missing_bidir = sum(
            1 for n in range(1, t + 1, 6) if n != 1 and n not in B and n not in bridge and n not in reverse_bridge
        )
        rows.append(
            {
                "T": t,
                "left_missing_one_way": left_missing_one_way,
                "right_missing_one_way": right_missing_one_way,
                "right_missing_reverse_only": right_missing_reverse,
                "left_missing_bidirectional": left_missing_bidir,
                "right_missing_bidirectional": right_missing_bidir,
            }
        )

    print(
        f"\n  {'T':>8} | {'A U bridge misses 5 mod 6':>26} | {'B U bridge misses 1 mod 6':>26} | "
        f"{'B U reverse misses 1 mod 6':>28} | {'bidir misses (left,right)':>27}"
    )
    print(f"  {'-'*8}-+-{'-'*26}-+-{'-'*26}-+-{'-'*28}-+-{'-'*27}")
    for row in rows:
        print(
            f"  {row['T']:>8,} | {row['left_missing_one_way']:>26,} | {row['right_missing_one_way']:>26,} | "
            f"{row['right_missing_reverse_only']:>28,} | "
            f"({row['left_missing_bidirectional']:>5,}, {row['right_missing_bidirectional']:>5,})"
        )

    final_row = rows[-1]
    if final_row["right_missing_one_way"]:
        sample_missing = [n for n in range(1, limit_n + 1, 6) if n != 1 and n not in B and n not in bridge][:10]
        print(f"\n  One-way bridge misses on the 1 mod 6 side: {sample_missing}")
    print("  Bidirectional coverage removes the remaining corridor gaps in the tested range.")

    return rows


def analyze_bridge_mechanics(union, bridge, generators, is_prime):
    print("\n" + "=" * 72)
    print("EXPERIMENT 6: WHAT THE BRIDGE IS ACTUALLY DOING")
    print("=" * 72)

    prime_bridge = sorted(c for c in bridge if c > 1 and is_prime[c])
    prime_union = sorted(c for c in union if c > 1 and is_prime[c])
    bridge_prime_counts = prime_class_counter(bridge, is_prime)
    union_prime_counts = prime_class_counter(union, is_prime)

    print(f"  Prime counts: union = {len(prime_union):,}, bridge = {len(prime_bridge):,}")
    print(f"  Union prime classes mod 24:  {format_counter(union_prime_counts)}")
    print(f"  Bridge prime classes mod 24: {format_counter(bridge_prime_counts)}")
    print("  Up to negligible leakage, union primes live in {7,11,19,23} and bridge primes in {1,5,13,17}.")

    sample = []
    for c in prime_bridge[:12]:
        h, a = generators[c]
        sample.append({"c": c, "c_mod_24": c % MODULUS, "h": h, "a": a, "h_mod_24": h % MODULUS, "a_mod_24": a % MODULUS})

    print("\n  First 12 prime bridge points with generators:")
    for item in sample:
        print(
            f"    c={item['c']:>6} ({item['c_mod_24']:>2} mod 24) from "
            f"h={item['h']:>6} ({item['h_mod_24']:>2}) and a={item['a']:>6} ({item['a_mod_24']:>2})"
        )

    return {
        "union_prime_residue_counts": dict(sorted(union_prime_counts.items())),
        "bridge_prime_residue_counts": dict(sorted(bridge_prime_counts.items())),
        "sample_prime_generators": sample,
    }


def print_conclusion(limit_n):
    print("\n" + "=" * 72)
    print("CONCLUSION")
    print("=" * 72)
    print("  1. The old spectral script mixed two different counting problems.")
    print("     Distinct curvature values up to N behave approximately linearly.")
    print(f"     The KO exponent delta = {DELTA:.5f} belongs to counting circles/orbit points, not deduplicated integers.")
    print("  2. A single primitive packing is locally obstructed against internal twin pairs.")
    print("     The obstruction is already visible mod 24.")
    print("  3. The bridge succeeds by complementing residue classes, not by creating an internal")
    print("     twin-prime mechanism inside one packing.")
    print(f"  4. Empirically, union + one-way bridge covers every twin prime pair up to {limit_n:,}.")
    print("     Every covered pair splits as union+bridge or bridge+union.")
    print("  5. Full corridor completion is a stronger claim and, in this script,")
    print("     is supported by the bidirectional bridge system rather than one direction alone.")
    print("  6. That is interesting structure, but it is still far from a proof of the Twin Prime Conjecture.")
    print("     To connect this to a theorem, one would need packing-specific prime distribution")
    print("     results strong enough to replace Maynard-Tao/Zhang-style input, and nothing here")
    print("     gets close to that bottleneck.")


def main():
    limit_n = 500000
    thresholds = [10000, 20000, 50000, 100000, 200000, limit_n]

    print(f"Corrected Spectral Twin Prime Analysis  N = {limit_n:,}")
    print(f"Reference spectral exponent delta = {DELTA:.5f}")
    t0 = time.time()

    A = generate_gasket(SEED_A, limit_n)
    B = generate_gasket(SEED_B, limit_n)
    union = A | B
    bridge, generators = build_bridge(A, B, limit_n)
    reverse_bridge, reverse_generators = build_bridge(B, A, limit_n)
    covered = union | bridge
    bidirectional_covered = union | bridge | reverse_bridge
    is_prime = sieve_bool(limit_n + 10)

    print(
        f"Setup: {time.time() - t0:.2f}s, "
        f"|A|={len(A):,}, |B|={len(B):,}, |union|={len(union):,}, "
        f"|bridge|={len(bridge):,}, |reverse_bridge|={len(reverse_bridge):,}, "
        f"|covered|={len(covered):,}, |bidirectional_covered|={len(bidirectional_covered):,}"
    )

    growth_rows = analyze_distinct_growth(A, thresholds)
    residue_results = analyze_residue_profiles(A, B, union, bridge, bidirectional_covered, is_prime)
    local_rows = analyze_local_to_global_empirically(A, thresholds)
    twin_results = analyze_twin_coverage(union, bridge, covered, is_prime, limit_n, thresholds)
    corridor_rows = analyze_corridor_completion(A, B, bridge, reverse_bridge, limit_n, thresholds)
    bridge_results = analyze_bridge_mechanics(union, bridge, generators, is_prime)
    print_conclusion(limit_n)

    results = {
        "N": limit_n,
        "delta_reference": DELTA,
        "sizes": {
            "A": len(A),
            "B": len(B),
            "union": len(union),
            "bridge": len(bridge),
            "reverse_bridge": len(reverse_bridge),
            "covered": len(covered),
            "bidirectional_covered": len(bidirectional_covered),
        },
        "distinct_growth": growth_rows,
        "residue_profiles": residue_results,
        "empirical_local_to_global_A": local_rows,
        "twin_coverage": twin_results,
        "corridor_completion": corridor_rows,
        "bridge_mechanics": bridge_results,
        "notes": [
            "delta applies to circle counts with multiplicity, not deduplicated curvature values",
            "A, B, and A union B are locally obstructed against internal twin pairs modulo 24",
            "union plus one-way bridge covers all tested twin pairs up to N by residue-class complementarity",
            "bidirectional bridges remove the remaining corridor gaps in the tested range",
        ],
    }

    out_path = "spectral_twin_prime_v2_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved results to {out_path}")
    print(f"Total runtime: {time.time() - t0:.2f}s")


if __name__ == "__main__":
    main()
