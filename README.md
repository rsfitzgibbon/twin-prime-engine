# Twin Prime Search Engine

High-performance twin prime finder. The **Rust v5.2 engine** discovers **1000-digit twin primes in ~2 seconds** and **2000-digit twin primes in ~153 seconds** — up to **339× faster** than the Python v2 engine.

**License:** MIT

---

## Repository Roles

This repo is the **curated release-facing engine branch**.

Use this repo for:

- canonical engine code
- public benchmark claims
- release notes
- the canonical paper copy

The active experimental branch lives at:

- `Newtheory` (separate local research branch)

That branch owns exploratory theory work and prototypes until they are promoted here.

---

## Performance

### Rust v5.2 (GMP + 3-Tier Sieve + Precomputed Mod)

| Digits | Time | vs Python v2 | vs Rust v3 | Pipeline |
|--------|---------|--------------|------------|----------|
| 100 | **0.01s** | 339× faster | 4.6× faster | 3.6M raw -> 47K sieved -> 176 tested |
| 500 | **0.30s** | 40× faster | 7.1× faster | 14.3M raw -> 106K sieved -> 608 tested |
| 1,000 | **2.09s** | 36× faster | 6.6× faster | 25.6M raw -> 189K sieved -> 552 tested |
| 2,000 | **~153s** | 10× faster | 4.1× faster | 10.2M raw -> 70K sieved -> 1.6K tested |

**v5.2 improvements over v5.1:** Precomputed `m_low mod p` for all 11M sieve primes eliminates GMP `mod_u` from the sieve hot path. Per-batch sieve uses only u64 arithmetic. **3.6× faster at 500d, 5.3× faster at 1000d.** At 2000d, SPRP testing dominates so sieve improvements have diminishing returns.

### v5.2 vs Previous Rust Versions

| Digits | v3 (num-bigint) | v4 (GMP) | v5.2 (optimized) | v5.2 speedup |
|--------|-----------------|----------|------------------|--------------|
| 100 | 0.03s | ~0.03s | **0.01s** | 3-5× |
| 500 | 2.13s | ~1.0s | **0.30s** | **3-7×** |
| 1,000 | 13.74s | ~10s | **2.09s** | **5-7×** |
| 2,000 | 631s | 882s | **~153s** | **4.1-5.8×** |

### Python v2 (Reference)

| Digits | Time | Sieved Candidates |
|--------|---------|-------------------|
| 100 | 2.23s | 36,817 |
| 500 | 11.98s | 36,790 |
| 1,000 | 75.28s | 36,813 |
| 2,000 | 1,554s | 36,880 |

## Architecture

Both engines use a multi-stage pipeline for twin primes of the form `(6m-1, 6m+1)`:

### 1. Three-Tier Algebraic Sieve (2×10^8 primes)

- **Base sieve** (primes to 10^6): Per-prime loop eliminates candidates where `6m ± 1 ≡ 0 (mod p)`
- **Extended sieve** (primes 10^6 to 10^8): Eliminates remaining composites
- **Deep sieve** (primes 10^8 to 2×10^8): Additional 5.3M primes for ≥1500-digit targets
- Survival rate: ~0.68% of candidates pass the deep sieve (vs ~0.74% without deep tier)

### 2. Multi-Stage Primality Testing

- **Multi-base SPRP** (bases 2, 3, 5, 7) — fast composite filter, eliminates ~99.999% of remaining candidates
- **BPPSW confirmation** (SPRP(2) + Lucas test) — deterministic for all known composites

### 3. Parallel Execution

- **Rust v5.2:** Rayon work-stealing across all CPU cores, non-overlapping batch placement
- **Python v2:** ThreadPool x12 (gmpy2 releases GIL)

### 4. Adaptive Parameters

- Batch size and sieve depth scale with digit target
- Low digits (<=150): base sieve only, 512K batch, 62KB bitset
- High digits (>1500): full sieve, 10-20M batch, 1.2-2.5MB bitset

## v5.2 Optimizations

| Optimization | Impact | Details |
|-------------|--------|---------|
| **Precomputed mod_u** | **3.6-5.3× at 500-1000d** | One-time GMP `mod_u` per search; per-batch sieve uses u64 arithmetic only |
| **Packed bitset sieve** | 8x smaller working set | Bool array -> u64 bitset; fits L2 cache at 1000 digits |
| **In-place SPRP** | Zero per-test allocation | Pre-allocated SprpCtx/TestCtx buffers reused across all candidates |
| **GMP Montgomery modpow** | Hardware-accelerated | rug crate wraps GMP's hand-tuned assembly for modular exponentiation |
| **BPPSW-only confirmation** | Fewer redundant tests | is_probably_prime(0) vs (25): same accuracy, no extra Miller-Rabin rounds |
| **Non-overlapping batches** | No wasted coverage | Sequential stride from random base eliminates batch overlap |
| **3-tier deep sieve** | 8% fewer survivors | Deep tier (10^8-2×10^8) adds 5.3M primes for ≥1500-digit targets |
| **u32 compact cache** | 50% less cache memory | 84.5 MB vs ~169 MB; sieve build 35% faster (2.9s vs 4.5s) |
| **Inline survivor testing** | Zero Vec allocation | Iterate bitset directly with bit tricks + hardware popcount |

## Why Rust v5.2 Is Faster

| Factor | Python v2 | Rust v5.2 |
|--------|-----------|-----------|
| Sieve data structure | bool array (10MB) | Packed bitset (1.25MB, fits L2 cache) |
| Sieve hot path | GMP mod_u per batch per prime | Precomputed mod + u64 arithmetic only |
| Sieve depth | ~5×10^7 primes | 2×10^8 primes (3 tiers, adaptive) |
| Sieve loop | ~500ns/iteration (interpreter) | ~1ns/iteration (native + cache-friendly) |
| Parallelism | ThreadPool (GIL contention) | Rayon (true shared-memory, non-overlapping) |
| BigInt | gmpy2 -> C FFI | GMP via rug (in-place ops, zero allocation) |
| Primality testing | Allocates per-test | Pre-allocated buffers, in-place pow_mod_mut |
| Sieve build | 5-11s | 2.9s (bitset Eratosthenes + u32 cache) |
| Prime cache | N/A | 84.5 MB (u32 compact tables) |
| Memory per worker | ~10MB (bool array + objects) | ~1.3MB (bitset + 7 Integer buffers) |

The **339x speedup at 100 digits** reflects sieve dominance where the bitset cache advantage is largest. The **5.3x speedup at 1000 digits** (v5.2 vs v5.1) demonstrates that eliminating GMP `mod_u` from the sieve hot path removes the last remaining per-batch overhead. At 2000+ digits, SPRP modpow cost dominates — GPU acceleration (CGBN) is the path to the next breakthrough.

## Installation

### Rust Engine (Recommended)

```bash
cd rust-engine
cargo build --release
./target/release/twin_prime_engine
```

Requires Rust 1.70+ and MSYS2/MinGW (Windows) or GCC (Linux/macOS).
On Windows, GMP is provided via MSYS2: `pacman -S mingw-w64-x86_64-gmp`
Also ensure `C:\msys64\mingw64\bin` is on `PATH` when building the GNU target,
or build from an MSYS2 MinGW shell so tools like `dlltool.exe` are visible.

Current repo configuration targets GNU Windows by default via `.cargo/config.toml`.
If your Rust installation only has MSVC, install the GNU target first:

```bash
rustup target add x86_64-pc-windows-gnu
```

### Python Engine

```bash
pip install gmpy2 numpy
python twin_prime_engine.py
```

Requires Python 3.10+, gmpy2, NumPy.

### Practical Finder

The repo also carries the promoted finite-range finder path:

```bash
python twin_prime_finder.py --limit 1000000 --mode high_precision --preview 20
python twin_prime_finder.py --limit 1000000 --mode high_precision --engine rust --preview 20
python twin_prime_finder.py --limit 1000000 --mode high_precision --engine rust --score --top-k 4000 --preview 20
```

This path is aimed at candidate generation, ranking, and finite-range benchmarking,
not the large-digit search engine pipeline above. The promoted Rust finder now
covers hard-mask generation, soft rerank, and sieve-score ranking; only geometry
annotations remain Python-only.

### Fixed-n k-Space Search

The repo also includes a fixed-`n` twin-prime search path for pairs of the form
`k*2^n +/- 1`, with `n` fixed and `k` varying.

```bash
python twin_prime_k2n_hunter.py --n 1000 --k-start 3 --k-batch-size 10000 --sieve-limit 50000 --max-seconds 10
fixed_n_k_search.exe --n 1000 --k-start 3 --k-batch-size 10000 --sieve-limit 50000 --max-seconds 10
```

The Rust runner adds exact validation, continuous checkpointed search, and
continuous survivor export for external large-PRP backends. That export mode is
the practical route for million-digit-scale work.

### OpenPFGW Bridge

For large fixed-`n` runs, the recommended path is now:

1. run `fixed_n_campaign.exe` with `--backend compact_export`
2. add `--post-sieve-limit` to cut survivors further with an exact second-stage sieve
3. feed the compact batches into `openpfgw_batch_runner.py`

Example:

```bash
fixed_n_campaign.exe --n 1288907 --backend compact_export --export-dir campaign_388k_exact_post200k --k-start 3 --k-batch-size 4096 --sieve-limit 50000 --post-sieve-limit 200000 --max-batches 4 --json-out campaign_388k_exact_post200k_summary.json
python openpfgw_batch_runner.py --compact-dir campaign_388k_exact_post200k --output-dir openpfgw_jobs_388k --prepare-only --json-out openpfgw_jobs_388k_summary.json
```

This reconstructs OpenPFGW-oriented `plus.abcd`, `minus.abcd`, `plus.txt`, and
`minus.txt` files from the compact survivor batches. If a worker has `pfgw.exe`
installed, the same bridge script can invoke it directly.

For a persistent worker, use:

```bash
python openpfgw_worker.py --compact-dir campaign_388k_exact_post200k --output-dir openpfgw_jobs_388k --pfgw-exe C:\path\to\pfgw64.exe --state-file worker_state.json --result-log worker_results.jsonl --poll-seconds 30
```

## Key Insight: Selberg Integration

The companion benchmark (`twin_prime_gasket_test.py`) tested three strategies:

| Strategy | Approach | Result |
|----------|----------|--------|
| **Baseline** | Sieve 5×10^7 + parallel SPRP | Reference timing |
| **Selberg-scored** | +200 extra primes scoring | 15× fewer SPRP tests, but scoring overhead negates benefit |
| **Gasket-aligned** | Residue filter mod 2520 | Covers 2520/2520 = 100% (universal covering, no-op) |

**Solution:** Integrate Selberg's extra trial division directly INTO the sieve (extending from 5×10^7 to 10^8). This captures the same benefit with zero per-candidate overhead.

## The Gapless Gasket

**"The Gapless Gasket: Universal Residue Coverage via Apollonian Packing Pairs and Bridge Complements"** — March 2026

Two primitive Apollonian circle packings with seeds (−1,2,2,3) and (−2,3,6,7), combined with a bidirectional bridge construction, produce a set covering **all 2520 residue classes modulo 2520**. Every twin-prime pair (p, p+2) up to N=500,000 is represented — **4,565/4,565 pairs, zero exceptions**.

### Twin Prime Coverage

| N | Covered | Total | Coverage |
|---|---------|-------|----------|
| 10,000 | 205 | 205 | 100% |
| 50,000 | 705 | 705 | 100% |
| 100,000 | 1,224 | 1,224 | 100% |
| 500,000 | 4,565 | 4,565 | 100% |

## Files

| File | Description |
|------|-------------|
| `rust-engine/` | Rust v5.2 engine (GMP + 3-tier sieve + precomputed mod + in-place SPRP + Rayon) |
| `twin_prime_engine.py` | Python v2 engine (gmpy2 + NumPy) |
| `twin_prime_gasket_test.py` | Benchmark: Baseline vs Selberg vs Gasket |
| `spectral_twin_prime_v2.py` | Apollonian gasket generation & coverage analysis |
| `twin_prime_finder.py` | Practical finite-range corridor + arithmetic finder |
| `rust_twin_prime_finder/` | Rust acceleration for the practical finder path |
| `twin_prime_k2n_hunter.py` | Python fixed-`n` k-space twin-prime hunter |
| `rust-engine/src/bin/fixed_n_k_search.rs` | Rust fixed-`n` k-space search, export, and checkpoint runner |
| `twin_prime_finder_rust_benchmark.py` | Python-vs-Rust benchmark harness for the finder |
| `paper/gapless_gasket.tex` | Full LaTeX research paper |
| `paper/supplementary_data.json` | All numerical evidence (JSON) |

## References

- Kontorovich, A. and Oh, H. (2011). "Apollonian circle packings and closed horospheres on hyperbolic 3-manifolds."
- Bourgain, J. and Fuchs, E. (2011). "A proof of the positive density conjecture for integer Apollonian circle packings."
- Hardy, G.H. and Littlewood, J.E. (1923). "Some problems of 'Partitio Numerorum'; III."
- Zhang, Y. (2014). "Bounded gaps between primes."
- Maynard, J. (2015). "Small gaps between primes."

## Citation

```bibtex
@article{fitzgibbon2026gapless,
  title={The Gapless Gasket: Universal Residue Coverage via Apollonian Packing Pairs and Bridge Complements},
  author={{Twin Prime Engine Project}},
  year={2026}
}
```

## License

MIT
