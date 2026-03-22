# Twin Prime Search Engine

High-performance twin prime finder. The **Rust v5 engine** discovers **2000-digit twin primes in ~66 seconds** — up to **192× faster** than the Python v2 engine.

**Author:** Twin Prime Engine Project
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

- `Newtheory`

That branch owns exploratory theory work and prototypes until they are promoted here.

---

## Performance

### Rust v5 (GMP + In-place Ops + Bitset Sieve)

| Digits | Time | vs Python v2 | vs Rust v3 | Pipeline |
|--------|---------|--------------|------------|----------|
| 100 | **0.01s** | 192× faster | 2.6× faster | 6.1M raw -> 81K sieved -> 252 tested |
| 500 | **1.06s** | 11.3× faster | 2.0× faster | 8.2M raw -> 61K sieved -> 78 tested |
| 1,000 | **13.30s** | 5.7× faster | 1.0× faster | 61M raw -> 452K sieved -> 8K tested |
| 2,000 | **65.96s** | 23.6× faster | 9.6× faster | 82M raw -> 604K sieved -> 5.3K tested |

### v5 vs Previous Rust Versions

| Digits | v3 (num-bigint) | v4 (GMP) | v5 (optimized) | v5 speedup |
|--------|-----------------|----------|----------------|------------|
| 100 | 0.03s | ~0.03s | **0.01s** | 3× |
| 500 | 2.13s | ~1.0s | **1.06s** | 1-2× |
| 1,000 | 13.74s | ~10s | **13.30s** | ~1× |
| 2,000 | 631s | 882s | **65.96s** | **9.6-13.4×** |

### Python v2 (Reference)

| Digits | Time | Sieved Candidates |
|--------|---------|-------------------|
| 100 | 2.23s | 36,817 |
| 500 | 11.98s | 36,790 |
| 1,000 | 75.28s | 36,813 |
| 2,000 | 1,554s | 36,880 |

## Architecture

Both engines use a multi-stage pipeline for twin primes of the form `(6m-1, 6m+1)`:

### 1. Two-Tier Algebraic Sieve (10^8 primes)

- **Base sieve** (primes to 10^6): Per-prime loop eliminates candidates where `6m ± 1 ≡ 0 (mod p)`
- **Extended sieve** (primes 10^6 to 10^8): Eliminates remaining composites
- Survival rate: ~1.26% of candidates pass the sieve

### 2. Multi-Stage Primality Testing

- **Multi-base SPRP** (bases 2, 3, 5, 7) — fast composite filter, eliminates ~99.999% of remaining candidates
- **BPPSW confirmation** (SPRP(2) + Lucas test) — deterministic for all known composites

### 3. Parallel Execution

- **Rust v5:** Rayon work-stealing across all CPU cores, non-overlapping batch placement
- **Python v2:** ThreadPool x12 (gmpy2 releases GIL)

### 4. Adaptive Parameters

- Batch size and sieve depth scale with digit target
- Low digits (<=150): base sieve only, 512K batch, 62KB bitset
- High digits (>1500): full sieve, 10-20M batch, 1.2-2.5MB bitset

## v5 Optimizations

| Optimization | Impact | Details |
|-------------|--------|---------|
| **Packed bitset sieve** | 8x smaller working set | Bool array -> u64 bitset; fits L2 cache at 1000 digits |
| **In-place SPRP** | Zero per-test allocation | Pre-allocated SprpCtx/TestCtx buffers reused across all candidates |
| **GMP Montgomery modpow** | Hardware-accelerated | rug crate wraps GMP's hand-tuned assembly for modular exponentiation |
| **BPPSW-only confirmation** | Fewer redundant tests | is_probably_prime(0) vs (25): same accuracy, no extra Miller-Rabin rounds |
| **Non-overlapping batches** | No wasted coverage | Sequential stride from random base eliminates batch overlap |

## Why Rust v5 Is Faster

| Factor | Python v2 | Rust v5 |
|--------|-----------|---------|
| Sieve data structure | bool array (10MB) | Packed bitset (1.25MB, fits L2 cache) |
| Sieve loop | ~500ns/iteration (interpreter) | ~1ns/iteration (native + cache-friendly) |
| Parallelism | ThreadPool (GIL contention) | Rayon (true shared-memory, non-overlapping) |
| BigInt | gmpy2 -> C FFI | GMP via rug (in-place ops, zero allocation) |
| Primality testing | Allocates per-test | Pre-allocated buffers, in-place pow_mod_mut |
| Sieve build | 5-11s | 2.3s |
| Memory per worker | ~10MB (bool array + objects) | ~1.3MB (bitset + 7 Integer buffers) |

The **192x speedup at 100 digits** reflects sieve dominance where the bitset cache advantage is largest. The **13.4x speedup at 2000 digits** (v5 vs v4) demonstrates that in-place GMP operations eliminate the allocation bottleneck that plagued v4, where millions of temporary Integer objects caused cache thrashing and memory pressure.

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
```

This path is aimed at candidate generation, ranking, and finite-range benchmarking,
not the large-digit search engine pipeline above.

## Key Insight: Selberg Integration

The companion benchmark (`twin_prime_gasket_test.py`) tested three strategies:

| Strategy | Approach | Result |
|----------|----------|--------|
| **Baseline** | Sieve 5×10^7 + parallel SPRP | Reference timing |
| **Selberg-scored** | +200 extra primes scoring | 15× fewer SPRP tests, but scoring overhead negates benefit |
| **Gasket-aligned** | Residue filter mod 2520 | Covers 2520/2520 = 100% (universal covering, no-op) |

**Solution:** Integrate Selberg's extra trial division directly INTO the sieve (extending from 5×10^7 to 10^8). This captures the same benefit with zero per-candidate overhead.

## The Gapless Gasket

**"The Gapless Gasket: Universal Residue Coverage via Apollonian Packing Pairs and Bridge Complements"** — Twin Prime Engine Project, March 2026

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
| `rust-engine/` | Rust v5 engine (GMP + in-place SPRP + bitset sieve + Rayon) |
| `twin_prime_engine.py` | Python v2 engine (gmpy2 + NumPy) |
| `twin_prime_gasket_test.py` | Benchmark: Baseline vs Selberg vs Gasket |
| `spectral_twin_prime_v2.py` | Apollonian gasket generation & coverage analysis |
| `twin_prime_finder.py` | Practical finite-range corridor + arithmetic finder |
| `rust_twin_prime_finder/` | Rust acceleration for the practical finder path |
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
  author={Twin Prime Engine Project},
  year={2026}
}
```

## License

MIT
