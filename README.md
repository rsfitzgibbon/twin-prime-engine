# Twin Prime Search Engine

High-performance twin prime finder. The **Rust v3 engine** discovers **1000-digit twin primes in ~14 seconds** — up to **76× faster** than the Python v2 engine.

**Author:** Twin Prime Engine Project
**License:** MIT

---

## Performance

### Rust v3 (Native)

| Digits | Time | vs Python v2 | Pipeline |
|--------|---------|--------------|----------|
| 100 | **0.03s** | 76× faster | 2.5M raw → 32K sieved → 112 tested |
| 500 | **2.13s** | 5.6× faster | 12M raw → 88K sieved → 165 tested |
| 1,000 | **13.74s** | 5.5× faster | 40M raw → 294K sieved → 2K tested |
| 2,000 | **631s** | 2.5× faster | 80M raw → 588K sieved → 15K tested |

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

- **Rust v3:** Rayon work-stealing across all CPU cores, zero overhead shared memory
- **Python v2:** ThreadPool ×12 (gmpy2 releases GIL)

### 4. Adaptive Parameters

- Batch size and sieve depth scale with digit target
- Low digits (≤150): base sieve only, 500K batch
- High digits (>1500): full sieve, 10M batch

## Why Rust Is Faster

| Factor | Python v2 | Rust v3 |
|--------|-----------|---------|
| Sieve loop | ~500ns/iteration (interpreter) | ~1ns/iteration (native) |
| Parallelism | ThreadPool (GIL contention) | Rayon (true shared-memory threads) |
| BigInt | gmpy2 → C FFI calls | num-bigint (pure Rust, no FFI) |
| Sieve build | 5–11s | 1.88s |
| Memory | Per-worker pickle overhead | Zero-copy shared data |

The 76× speedup at 100 digits reflects sieve dominance (tight loops). At 2000+ digits, the primality test dominates and both engines call equivalent algorithms, narrowing the gap to 2.5×.

## Installation

### Rust Engine (Recommended)

```bash
cd rust-engine
cargo build --release
./target/release/twin_prime_engine
```

Requires Rust 1.70+ and MSVC Build Tools (Windows) or GCC (Linux/macOS).

### Python Engine

```bash
pip install gmpy2 numpy
python twin_prime_engine.py
```

Requires Python 3.10+, gmpy2, NumPy.

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
| `rust-engine/` | Rust v3 engine (Rayon + num-bigint + BPPSW) |
| `twin_prime_engine.py` | Python v2 engine (gmpy2 + NumPy) |
| `twin_prime_gasket_test.py` | Benchmark: Baseline vs Selberg vs Gasket |
| `spectral_twin_prime_v2.py` | Apollonian gasket generation & coverage analysis |
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
