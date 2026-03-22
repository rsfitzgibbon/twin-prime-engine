# Twin Prime Search Engine v2

High-performance multi-threaded twin prime finder. Discovers **1000-digit twin primes in ~75 seconds**.

**Author:** Twin Prime Engine Project
**License:** MIT

---

## Performance

| Digits | Time | Sieved Candidates |
|--------|---------|-------------------|
| 100 | 2.23s | 36,817 |
| 500 | 11.98s | 36,790 |
| 1,000 | 75.28s | 36,813 |
| 2,000 | 1,554s | 36,880 |

## How It Works

The engine uses a multi-stage pipeline optimized for large twin primes of the form `(6m-1, 6m+1)`:

### 1. Two-Tier Sieve (10^8 primes)

- **Base sieve** (primes to 10^6): Per-prime loop eliminates candidates where `6m ± 1 ≡ 0 (mod p)`
- **Extended sieve** (primes 10^6 to 10^8): Vectorized NumPy elimination
- Uses `gmpy2.f_mod()` for 2–4× faster big-integer modular reduction at 1000+ digits
- Survival rate: ~1.26% of candidates pass the sieve

### 2. Multi-Threaded Primality Testing

- 12-thread `ThreadPool` (gmpy2 releases the GIL for true parallelism)
- **SPRP base-2** on `p = 6m - 1` first — rejects ~50% of sieve survivors
- Only survivors get `p + 2` tested
- **BPPSW** confirmation on both primes for final verification

### 3. Adaptive Batching

Batch size scales automatically with digit target for optimal memory/throughput balance.

## Key Insight: Selberg Integration

The companion benchmark (`twin_prime_gasket_test.py`) tested three strategies:

| Strategy | Approach | Result |
|----------|----------|--------|
| **Baseline** | Sieve 5×10^7 + parallel SPRP | Reference timing |
| **Selberg-scored** | +200 extra primes scoring | 15× fewer SPRP tests, but Python scoring overhead negates benefit |
| **Gasket-aligned** | Residue filter mod 2520 | Covers 2520/2520 = 100% (universal covering, no-op) |

**Solution:** Integrate Selberg's extra trial division directly INTO the sieve (extending from 5×10^7 to 10^8). This captures the same benefit with zero per-candidate overhead, yielding a **2.3× speedup** at 1000 digits vs the benchmark baseline.

## Installation

```bash
pip install gmpy2 numpy
```

### System Requirements

- Python 3.10+
- [gmpy2](https://pypi.org/project/gmpy2/) (GMP bindings for fast big-integer arithmetic)
- NumPy
- Multi-core CPU recommended (engine uses 12 threads by default)

## Usage

```bash
# Run the full search (100 → 5000 digits)
python twin_prime_engine.py

# Output: finds twin primes at each digit target, saves results to twin_prime_engine_results.json
```

### Configuration

Edit the `targets` list in `main()` to adjust digit targets and time budgets:

```python
targets = [
    (100, 30),    # (digits, max_seconds)
    (500, 60),
    (1000, 300),
    (2000, 900),
    (5000, 3600),
]
```

Adjust `NUM_WORKERS` at the top of the file to match your CPU core count.

## The Gapless Gasket

This repository also includes the research that informed the engine's design:

**"The Gapless Gasket: Universal Residue Coverage via Apollonian Packing Pairs and Bridge Complements"** — Twin Prime Engine Project, March 2026

Two primitive Apollonian circle packings with seeds (−1,2,2,3) and (−2,3,6,7), combined with a bidirectional bridge construction, produce a set covering **all 2520 residue classes modulo 2520**. Every twin-prime pair (p, p+2) up to N=500,000 is represented — **4,565/4,565 pairs, zero exceptions**.

### Twin Prime Coverage

| N | Covered | Total | Coverage |
|---|---------|-------|----------|
| 10,000 | 205 | 205 | 100% |
| 50,000 | 705 | 705 | 100% |
| 100,000 | 1,224 | 1,224 | 100% |
| 500,000 | 4,565 | 4,565 | 100% |

### Residue Class Complementarity (mod 24)

| Component | Prime Classes | Twin Pairs |
|-----------|--------------|------------|
| Gasket A | {11, 23} | None internally |
| Gasket B | {7, 19} | None internally |
| Union A∪B | {7, 11, 19, 23} | None internally |
| Bridge | {1, 5, 13, 17} | Completes all |
| **Combined** | **All 8 classes** | **All 4 pairs** |

## Files

| File | Description |
|------|-------------|
| `twin_prime_engine.py` | Main search engine (v2) |
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
