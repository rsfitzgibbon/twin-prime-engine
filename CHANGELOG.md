# Twin Prime Engine — Development Log

**Author:** Twin Prime Engine Project
**Repository:** twin-prime-engine

---

## v5.1.0 — 2026-03-22

### Three-tier sieve + inline survivor testing

**Key changes:**
- Compact sieve-tier cache uses `u32` prime/inverse tables instead of `u64`
  - Cuts resident reusable prime-plan memory roughly in half
  - Measured cache footprint: ~84.5 MB instead of ~169 MB for the same tiers
- Bitset Eratosthenes sieve (8x less memory for prime generation: 25MB vs 200MB)
- Three-tier algebraic sieve: base (10^6) + extended (10^8) + deep (2×10^8)
  - Deep tier adds ~5.2M primes for targets ≥1500 digits
  - Reduces survivor count by ~7-10% at high digit counts
- Inline survivor testing: iterate bitset directly using bit tricks
  - Eliminates Vec<u64> allocation (saved ~5MB per batch at 2000 digits)
  - Hardware popcount for survivor counting
- Sieve depth is now adaptive: base-only for ≤150d, extended for 151-1500d, deep for ≥1501d

**Benchmarks:** *(pending — benchmark running)*

---

## v5.0.0 — 2026-03-22 (commit b0fdef6)

### In-place GMP ops + bitset sieve

**Benchmarks (12-core Rayon, Windows 11):**

| Digits | Time | vs Python v2 | vs Rust v3 | vs v4 |
|--------|------|-------------|------------|-------|
| 100 | 0.01s | 270x | 3.6x | ~3x |
| 500 | 0.88s | 13.6x | 2.4x | ~1x |
| 1,000 | 5.38s | 14.0x | 2.6x | ~1x |
| 2,000 | 66-263s | 6-24x | 2.4-9.6x | 3.4-13.4x |

Per-test speed at 2000 digits: ~15ms (4x faster than v4's ~60ms).

**Key changes:**
- Packed u64 bitset sieve (8x smaller working set, fits L2 cache)
- In-place rug::Integer operations via SprpCtx/TestCtx (zero per-test allocation)
- BPPSW-only final confirmation (is_probably_prime(0), no extra MR rounds)
- Non-overlapping sequential batch placement per round
- Deadline check inside inner testing loop (fixes 5000-digit timeout bug)

**Root cause of v4 regression:** At 2000 digits, v4 created millions of temporary
Integer objects per batch (3 per candidate construction + ~5 per SPRP call x 8 calls),
causing GMP allocator thrashing and L2/L3 cache eviction. In-place operations
eliminated this entirely.

---

## v4.0.0 — 2026-03-22 (commit 1e0640d)

### GMP-accelerated engine via rug crate

- Switched from num-bigint (pure Rust) to rug/GMP for all big-integer arithmetic
- Montgomery modpow, Toom-Cook multiplication via GMP's hand-tuned assembly
- mpz_fdiv_ui for single-limb mod in sieve phase
- Build target: x86_64-pc-windows-gnu (MSYS2/MinGW for GMP)
- Regressed at 2000 digits (882s vs v3's 631s) due to allocation overhead

---

## v3.0.0 — 2026-03-21 (commit bb6d888)

### Rust engine with num-bigint

- First Rust implementation using num-bigint (pure Rust BigUint)
- Rayon work-stealing parallelism across all CPU cores
- Two-tier algebraic sieve (base to 10^6, extended to 10^8)
- Custom SPRP(2,3,5,7) + BPPSW primality cascade
- 76x faster than Python at 100 digits, 5.5x at 1000 digits

---

## v2.0.0 — 2026-03-21 (commit e6d5c50)

### Python reference engine

- gmpy2 + NumPy implementation
- ThreadPool x12 (gmpy2 releases GIL)
- Same two-tier sieve + SPRP + BPPSW pipeline
- Baseline: 2.23s (100-digit), 75.28s (1000-digit), 1554s (2000-digit)

---

## Research Completed (2026-03-22)

### GPU Acceleration Research
- **CGBN** (NVIDIA CUDA big-integer): 21x throughput over 20-core GMP at 2048 bits
  - Supports up to 32,768 bits (covers 9,900-digit numbers)
  - Cooperative groups: 16-32 GPU threads per big number
  - Integration path: cudarc (Rust) + CGBN (CUDA C++)
- **Turkish Sieve Engine**: GPU twin prime sieve, 1.07 T-items/s on RTX 5090
  - N/6 bit structure (same as our engine), tiered prime handling
- **wgpu**: Cross-platform GPU compute for sieve marking kernel (WGSL shaders)
- **GpuOwl/PRPLL**: Found M136279841 (first GPU Mersenne prime, Oct 2024)
  - OpenCL + NTT (Number Theoretic Transform) for GPU modpow

### Physics & Number Theory Research
- **Expansion/collapse** = constructive/destructive interference of Riemann zeta-zero
  oscillations (von Mangoldt explicit formula: psi(x) as superposition of waves)
- Twin primes are "more random than primes" (Brent, ANU) — validates random search
- Hardy-Littlewood density estimate confirmed near-optimal by RMT connection
- Berry-Keating 2024: Hamiltonian with eigenvalues = Riemann zeros shown self-adjoint
- Lemke Oliver-Soundararajan 2016: consecutive primes have last-digit biases
- Chebyshev bias: primes 4k+3 outnumber 4k+1 with 99.59% logarithmic density
- SSoZ (Segmented Sieve of Zakiya): primorial 30030 wheel reduces candidates to <10%

---

## Next Steps — Optimization Roadmap

### Priority 1: Mod-30 Wheel Pre-filter (v5.1)
- **Impact:** ~40% fewer sieve iterations
- **Approach:** Skip m values where m mod 5 in {1, 4} (6m-1 or 6m+1 divisible by 5)
- **Difficulty:** Low — compress bitset to only valid m positions

### Priority 2: GPU SPRP via CGBN (v6)
- **Impact:** 12-20x throughput at 2000+ digits
- **Approach:** Batch upload survivor p1/p2 values to GPU, run SPRP(2) on all simultaneously
- **Requirements:** NVIDIA GPU + CUDA toolkit + cudarc crate
- **Architecture:** Double-buffered pipeline — GPU tests batch N while CPU prepares batch N+1

### Priority 3: GPU Sieve Marking via wgpu (v6)
- **Impact:** 10-50x on sieve marking phase
- **Approach:** CPU computes mod_u offsets, GPU kernel marks bitset with atomic ops
- **Requirements:** Any GPU (cross-platform via WebGPU)

### Priority 4: Deeper Sieve for High Digits (v5.2)
- **Impact:** Fewer survivors at 5000+ digits (where SPRP dominates)
- **Approach:** Streaming sieve to 10^9 (50M additional primes) for targets > 3000 digits
- **Trade-off:** 26s extra sieve time per batch but saves thousands of seconds in testing

### Priority 5: Primorial 30030 Wheel (v6+)
- **Impact:** ~90% candidate reduction
- **Approach:** Full SSoZ-style wheel with 114 twin-eligible positions per cycle of 2310
- **Difficulty:** High — requires restructuring sieve inner loops

### Stretch Goal: 12x Faster Than PrimeGrid
- PrimeGrid's twin prime record: 388,342 digits (2016) using LLR + NewPGen
- Their approach: specialized k*2^n+-1 form (faster modular arithmetic)
- Our approach: general 10^d-range twin primes (harder but more general)
- Path: GPU SPRP (CGBN) + GPU sieve (wgpu) + deeper sieve + wheel optimization
