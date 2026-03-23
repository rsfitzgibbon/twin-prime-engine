# OpenPFGW Run Results

## Status

The repo is now connected to a real local OpenPFGW install:

- executable: `C:\Users\rsfit\Downloads\pfgw_win_4.1.7\distribution\pfgw64.exe`

The bridge script
[openpfgw_batch_runner.py](C:/Users/rsfit/twin-prime-engine/openpfgw_batch_runner.py)
was validated in three stages.

## 1. Trivial twin sanity check

Using a compact batch for:

- `n = 2`
- `k = 3`

OpenPFGW correctly detected:

- `3*2^2-1 = 11`
- `3*2^2+1 = 13`

Recorded in:

- [pfgw_tiny_jobs_summary.json](C:/Users/rsfit/twin-prime-engine/pfgw_tiny_jobs_summary.json)

## 2. Real nontrivial batch at n = 1000

Input campaign:

- `n = 1000`
- `k_batch_size = 1000`
- base sieve `50,000`
- exact post-sieve `200,000`

Front-end summary:

- base survivors: `22`
- post-sieved survivors: `17`

OpenPFGW completed the batch and found two twin PRP hits:

- `915*2^1000-1` and `915*2^1000+1`
- `1197*2^1000-1` and `1197*2^1000+1`

Recorded in:

- [pfgw_n1000_compact_summary.json](C:/Users/rsfit/twin-prime-engine/pfgw_n1000_compact_summary.json)
- [pfgw_n1000_jobs_summary.json](C:/Users/rsfit/twin-prime-engine/pfgw_n1000_jobs_summary.json)

## 3. 388k-digit probe

Input campaign:

- `n = 1,288,907`
- `k_batch_size = 64`
- base sieve `50,000`
- exact post-sieve `200,000`

Front-end summary:

- base survivors: `2`
- post-sieved survivors: `2`

The front-end export completed instantly.

The OpenPFGW run on this machine did **not** complete within `300s`, even for
that 2-survivor probe batch.

That means:

- the backend integration is working
- but `388k`-digit PRP work is still too heavy for this machine as a local worker

Recorded in:

- [pfgw_388k_tiny_compact_summary.json](C:/Users/rsfit/twin-prime-engine/pfgw_388k_tiny_compact_summary.json)

## Conclusion

The architecture is now correct:

1. exact Rust front-end
2. exact post-sieve
3. compact export
4. OpenPFGW backend bridge

The remaining blocker is no longer software wiring.
It is raw backend compute for very large PRP testing.

So the next real step toward the next largest twin prime is:

- run this backend path on stronger or distributed workers
- keep the exact post-sieve enabled to reduce backend load before export
