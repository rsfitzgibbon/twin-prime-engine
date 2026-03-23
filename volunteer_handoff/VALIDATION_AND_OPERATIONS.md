# Validation And Operations

## Recommended validation policy

- send each work unit to 2 hosts
- accept only if both result payloads match on:
  - `survivors`
  - `tested`
  - `found_pairs`
- if there is disagreement, send the same unit to a third host

## Recommended work-unit size

Choose units small enough that:

- failures are cheap
- timeouts are bounded
- duplicate validation is affordable

For a first external pilot:

- 100 to 1000 digits:
  - `batch_size = 100,000` to `500,000`
- 1000 to 5000 digits:
  - `batch_size = 200,000` to `2,000,000`
- 5000+ digits:
  - start smaller and benchmark per platform

## Deployment path

1. Private coordinator with trusted hosts
2. Duplicate validation turned on by default
3. Package signed worker binaries
4. External pilot with a small user group
5. Public BOINC-style deployment

## Important integration note

The current reference worker rebuilds sieve tables per run.
That is acceptable for demonstration and small pilots.

The resident Rust worker added to this bundle already reuses:

- small-prime tables
- extended prime tables
- inverse tables

That avoids paying setup cost for every job in the same process.

For production deployment, the next step is persistent checkpoint/restart on top of that resident worker.
