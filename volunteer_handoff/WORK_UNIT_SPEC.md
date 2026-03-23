# Work Unit Spec

## Work Unit

Each work unit is one deterministic contiguous band of corridor indices `m`, where candidate pairs are:

- `p = 6m - 1`
- `q = 6m + 1`

Required fields:

- `version`
- `job_id`
- `target_digits`
- `m_start`
  - decimal string, so very large corridor indices are representable exactly
- `batch_size`
- `sieve_depth`
- `sieve_limit`
- `created_utc`

Optional fields:

- `notes`
- `priority`

## Result

Required result fields:

- `job_id`
- `status`
- `elapsed_seconds`
- `total_raw`
- `survivors`
- `tested`
- `found_count`
- `found_pairs`

Each found pair entry should include:

- `m`
- `p`
- `q`
- `digits`

## Determinism rules

- same work unit must produce the same survivor set
- same work unit must produce the same candidate pairs
- output order should be sorted by `m`
- no random sampling inside assigned work

## Validation model

- duplicate each work unit on at least two hosts
- compare:
  - survivor count
  - tested count
  - found pair list
- if mismatch:
  - resend to a third host
  - quarantine the inconsistent result
