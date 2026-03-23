import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


SIEVE_LIMITS = {
    "base": 1_000_000,
    "extended": 100_000_000,
}


def main():
    parser = argparse.ArgumentParser(description="Create a deterministic twin-prime work unit.")
    parser.add_argument("--digits", type=int, required=True, help="Target digit count.")
    parser.add_argument("--m-start", type=int, required=True, help="Starting corridor index m.")
    parser.add_argument("--batch-size", type=int, required=True, help="Number of consecutive m values.")
    parser.add_argument("--sieve-depth", choices=["base", "extended"], default="extended")
    parser.add_argument("--job-id", type=str, default=None, help="Optional explicit job id.")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON path.")
    args = parser.parse_args()

    if args.digits < 2:
        raise SystemExit("digits must be >= 2")
    if args.m_start < 1:
        raise SystemExit("m-start must be >= 1")
    if args.batch_size < 1:
        raise SystemExit("batch-size must be >= 1")

    if args.job_id is None:
        job_id = f"tp-{args.digits}d-m{args.m_start}-n{args.batch_size}"
    else:
        job_id = args.job_id

    payload = {
        "version": 1,
        "job_id": job_id,
        "target_digits": args.digits,
        "m_start": str(args.m_start),
        "batch_size": args.batch_size,
        "sieve_depth": args.sieve_depth,
        "sieve_limit": SIEVE_LIMITS[args.sieve_depth],
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"Saved work unit to {args.out}")


if __name__ == "__main__":
    main()
