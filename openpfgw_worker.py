"""
Continuous worker that consumes compact fixed-n batches and feeds them into
OpenPFGW one batch at a time with a persistent checkpoint.
"""

import argparse
import json
import time
from pathlib import Path

from openpfgw_batch_runner import discover_meta_files, process_batch


def load_state(path: Path):
    if not path.exists():
        return {
            "processed": [],
            "runs": 0,
            "total_survivors": 0,
            "total_twin_hits": 0,
            "last_batch": None,
        }
    return json.loads(path.read_text(encoding="utf-8"))


def save_state(path: Path, state):
    path.write_text(json.dumps(state, indent=2), encoding="utf-8")


def append_jsonl(path: Path, item):
    with path.open("a", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(item))
        f.write("\n")


def pending_meta_files(compact_dir: Path, processed: set[str]):
    return [meta for meta in discover_meta_files(compact_dir) if meta.name not in processed]


def process_pending(args, state):
    processed = set(state["processed"])
    pending = pending_meta_files(args.compact_dir, processed)
    if args.max_batches_per_pass is not None:
        pending = pending[: args.max_batches_per_pass]

    results = []
    for meta_path in pending:
        result = process_batch(
            batch_meta=meta_path,
            output_dir=args.output_dir,
            pfgw_exe=args.pfgw_exe,
            plus_template=args.plus_template,
            minus_template=args.minus_template,
            prepare_only=args.prepare_only,
        )
        results.append(result)
        state["processed"].append(meta_path.name)
        state["runs"] += 1
        state["total_survivors"] += result["survivor_count"]
        state["total_twin_hits"] += result.get("twin_hit_count", 0)
        state["last_batch"] = meta_path.name
        save_state(args.state_file, state)
        if args.result_log is not None:
            append_jsonl(args.result_log, result)
    return results


def main():
    parser = argparse.ArgumentParser(description="Continuous OpenPFGW worker over compact fixed-n survivor batches.")
    parser.add_argument("--compact-dir", type=Path, required=True,
                        help="Directory containing compact_export batches.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory for reconstructed PFGW jobs and logs.")
    parser.add_argument("--pfgw-exe", type=Path, required=True,
                        help="Path to pfgw executable.")
    parser.add_argument("--state-file", type=Path, required=True,
                        help="Checkpoint/state JSON file.")
    parser.add_argument("--result-log", type=Path, default=None,
                        help="Optional JSONL append log of processed batch results.")
    parser.add_argument("--poll-seconds", type=float, default=30.0,
                        help="Sleep interval when no new batches are found.")
    parser.add_argument("--once", action="store_true",
                        help="Process current pending batches once and exit.")
    parser.add_argument("--prepare-only", action="store_true",
                        help="Prepare jobs but do not execute PFGW.")
    parser.add_argument("--max-batches-per-pass", type=int, default=None,
                        help="Optional cap on batches processed per scan loop.")
    parser.add_argument("--plus-template", type=str, default="\"{exe}\" \"{input}\"",
                        help="Command template for OpenPFGW invocation.")
    parser.add_argument("--minus-template", type=str, default="\"{exe}\" \"{input}\"",
                        help="Kept for interface compatibility; twin jobs use one invocation.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.state_file.parent.mkdir(parents=True, exist_ok=True)
    if args.result_log is not None:
        args.result_log.parent.mkdir(parents=True, exist_ok=True)

    state = load_state(args.state_file)

    while True:
        results = process_pending(args, state)
        if results:
            print(json.dumps({
                "processed_now": len(results),
                "total_runs": state["runs"],
                "total_survivors": state["total_survivors"],
                "total_twin_hits": state["total_twin_hits"],
                "last_batch": state["last_batch"],
            }, indent=2))
        elif args.once:
            print("No pending compact batches.")
            break
        else:
            time.sleep(max(0.1, args.poll_seconds))
            continue

        if args.once:
            break


if __name__ == "__main__":
    main()
