"""
Shard compact fixed-n export batches into deterministic worker directories.

This reassigns each *.meta.json plus its sibling *.kdelta.bin payload to
`parts` destination directories based on batch index modulo `parts`.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Shard compact export batches into worker directories.")
    parser.add_argument("--source-root", action="append", required=True,
                        help="Root directory to scan recursively for *.meta.json batches. Repeatable.")
    parser.add_argument("--dest-root", required=True, help="Destination root for worker part directories.")
    parser.add_argument("--parts", type=int, required=True, help="Number of destination parts.")
    parser.add_argument("--part-prefix", default="part", help="Directory prefix for parts.")
    parser.add_argument("--move", action="store_true", help="Move files instead of copying them.")
    parser.add_argument("--clean-dest", action="store_true", help="Delete destination root before sharding.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dest_root = Path(args.dest_root)
    if args.clean_dest and dest_root.exists():
        shutil.rmtree(dest_root)
    dest_root.mkdir(parents=True, exist_ok=True)

    meta_paths: list[Path] = []
    for root in args.source_root:
        meta_paths.extend(sorted(Path(root).rglob("*.meta.json")))

    assignments = []
    for meta in meta_paths:
        try:
            payload = json.loads(meta.read_text(encoding="utf-8"))
        except Exception:
            continue
        batch_index = int(payload["batch_index"])
        part = batch_index % args.parts
        part_dir = dest_root / f"{args.part_prefix}{part:03d}"
        part_dir.mkdir(parents=True, exist_ok=True)
        payload_name = payload["encoding"]["payload"]
        payload_path = meta.with_name(payload_name)

        meta_dest = part_dir / meta.name
        payload_dest = part_dir / payload_path.name
        op = shutil.move if args.move else shutil.copy2
        op(str(meta), str(meta_dest))
        if payload_path.exists():
            op(str(payload_path), str(payload_dest))

        assignments.append((batch_index, str(meta_dest)))

    summary = {
        "parts": args.parts,
        "dest_root": str(dest_root),
        "count": len(assignments),
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
