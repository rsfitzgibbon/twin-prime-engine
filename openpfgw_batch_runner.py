"""
Bridge compact fixed-n survivor batches into OpenPFGW-oriented inputs.

This script consumes the repo's compact export format, reconstructs exact
survivor k-values, writes OpenPFGW-friendly ABCD/text inputs, and can optionally
invoke an external PFGW executable when available.
"""

import argparse
import json
import re
import subprocess
from pathlib import Path


FORMULA_RE = re.compile(r"(\d+)\*2\^(\d+)([+-]1)")


def decode_varints(data):
    i = 0
    out = []
    while i < len(data):
        shift = 0
        value = 0
        while True:
            if i >= len(data):
                raise ValueError("truncated varint payload")
            byte = data[i]
            i += 1
            value |= (byte & 0x7F) << shift
            if byte < 0x80:
                break
            shift += 7
        out.append(value)
    return out


def load_compact_batch(meta_path: Path):
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    payload_name = meta["encoding"]["payload"]
    payload_path = meta_path.with_name(payload_name)
    raw = payload_path.read_bytes()
    if len(raw) < 8:
        raise ValueError(f"payload too short: {payload_path}")
    survivor_count = int.from_bytes(raw[:8], "little")
    deltas = decode_varints(raw[8:])
    if len(deltas) != survivor_count:
        raise ValueError(
            f"payload survivor count mismatch: header={survivor_count}, decoded={len(deltas)}"
        )

    batch_k0 = int(meta["batch_k0"])
    steps = []
    running = 0
    for delta in deltas:
        running += delta
        steps.append(running)
    survivors = [batch_k0 + 6 * step for step in steps]

    return {
        "meta": meta,
        "meta_path": meta_path,
        "payload_path": payload_path,
        "survivors": survivors,
    }


def write_formula_file(path: Path, survivors, n: int, sign: str):
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for k in survivors:
            f.write(f"{k}*2^{n}{sign}1\n")


def write_abcd_file(path: Path, survivors, n: int, sign: str):
    with path.open("w", encoding="utf-8", newline="\n") as f:
        if not survivors:
            return
        f.write(f"ABCD $a*2^{n}{sign}1 [{survivors[0]}]\n")
        prev = survivors[0]
        for k in survivors[1:]:
            f.write(f"{k - prev}\n")
            prev = k


def write_twin_abcd_file(path: Path, survivors, n: int):
    with path.open("w", encoding="utf-8", newline="\n") as f:
        if not survivors:
            return
        f.write(f"ABCD $a*2^{n}-1 & $a*2^{n}+1 [{survivors[0]}]\n")
        prev = survivors[0]
        for k in survivors[1:]:
            f.write(f"{k - prev}\n")
            prev = k


def prepare_openpfgw_files(batch, output_dir: Path):
    meta = batch["meta"]
    survivors = batch["survivors"]
    n = int(meta["n"])
    stem = batch["meta_path"].stem.replace(".meta", "")
    batch_dir = output_dir / stem
    batch_dir.mkdir(parents=True, exist_ok=True)

    plus_txt = batch_dir / "plus.txt"
    minus_txt = batch_dir / "minus.txt"
    plus_abcd = batch_dir / "plus.abcd"
    minus_abcd = batch_dir / "minus.abcd"
    twin_abcd = batch_dir / "twin.abcd"
    manifest = batch_dir / "manifest.json"

    write_formula_file(plus_txt, survivors, n, "+")
    write_formula_file(minus_txt, survivors, n, "-")
    write_abcd_file(plus_abcd, survivors, n, "+")
    write_abcd_file(minus_abcd, survivors, n, "-")
    write_twin_abcd_file(twin_abcd, survivors, n)

    manifest_payload = {
        "format": "openpfgw_batch_v1",
        "source_format": meta["format"],
        "family": "k*2^n±1",
        "n": n,
        "batch_index": meta["batch_index"],
        "batch_k0": meta["batch_k0"],
        "survivor_count": len(survivors),
        "files": {
            "plus_txt": plus_txt.name,
            "minus_txt": minus_txt.name,
            "plus_abcd": plus_abcd.name,
            "minus_abcd": minus_abcd.name,
            "twin_abcd": twin_abcd.name,
        },
    }
    manifest.write_text(json.dumps(manifest_payload, indent=2), encoding="utf-8")

    return {
        "batch_dir": batch_dir,
        "plus_txt": plus_txt,
        "minus_txt": minus_txt,
        "plus_abcd": plus_abcd,
        "minus_abcd": minus_abcd,
        "twin_abcd": twin_abcd,
        "manifest": manifest,
        "survivors": survivors,
        "n": n,
    }


def render_command(template: str, exe: Path, input_path: Path, output_path: Path):
    return template.format(
        exe=str(exe),
        input=str(input_path),
        output=str(output_path),
    )


def parse_pfgw_hits(text: str, expected_n: int, sign: str):
    hits = {}
    for line in text.splitlines():
        lower = line.lower()
        # PFGW emits many progress lines of the form "PRP: ...".
        # Only lines with an explicit final result ("is ...") count as hits.
        if " is " not in lower:
            continue
        if "prp" not in lower and "prime" not in lower:
            continue
        match = FORMULA_RE.search(line)
        if not match:
            continue
        k, n_str, line_sign = match.groups()
        if int(n_str) != expected_n or line_sign != f"{sign}1":
            continue
        if "composite" in lower:
            continue
        hits[int(k)] = line.strip()
    return hits


def run_pfgw(template: str, exe: Path, input_path: Path, output_path: Path):
    cmd = render_command(template, exe, input_path, output_path)
    proc = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    combined = (proc.stdout or "") + ("\n" + proc.stderr if proc.stderr else "")
    output_path.write_text(combined, encoding="utf-8")
    return {
        "returncode": proc.returncode,
        "command": cmd,
        "output_path": str(output_path),
        "text": combined,
    }


def process_batch(batch_meta: Path, output_dir: Path, pfgw_exe: Path | None,
                  plus_template: str, minus_template: str, prepare_only: bool):
    batch = load_compact_batch(batch_meta)
    prepared = prepare_openpfgw_files(batch, output_dir)

    result = {
        "batch_meta": batch_meta.name,
        "batch_dir": prepared["batch_dir"].name,
        "n": prepared["n"],
        "survivor_count": len(prepared["survivors"]),
        "artifacts": {
            "manifest": prepared["manifest"].name,
            "plus_txt": prepared["plus_txt"].name,
            "minus_txt": prepared["minus_txt"].name,
            "plus_abcd": prepared["plus_abcd"].name,
            "minus_abcd": prepared["minus_abcd"].name,
            "twin_abcd": prepared["twin_abcd"].name,
        },
    }

    if prepare_only or pfgw_exe is None:
        result["mode"] = "prepare_only"
        return result

    twin_log = prepared["batch_dir"] / "twin.pfgw.log"
    twin_run = run_pfgw(plus_template, pfgw_exe, prepared["twin_abcd"], twin_log)

    plus_hits = parse_pfgw_hits(twin_run["text"], prepared["n"], "+")
    minus_hits = parse_pfgw_hits(twin_run["text"], prepared["n"], "-")
    twin_hits = []
    for k in sorted(set(plus_hits) & set(minus_hits)):
        twin_hits.append(
            {
                "k": str(k),
                "formula_minus": f"{k}*2^{prepared['n']}-1",
                "formula_plus": f"{k}*2^{prepared['n']}+1",
                "minus_line": minus_hits[k],
                "plus_line": plus_hits[k],
            }
        )

    result.update(
        {
            "mode": "executed",
            "twin_run": {
                "returncode": twin_run["returncode"],
                "command": twin_run["command"],
                "output_path": twin_log.name,
                "hit_count": len(plus_hits),
                "plus_hit_count": len(plus_hits),
                "minus_hit_count": len(minus_hits),
            },
            "twin_hit_count": len(twin_hits),
            "twin_hits": twin_hits,
        }
    )
    return result


def discover_meta_files(compact_dir: Path):
    return sorted(compact_dir.glob("*.meta.json"))


def main():
    parser = argparse.ArgumentParser(description="Prepare or run OpenPFGW jobs from compact fixed-n batches.")
    parser.add_argument("--compact-dir", type=Path, required=True,
                        help="Directory containing compact_export batch files (*.meta.json, *.kdelta.bin).")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory for reconstructed OpenPFGW job files.")
    parser.add_argument("--batch-limit", type=int, default=None,
                        help="Optional limit on number of batches to process.")
    parser.add_argument("--prepare-only", action="store_true",
                        help="Only reconstruct OpenPFGW inputs; do not execute PFGW.")
    parser.add_argument("--pfgw-exe", type=Path, default=None,
                        help="Path to pfgw executable. If omitted, runs in prepare-only mode.")
    parser.add_argument("--plus-template", type=str, default="\"{exe}\" \"{input}\"",
                        help="Command template for plus-side runs. Placeholders: {exe} {input} {output}")
    parser.add_argument("--minus-template", type=str, default="\"{exe}\" \"{input}\"",
                        help="Command template for minus-side runs. Placeholders: {exe} {input} {output}")
    parser.add_argument("--json-out", type=Path, default=None,
                        help="Optional JSON summary output.")
    args = parser.parse_args()

    meta_files = discover_meta_files(args.compact_dir)
    if args.batch_limit is not None:
        meta_files = meta_files[: args.batch_limit]

    if not meta_files:
        raise SystemExit(f"No compact batch metadata found in {args.compact_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = [
        process_batch(
            batch_meta=meta_path,
            output_dir=args.output_dir,
            pfgw_exe=args.pfgw_exe,
            plus_template=args.plus_template,
            minus_template=args.minus_template,
            prepare_only=args.prepare_only or args.pfgw_exe is None,
        )
        for meta_path in meta_files
    ]

    summary = {
        "compact_dir": str(args.compact_dir),
        "output_dir": str(args.output_dir),
        "batch_count": len(results),
        "total_survivors": sum(item["survivor_count"] for item in results),
        "prepare_only": args.prepare_only or args.pfgw_exe is None,
        "results": results,
    }

    if args.json_out:
        args.json_out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"Saved JSON to {args.json_out}")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
