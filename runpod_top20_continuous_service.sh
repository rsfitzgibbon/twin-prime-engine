#!/usr/bin/env bash
set -euo pipefail

ROOT="/workspace/twin-prime-engine"
PFGW="/workspace/openpfgw/pfgw64"
PARTS_ROOT="/workspace/top20_parts096"
JOBS_ROOT="/workspace/top20_jobs96"
STAGE_ROOT="/workspace/top20_stage"

pkill -f "${ROOT}/openpfgw_worker.py --compact-dir ${PARTS_ROOT}/" || true
pkill -f "${PFGW} ${JOBS_ROOT}/" || true
pkill -f "${ROOT}/rust-engine/target/x86_64-unknown-linux-gnu/release/fixed_n_campaign --n 240000" || true
pkill -f "shard_compact_batches.py --source-root ${STAGE_ROOT}" || true

rm -rf "${JOBS_ROOT}" "${STAGE_ROOT}"
rm -f /workspace/top20_stream.events.jsonl /workspace/top20_stream.checkpoint.json /workspace/top20_stream.stdout.log /workspace/top20_stream.pid
rm -f /workspace/top20_distributor.stdout.log /workspace/top20_distributor.pid

mkdir -p "${JOBS_ROOT}" "${STAGE_ROOT}"

for d in "${PARTS_ROOT}"/part*; do
  [[ -d "$d" ]] || continue
  name=$(basename "$d")
  count=$(find "$d" -maxdepth 1 -name '*.meta.json' | wc -l)
  [[ "$count" -gt 0 ]] || continue
  nohup python3 "${ROOT}/openpfgw_worker.py" \
    --compact-dir "$d" \
    --output-dir "${JOBS_ROOT}/${name}" \
    --pfgw-exe "${PFGW}" \
    --state-file "${JOBS_ROOT}/${name}.state.json" \
    --result-log "${JOBS_ROOT}/${name}.results.jsonl" \
    --poll-seconds 5 \
    >"${JOBS_ROOT}/${name}.stdout.log" 2>&1 &
  echo $! > "${JOBS_ROOT}/${name}.pid"
done

nohup "${ROOT}/rust-engine/target/x86_64-unknown-linux-gnu/release/fixed_n_campaign" \
  --n 240000 \
  --backend compact_export \
  --k-start 3 \
  --k-batch-size 4096 \
  --sieve-limit 50000 \
  --post-sieve-limit 200000 \
  --start-batch 2260 \
  --status-every 25 \
  --event-log /workspace/top20_stream.events.jsonl \
  --checkpoint-out /workspace/top20_stream.checkpoint.json \
  --export-dir "${STAGE_ROOT}" \
  --max-seconds 0 \
  >/workspace/top20_stream.stdout.log 2>&1 &
echo $! > /workspace/top20_stream.pid

nohup bash -lc "while true; do cd ${ROOT} && python3 shard_compact_batches.py --source-root ${STAGE_ROOT} --dest-root ${PARTS_ROOT} --parts 96 --move; sleep 5; done" \
  >/workspace/top20_distributor.stdout.log 2>&1 &
echo $! > /workspace/top20_distributor.pid

echo "workers=$(find "${JOBS_ROOT}" -maxdepth 1 -name 'part*.pid' | wc -l)"
echo -n "producer_pid="
cat /workspace/top20_stream.pid
echo -n "distributor_pid="
cat /workspace/top20_distributor.pid
