#!/usr/bin/env bash
# Stage 2 -- bulk backbone generation, 900 backbones on the alpha1/alpha2 crop.
#
# Authorised by Harjas 2026-11-07 ("prioritise achieving the objectives of the
# comp", n_backbones "whatever is good for efficiency and performance").
#
# 900, not the configured 3500: measured 64.9 s/backbone, so 3500 is ~63 GPU-h
# against a ~16 GPU-h block. 900 is ~16.2 GPU-h. Volume is not the binding
# constraint anyway -- 32/32 pilot backbones contacted the peptide and 47% hit
# all four hotspots.
#
# Generated in CHUNKS so the run is resumable: RFdiffusion writes each design as
# it finishes, so the chunk loop just counts what is already on disk and asks
# for the remainder. A crash or a VM restart costs at most one chunk.
#
# Length band 65-85 aa keeps every design inside the published 40-100 aa
# molecule category 1. A smaller out-of-band block is generated separately at
# the end for category hedging, per docs/PIPELINE.md Stage 2.
set -uo pipefail
cd "$(dirname "$0")/.."

TARGET="${1:-900}"
CHUNK="${2:-60}"
OUT=gen/rfdiff
mkdir -p "$OUT" logs

log() { echo -e "\n=== [$(date -u +%H:%M)] $* ===" ; }

have() { ls "$OUT"/g_*.pdb 2>/dev/null | wc -l; }

log "Stage 2 generation: target ${TARGET} backbones, chunk ${CHUNK}"
echo "already on disk: $(have)"

while [ "$(have)" -lt "$TARGET" ]; do
  avail=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  if [ "$avail" -lt 5 ]; then
    echo "ABORT: only ${avail}G free, stopping generation"; exit 1
  fi

  n=$(have)
  want=$(( TARGET - n ))
  [ "$want" -gt "$CHUNK" ] && want=$CHUNK

  echo "[$(date -u +%H:%M)] have ${n}/${TARGET}, generating ${want} more (${avail}G free)"
  s=$(date +%s)
  ./rfdiffusion.sh \
    inference.output_prefix="/work/${OUT}/g_${n}" \
    inference.input_pdb=/work/structures/crop/g12d_crop.pdb \
    'contigmap.contigs=[A1-180/0 C1-10/0 65-85]' \
    'ppi.hotspot_res=[C4,C5,C7,C8]' \
    inference.num_designs="$want" \
    hydra.run.dir=/work/.hydra >> logs/gen_stage2.log 2>&1
  rc=$?
  e=$(date +%s)
  echo "  chunk rc=${rc} wall=$((e-s))s  total now $(have)"

  if [ "$rc" -ne 0 ]; then
    echo "  chunk failed; retrying once after a pause"
    sleep 30
    continue
  fi
  # guard against a chunk that produces nothing, to avoid an infinite loop
  if [ "$(have)" -le "$n" ]; then
    echo "ABORT: chunk produced no new backbones"; exit 1
  fi
done

log "Stage 2 generation COMPLETE: $(have) backbones"
du -sh "$OUT"
