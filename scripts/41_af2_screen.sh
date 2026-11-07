#!/usr/bin/env bash
# AF2-initial-guess screen -- the PRIMARY design filter.
#
# WHY THIS IS PRIMARY AND BOLTZ IS NOT, measured on 32 pilot designs:
#
#   pae_interaction < 10        2/32
#   plddt_binder > 80          24/32
#   binder_aligned_rmsd < 1.5  12/32
#   target_aligned_rmsd < 5     2/32
#   ALL THREE                   1/32   (basic_p_3_s2, pae 5.21, tgt_rmsd 1.38)
#
# Median target_aligned_rmsd was 40.7 A: most designs FOLD correctly but
# RELOCATE the binder wholesale, and nothing in the Boltz pipeline detected
# that. Worse, the one design clearing AF2 scored only 1/2 on our calibrated
# Boltz gates and would have been rejected, while all five designs passing 2/2
# gates scored pae_interaction 25-28.
#
# pae_interaction < 10 is the canonical published criterion (Bennett et al.
# 2023) and the only one here with experimental backing. Boltz keeps its job --
# the three-state g12d/wt/empty differential, which AF2-ig structurally cannot
# provide since it sees one 2-chain complex -- but it no longer selects.
#
# BOTH pae_interaction AND target_aligned_rmsd are required. basic_p_5_s1
# scored pae 8.91 (a pass) with target_aligned_rmsd 39.3 -- confidently
# predicting an interface that is NOT the designed one. Confidence without
# pose agreement is a different complex that happens to score well.
#
# EVERY MPNN sequence is screened, not a pre-selected subset: MPNN score does
# not predict AF2 success, so pre-selecting one per backbone discarded
# diversity for nothing. At ~10 s/design this is affordable.
#
# predict_pdb.py sorts by length internally (one JAX compile per distinct
# length, ~2 min each, then ~10 s/design) and checkpoints every design, so this
# is resumable.
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh

MPNN_DIR="${1:-gen/mpnn}"
BACKBONES="${2:-gen/rfdiff}"
WORK="${3:-/mnt/comp5data/comp5/af2_screen}"
POOLS="${4:-neutral,basic}"

log() { echo -e "\n=== [$(date -u +%H:%M)] $* ==="; }

mkdir -p "$WORK"

log "1/3 threading every MPNN sequence onto its backbone"
conda activate boltz
if [ ! -f "$WORK/in/threaded_manifest.csv" ]; then
  python -I scripts/40_thread_af2.py --all \
    --mpnn "$MPNN_DIR" --backbones "$BACKBONES" --out "$WORK/in" --pools "$POOLS"
else
  echo "already threaded: $(ls "$WORK"/in/*.pdb 2>/dev/null | wc -l) designs"
fi

N=$(ls "$WORK"/in/*.pdb 2>/dev/null | wc -l)
if [ "$N" -eq 0 ]; then echo "ABORT: nothing threaded"; exit 1; fi
log "2/3 AF2-initial-guess on ${N} designs (resumable)"
conda activate af2ig
python scripts/af2ig/predict_pdb.py \
  -pdbdir "$WORK/in" \
  -outpdbdir "$WORK/out" \
  -scorefilename "$WORK/scores.sc" \
  -checkpoint_name "$WORK/check.point" 2>&1 | grep -vE "WARNING:tensorflow|Instructions for updating|^$"

log "3/3 parsing scores and applying the canonical filter"
conda activate boltz
python -I scripts/42_af2_filter.py --scores "$WORK/scores.sc" \
  --threaded "$WORK/in/threaded_manifest.csv" \
  --out results/processed/af2_screen.csv

echo -e "\n=== AF2 screen done $(date -u +%H:%M) ==="
