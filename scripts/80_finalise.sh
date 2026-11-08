#!/usr/bin/env bash
# Endgame: Boltz three-state on survivors -> novelty -> rank -> submission CSV.
#
# NEVER SUBMITS. configs/challenge.yaml:submission.auto_submit is false and
# CLAUDE.md section 2.5 reserves submission for Harjas. This produces
# submission/comp5_track3.csv plus submission/AUDIT.md for his approval.
#
# Each stage is skipped if its output exists, so this is safe to re-run.
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh
log() { echo -e "\n=== [$(date -u +%H:%M)] $* ==="; }

TOPN="${1:-90}"

if [ ! -f /mnt/comp5data/comp5/boltz_final/final_three_state.csv ]; then
  log "1/4 three-state Boltz-2 on AF2 survivors (5 averaged samples)"
  bash scripts/50_boltz_survivors.sh "$TOPN"
else
  log "1/4 Boltz already done, skipping"
fi

conda activate boltz
log "2/4 de novo novelty screen"
python -I scripts/61_novelty.py

log "3/4 hard gates, diversity cap, lexicographic ranking"
python -I scripts/60_rank.py

log "4/4 submission CSV + audit (does NOT submit)"
python -I scripts/70_submission.py
rc=$?

log "FINALISE COMPLETE (submission audit exit ${rc}; non-zero means blocking problems)"
echo "Harjas must review submission/AUDIT.md before submitting."
