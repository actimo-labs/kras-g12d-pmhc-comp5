#!/usr/bin/env bash
# Three-state Boltz-2 on the AF2 survivors -- the only stage that scores the
# competition's actual objectives.
#
# DIVISION OF LABOUR. AF2-ig answered "does this sequence fold into this
# backbone and hold this interface", on a 2-chain crop. It cannot answer the
# ranking questions, because it never sees the wild-type peptide or the empty
# groove. Boltz does: g12d / wt / empty, full 4-chain complex, which is where
# iptm_margin_wt (objective 1) and iptm_margin_empty (objective 2) come from.
#
# FIVE DIFFUSION SAMPLES, NOT ONE. Measured: at one sample iptm_margin_wt has a
# between-seed Spearman rho of +0.239 -- the ranking is near-random under a
# reseed. At five averaged samples rho is +0.846 and margin_empty +0.881.
# Single-sample margins in this repo are noise and must never be ranked on.
# Cost is 2.41x, not 5x, because Boltz shares the trunk and repeats only the
# diffusion module.
#
# Survivors are capped by AF2 quality before folding. With ~135+ unique hits and
# 20 submission slots there is no need to fold all of them, and the cap keeps
# this stage inside its block.
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh
conda activate boltz

TOPN="${1:-90}"
WORK="${2:-/mnt/comp5data/comp5/boltz_final}"
log() { echo -e "\n=== [$(date -u +%H:%M)] $* ==="; }

log "1/3 collecting AF2 survivors across all arms, deduplicated"
python -I scripts/51_collect_survivors.py --top "$TOPN" --out "$WORK/seqs"
N=$(ls "$WORK"/seqs/*.fa 2>/dev/null | wc -l)
if [ "$N" -eq 0 ]; then echo "ABORT: no survivors"; exit 1; fi
echo "to fold: $N designs x 3 states x 5 samples"

log "2/3 three-state Boltz-2, 5 averaged diffusion samples"
python -I scripts/05_boltz3.py --fasta-dir "$WORK/seqs" \
  --out "$WORK" --tag final \
  --diffusion-samples 5 --average-samples

log "3/3 interface metrics on the folded g12d complexes"
python -I scripts/03_interface_metrics.py \
  --glob "$WORK/final/g12d/pred/**/*_model_0.pdb" \
  --binder-chain X --out results/processed/interface_final.csv

python -I scripts/06_manifest.py \
  --boltz "$WORK/final_three_state.csv" \
  --interface results/processed/interface_final.csv \
  --out results/processed/manifest_final --track track3

echo -e "\n=== Boltz on survivors done $(date -u +%H:%M) ==="
