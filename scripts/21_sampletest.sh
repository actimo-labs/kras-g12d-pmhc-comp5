#!/usr/bin/env bash
# Does averaging over diffusion samples rescue the ipTM margins?
#
# THE PROBLEM. Measured at one sample per design, n=32 (RUNLOG Stage 08):
#     iptm_margin_wt     S/N 1.16x   Spearman rho between seeds +0.239
#     iptm_margin_empty  S/N 2.37x   Spearman rho between seeds +0.853
# We rank designs and never use the absolute value, so rho is what decides
# usability. margin_wt's ranking is near-random under a reseed; S/N alone
# (1.16x) made it look borderline-survivable and was misleading.
#
# This matters more than any other pilot result because mutation specificity is
# objective 1 and ranks FIRST in the official lexicographic order.
#
# Two ways out; this measures the cheap one:
#   (a) average n diffusion samples in ONE run -- model loads once, ~sqrt(n)
#       noise cut if the samples are independent. Cheap.
#   (b) n independent seeded runs -- n x cost. At n=5 that is ~42 GPU-h for
#       Stage 3 alone. Infeasible.
#
# Runs 12 designs spanning the margin range, 3 states, diffusion_samples=N,
# TWICE with different seeds, so the averaged estimator's noise is measured the
# same way the single-sample noise was.
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh
conda activate boltz

N_SAMP="${1:-5}"
echo "=== [$(date -u +%H:%M)] diffusion-sample averaging test, n=${N_SAMP} ==="

SEL=$(python -I scripts/21c_pick_span.py)
echo "designs: $SEL"

if [ ! -f pilot/boltz_samp/samp5_three_state.csv ]; then
  python -I scripts/05_boltz3.py \
    --fasta-dir pilot/mpnn/selected/neutral/seqs pilot/mpnn/selected/basic/seqs \
    --out pilot/boltz_samp --tag samp5 --only "$SEL" \
    --diffusion-samples "$N_SAMP" --average-samples --seed 99
fi

if [ ! -f pilot/boltz_samp/samp5b_three_state.csv ]; then
  python -I scripts/05_boltz3.py \
    --fasta-dir pilot/mpnn/selected/neutral/seqs pilot/mpnn/selected/basic/seqs \
    --out pilot/boltz_samp --tag samp5b --only "$SEL" \
    --diffusion-samples "$N_SAMP" --average-samples --seed 4321
fi

python -I scripts/21b_sample_analysis.py

echo -e "\n=== sample test done $(date -u +%H:%M) ==="
