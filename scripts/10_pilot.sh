#!/usr/bin/env bash
# Stage 1 pilot, end to end: decoy generation -> MPNN -> select -> three-state
# Boltz-2 -> interface metrics -> manifest.
#
# Gate 1 asks for this to run with no manual intervention, so it is one script.
# Peptide-targeted backbones are generated separately (they take ~35 min and
# were already running when this was written); this picks up from there and
# generates the framework-binder DECOY arm itself.
#
# The decoy arm is not optional. Gate 1 criterion 3 asks whether the
# g12d-minus-empty spread is wider for high-peptide-contact designs than for
# framework-dominated ones. That is a comparison; without the decoys there is
# nothing to compare against and the criterion cannot be evaluated.
#
# Usage: bash scripts/10_pilot.sh [N_DECOY]
set -euo pipefail
cd "$(dirname "$0")/.."

N_DECOY="${1:-8}"
CONDA="/home/HarjasG/miniforge3/etc/profile.d/conda.sh"
# shellcheck disable=SC1090
source "$CONDA"

log() { echo -e "\n=== $* ===" ; }
t0=$(date +%s)

# ---------------------------------------------------------------- 1. decoys
if [ "$(ls pilot/rfd_decoy/d_*.pdb 2>/dev/null | wc -l)" -lt "$N_DECOY" ]; then
  log "1/7 decoy generation (beta2m framework control, n=$N_DECOY)"
  mkdir -p pilot/rfd_decoy
  s=$(date +%s)
  ./rfdiffusion.sh \
    inference.output_prefix=/work/pilot/rfd_decoy/d \
    inference.input_pdb=/work/structures/crop/decoy_crop.pdb \
    'contigmap.contigs=[B1-99/0 65-85]' \
    'ppi.hotspot_res=[B44,B45,B47,B48]' \
    inference.num_designs="$N_DECOY" \
    hydra.run.dir=/work/.hydra > logs/rfd_decoy.log 2>&1
  echo "decoys: $(ls pilot/rfd_decoy/d_*.pdb | wc -l) in $(( $(date +%s)-s ))s"
else
  log "1/7 decoys already present, skipping"
fi

# ------------------------------------------------------------------ 2. MPNN
# CLAUDE.md section 1: every long run is checkpointed and resumable. Each stage
# below skips if its output already exists, so a crash in a later stage does not
# cost the GPU time already spent on an earlier one.
conda activate mpnn
if [ "$(ls pilot/mpnn/basic/seqs/*.fa 2>/dev/null | wc -l)" -gt 0 ]; then
  log "2/7 MPNN peptide-targeted arm already done, skipping"
else
  log "2/7 ProteinMPNN -- peptide-targeted arm"
  python -I scripts/04_mpnn.py --rfdiff pilot/rfd --out pilot/mpnn --num-seq 4
fi

if [ "$(ls pilot/mpnn_decoy/neutral/seqs/*.fa 2>/dev/null | wc -l)" -gt 0 ]; then
  log "2b/7 MPNN decoy arm already done, skipping"
else
  log "2b/7 ProteinMPNN -- decoy arm (neutral pool only; the R/K bias is aimed"
  echo "     at peptide p6, which a beta2m binder never sees)"
  python -I scripts/04_mpnn.py --rfdiff pilot/rfd_decoy --out pilot/mpnn_decoy \
    --num-seq 4 --pools neutral
fi

# ---------------------------------------------------------------- 3. select
conda activate boltz
if [ "$(ls pilot/mpnn/selected/basic/seqs/*.fa 2>/dev/null | wc -l)" -gt 0 ]; then
  log "3/7 selection already done, skipping"
else
  log "3/7 selection"
  python -I scripts/04b_select.py --mpnn pilot/mpnn --pools neutral,basic \
    --top-per-pool 16 --per-backbone 1
fi
if [ "$(ls pilot/mpnn_decoy/selected/neutral/seqs/*.fa 2>/dev/null | wc -l)" -eq 0 ]; then
  python -I scripts/04b_select.py --mpnn pilot/mpnn_decoy --pools neutral \
    --top-per-pool "$N_DECOY" --per-backbone 1
fi

# ------------------------------------------------------- 4. three-state fold
log "4/7 three-state Boltz-2 -- peptide-targeted arm"
python -I scripts/05_boltz3.py \
  --fasta-dir pilot/mpnn/selected/neutral/seqs pilot/mpnn/selected/basic/seqs \
  --out pilot/boltz --tag main

log "5/7 three-state Boltz-2 -- decoy arm"
python -I scripts/05_boltz3.py \
  --fasta-dir pilot/mpnn_decoy/selected/neutral/seqs \
  --out pilot/boltz --tag decoy

# --------------------------------------------------------------- 5. metrics
log "6/7 interface metrics on folded g12d complexes"
# Side chains exist now, so the area-based metrics are valid here -- unlike on
# the raw RFdiffusion backbones, where they are emitted as not_run.
python -I scripts/03_interface_metrics.py \
  --glob 'pilot/boltz/main/g12d/pred/**/*_model_0.pdb' \
  --binder-chain X --out results/processed/interface_main.csv
python -I scripts/03_interface_metrics.py \
  --glob 'pilot/boltz/decoy/g12d/pred/**/*_model_0.pdb' \
  --binder-chain X --out results/processed/interface_decoy.csv

# -------------------------------------------------------------- 6. manifest
log "7/7 manifest"
python -I scripts/06_manifest.py \
  --boltz pilot/boltz/main_three_state.csv \
  --interface results/processed/interface_main.csv \
  --out results/processed/manifest --track pilot_main
python -I scripts/06_manifest.py \
  --boltz pilot/boltz/decoy_three_state.csv \
  --interface results/processed/interface_decoy.csv \
  --out results/processed/manifest_decoy --track pilot_decoy_control

echo -e "\n=== PILOT COMPLETE in $(( ($(date +%s)-t0)/60 )) min ==="
