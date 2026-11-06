#!/usr/bin/env bash
# Overnight run, 2026-11-07. DIAGNOSIS AND CALIBRATION ONLY -- NOT Stage 2 bulk.
#
# Gate 1 FAILED two of four criteria (results/RUNLOG.md Stage 06). docs/PIPELINE.md
# Stage 1 is explicit: "If Gate 1 fails: stop, report, do not proceed to bulk."
# The brief also requires Harjas's go-ahead before bulk generation, and he is
# asleep. So this script runs only the work that diagnoses the two failures and
# calibrates the thresholds -- all of which Stage 4 requires anyway, so none of
# it is wasted whichever way the Stage 2 decision goes.
#
# What it does NOT do: no 3,500-backbone generation, no new model or framework,
# no writes outside this repo, no network calls, no git push.
#
# Budget: ~2 GPU-h of the ~71 h total. Every stage skips if its output exists.
set -uo pipefail          # NOT -e: one failed stage must not kill the night
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh

N_EXTRA_DECOY="${1:-16}"
log() { echo -e "\n=== [$(date -u +%H:%M)] $* ===" ; }
t0=$(date +%s)

# Guard: refuse to run if the disk is nearly gone. The box is at 100% because of
# another user's 1.6 T; we have ~15 G and this run needs ~2 G.
avail=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
if [ "$avail" -lt 4 ]; then
  echo "ABORT: only ${avail}G free, refusing to start"; exit 1
fi
echo "disk free at start: ${avail}G"

# ---------------------------------------------------------- 1. noise floor
# The single most important missing number. Pilot margins sit near zero
# (median +0.028) with sd 0.157; without a run-to-run variance estimate none of
# that spread can be called signal.
if [ ! -f pilot/boltz_rep/rep32_three_state.csv ]; then
  log "1/6 seed replicate, all 32 main designs, 3 states, seed 1234"
  conda activate boltz
  python -I scripts/05_boltz3.py \
    --fasta-dir pilot/mpnn/selected/neutral/seqs pilot/mpnn/selected/basic/seqs \
    --out pilot/boltz_rep --tag rep32 --seed 1234
else
  log "1/6 seed replicate already done, skipping"
fi

# ------------------------------------------------- 2. sequence negatives
if [ ! -d pilot/controls/shuffled/seqs ]; then
  log "2/6 building composition-matched sequence controls"
  conda activate boltz
  python -I scripts/07_controls.py --from results/processed/manifest.csv
else
  log "2/6 controls already built, skipping"
fi

for kind in shuffled reversed; do
  if [ ! -f "pilot/boltz/${kind}_three_state.csv" ]; then
    log "3/6 folding ${kind} controls, 3 states"
    conda activate boltz
    python -I scripts/05_boltz3.py --fasta-dir "pilot/controls/${kind}/seqs" \
      --out pilot/boltz --tag "$kind" --no-composition-filter
    python -I scripts/03_interface_metrics.py \
      --glob "pilot/boltz/${kind}/g12d/pred/**/*_model_0.pdb" \
      --binder-chain X --out "results/processed/interface_${kind}.csv"
    python -I scripts/06_manifest.py \
      --boltz "pilot/boltz/${kind}_three_state.csv" \
      --interface "results/processed/interface_${kind}.csv" \
      --out "results/processed/manifest_${kind}" --track "control_${kind}"
  else
    log "3/6 ${kind} controls already folded, skipping"
  fi
done

# ------------------------------------------- 4. more framework decoys
# Only 5 of 8 pilot decoys had zero peptide burial after folding, so the
# framework arm is calibrating on n=5. More decoys give that arm real power.
if [ "$(ls pilot/rfd_decoy2/e_*.pdb 2>/dev/null | wc -l)" -lt "$N_EXTRA_DECOY" ]; then
  log "4/6 generating ${N_EXTRA_DECOY} more beta2m framework decoys"
  mkdir -p pilot/rfd_decoy2
  ./rfdiffusion.sh \
    inference.output_prefix=/work/pilot/rfd_decoy2/e \
    inference.input_pdb=/work/structures/crop/decoy_crop.pdb \
    'contigmap.contigs=[B1-99/0 65-85]' \
    'ppi.hotspot_res=[B44,B45,B47,B48]' \
    inference.num_designs="$N_EXTRA_DECOY" \
    hydra.run.dir=/work/.hydra >> logs/rfd_decoy2.log 2>&1
fi

if [ ! -f pilot/boltz/decoy2_three_state.csv ]; then
  log "5/6 MPNN + fold the extra decoys"
  conda activate mpnn
  [ -d pilot/mpnn_decoy2/neutral/seqs ] || \
    python -I scripts/04_mpnn.py --rfdiff pilot/rfd_decoy2 --out pilot/mpnn_decoy2 \
      --num-seq 4 --pools neutral
  conda activate boltz
  [ -d pilot/mpnn_decoy2/selected/neutral/seqs ] || \
    python -I scripts/04b_select.py --mpnn pilot/mpnn_decoy2 --pools neutral \
      --top-per-pool "$N_EXTRA_DECOY" --per-backbone 1
  python -I scripts/05_boltz3.py --fasta-dir pilot/mpnn_decoy2/selected/neutral/seqs \
    --out pilot/boltz --tag decoy2
  python -I scripts/03_interface_metrics.py \
    --glob 'pilot/boltz/decoy2/g12d/pred/**/*_model_0.pdb' \
    --binder-chain X --out results/processed/interface_decoy2.csv
  python -I scripts/06_manifest.py --boltz pilot/boltz/decoy2_three_state.csv \
    --interface results/processed/interface_decoy2.csv \
    --out results/processed/manifest_decoy2 --track pilot_decoy_control
else
  log "5/6 extra decoys already folded, skipping"
fi

# ------------------------------------------------- 6. calibrate + re-gate
log "6/6 calibration and Gate 1 re-evaluation"
conda activate boltz

# merge both decoy manifests for the framework-negative family
python - <<'PY'
import csv, pathlib
out = pathlib.Path("results/processed/manifest_decoy_all.csv")
rows, hdr = [], None
for f in ("results/processed/manifest_decoy.csv",
          "results/processed/manifest_decoy2.csv"):
    p = pathlib.Path(f)
    if not p.exists():
        continue
    with open(p) as fh:
        r = csv.DictReader(fh)
        hdr = hdr or r.fieldnames
        rows += list(r)
if rows:
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader(); w.writerows(rows)
    print(f"merged {len(rows)} decoy rows -> {out}")
PY

python -I scripts/08_calibrate.py \
  --decoy results/processed/manifest_decoy_all.csv \
  --replicate pilot/boltz_rep/rep32_three_state.csv \
  2>&1 | tee logs/calibration.log

python -I scripts/11_gate1.py 2>&1 | tee logs/gate1_final.log

echo -e "\n=== OVERNIGHT COMPLETE in $(( ($(date +%s)-t0)/60 )) min ==="
echo "disk free at end: $(df -BG --output=avail / | tail -1 | tr -dc '0-9')G"
