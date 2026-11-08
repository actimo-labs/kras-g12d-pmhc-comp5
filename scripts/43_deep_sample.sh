#!/usr/bin/env bash
# Deep sequence sampling on backbones that have ALREADY produced an AF2 hit.
#
# WHY. 22 full-criteria hits in the main screen span only 12 distinct backbones
# out of 900. Once a backbone works it tends to work for several sequences,
# which is the documented behaviour of this protocol: the backbone fixes the
# binding geometry and sequence design then succeeds or fails against it. So
# the cheapest remaining yield is more sequences on proven geometry, not more
# geometry.
#
# WHY A HIGHER TEMPERATURE. The main run used sampling_temp 0.1, where
# ProteinMPNN is near-deterministic: at just 2 sequences per pool it already
# returned 6.4% duplicates, and the +1.2 positional R/K bias failed to change
# the sequence 12-14% of the time because it could not flip the argmax. Asking
# for 24 sequences at temp 0.1 would mostly return the same sequence 24 times.
# Deep sampling therefore runs at temp 0.25.
#
# The trade is explicit: higher temperature buys sequence diversity at the cost
# of per-sequence MPNN confidence. That is the right trade here because AF2 is
# the arbiter, not the MPNN score -- and the pilot showed MPNN score does not
# predict AF2 success at all. Designs are screened, not trusted.
#
# Hit backbones are read from the AF2 screen output at run time rather than
# hardcoded, so re-running after the beta and long arms picks up their hits too.
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh

NSEQ="${1:-24}"
TEMP="${2:-0.25}"
WORK="${3:-/mnt/comp5data/comp5/deep}"
log() { echo -e "\n=== [$(date -u +%H:%M)] $* ==="; }

mkdir -p "$WORK"
conda activate boltz

log "1/4 collecting hit backbones from every completed AF2 screen"
python -I scripts/44_hit_backbones.py --out "$WORK/hit_backbones.txt"
N_BB=$(wc -l < "$WORK/hit_backbones.txt")
if [ "$N_BB" -eq 0 ]; then echo "ABORT: no hit backbones yet"; exit 1; fi
echo "hit backbones: $N_BB"

log "2/4 staging those backbones"
mkdir -p "$WORK/bb"
while read -r line; do
  src="${line%%:*}"; bb="${line##*:}"
  [ -f "$src/$bb.pdb" ] && cp "$src/$bb.pdb" "$src/$bb.trb" "$WORK/bb/" 2>/dev/null
done < "$WORK/hit_backbones.txt"
echo "staged: $(ls "$WORK"/bb/*.pdb 2>/dev/null | wc -l)"

log "3/4 MPNN deep sample: ${NSEQ} seqs/pool at temp ${TEMP}"
conda activate mpnn
python -I scripts/04_mpnn.py --rfdiff "$WORK/bb" --out "$WORK/mpnn" \
  --num-seq "$NSEQ" --temp "$TEMP"

log "4/4 thread + AF2 screen"
conda activate boltz
python -I scripts/40_thread_af2.py --all --mpnn "$WORK/mpnn" \
  --backbones "$WORK/bb" --out "$WORK/in" \
  --pools "sol_neutral,sol_basic,van_neutral,van_basic"

conda activate af2ig
python scripts/af2ig/predict_pdb.py \
  -pdbdir "$WORK/in" -outpdbdir "$WORK/out" \
  -scorefilename "$WORK/scores.sc" -checkpoint_name "$WORK/check.point" \
  2>&1 | grep -vE "WARNING:tensorflow|Instructions for updating|^$"

conda activate boltz
python -I scripts/42_af2_filter.py --scores "$WORK/scores.sc" \
  --threaded "$WORK/in/threaded_manifest.csv" \
  --out results/processed/af2_deep.csv

echo -e "\n=== deep sampling done $(date -u +%H:%M) ==="
