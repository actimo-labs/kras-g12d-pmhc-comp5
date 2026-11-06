#!/usr/bin/env bash
# Does averaging over diffusion samples rescue iptm_margin_wt?
#
# THE PROBLEM this exists to answer. Measured on 8 designs refolded with a
# second seed:
#     margin_wt     mean |seedA - seedB| = 0.186   spread across 32 designs sd = 0.157
#     margin_empty  mean |seedA - seedB| = 0.173   spread across 32 designs sd = 0.322
# So single-sample margin_wt has signal-to-noise 0.84 -- the seed variance of ONE
# design is larger than the spread across ALL designs. One design swung
# -0.353 -> +0.185 on a seed change alone. Ranking on it ranks the RNG.
#
# This matters more than any other pilot result because mutation specificity is
# objective 1, it ranks FIRST in the competition, and docs/PIPELINE.md Stage 4
# rests entirely on iptm_margin_wt.
#
# Two ways out, and this measures the cheap one:
#   (a) average n diffusion samples in ONE run  -- model loads once, ~sqrt(n)
#       noise reduction if the samples are independent. Cheap.
#   (b) n separate seeded runs                  -- n x the cost. Expensive.
#
# If (a) works, Stage 3 stays affordable. If it does not, iptm_margin_wt has to
# be dropped as the objective-1 proxy and the ranking plan needs rethinking.
#
# Runs 12 designs spanning the margin range, 3 states, diffusion_samples=5.
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh
conda activate boltz

N_SAMP="${1:-5}"
echo "=== [$(date -u +%H:%M)] diffusion-sample averaging test, n=${N_SAMP} ==="

SEL=$(python - <<'PY'
import csv
rows=list(csv.DictReader(open('pilot/boltz/main_three_state.csv')))
def g(r,k):
    try: return float(r[k])
    except: return None
rows.sort(key=lambda r: g(r,'iptm_margin_wt_peptide') or 0)
idx=sorted({0,1,2,len(rows)//4,len(rows)//4+1,len(rows)//2,len(rows)//2+1,
            3*len(rows)//4,3*len(rows)//4+1,len(rows)-3,len(rows)-2,len(rows)-1})
print(','.join(rows[i]['candidate_id'] for i in idx))
PY
)
echo "designs: $SEL"

if [ ! -f pilot/boltz_samp/samp5_three_state.csv ]; then
  python -I scripts/05_boltz3.py \
    --fasta-dir pilot/mpnn/selected/neutral/seqs pilot/mpnn/selected/basic/seqs \
    --out pilot/boltz_samp --tag samp5 --only "$SEL" \
    --diffusion-samples "$N_SAMP" --average-samples --seed 99
fi

# Second independent run at the same settings, so we can measure the noise of
# the AVERAGED estimator the same way we measured the single-sample one.
if [ ! -f pilot/boltz_samp/samp5b_three_state.csv ]; then
  python -I scripts/05_boltz3.py \
    --fasta-dir pilot/mpnn/selected/neutral/seqs pilot/mpnn/selected/basic/seqs \
    --out pilot/boltz_samp --tag samp5b --only "$SEL" \
    --diffusion-samples "$N_SAMP" --average-samples --seed 4321
fi

python - <<'PY'
import csv, statistics as st
def load(p):
    try:
        return {r['candidate_id']: r for r in csv.DictReader(open(p))}
    except FileNotFoundError:
        return {}
def fl(r,k):
    try: return float(r[k])
    except: return None

a = load('pilot/boltz_samp/samp5_three_state.csv')
b = load('pilot/boltz_samp/samp5b_three_state.csv')
if not a or not b:
    raise SystemExit('sample-test outputs missing')

print(f"\n{'design':<22} {'A':>8} {'B':>8} {'|d|':>7}")
for key, label in (('iptm_margin_wt_peptide','margin_wt'),
                   ('iptm_margin_empty','margin_empty')):
    d=[]
    print(f"\n--- {label} (averaged over diffusion samples) ---")
    for cid in sorted(set(a) & set(b)):
        x, y = fl(a[cid], key), fl(b[cid], key)
        if x is None or y is None: continue
        d.append(abs(x-y))
        print(f"{cid:<22} {x:+8.3f} {y:+8.3f} {abs(x-y):7.3f}")
    if d:
        single = {'margin_wt':0.186, 'margin_empty':0.173}[label]
        spread = {'margin_wt':0.157, 'margin_empty':0.322}[label]
        print(f"  averaged noise mean |A-B| = {st.mean(d):.3f}")
        print(f"  single-sample noise       = {single:.3f}")
        print(f"  noise reduction           = {single/st.mean(d):.2f}x" if st.mean(d)>0 else "")
        print(f"  signal-to-noise: {spread/single:.2f}x -> {spread/st.mean(d):.2f}x"
              if st.mean(d)>0 else "")
        print(f"  VERDICT: {'usable (S/N > 2)' if st.mean(d)>0 and spread/st.mean(d)>2 else 'STILL NOT USABLE for ranking'}")
PY
echo -e "\n=== sample test done $(date -u +%H:%M) ==="
