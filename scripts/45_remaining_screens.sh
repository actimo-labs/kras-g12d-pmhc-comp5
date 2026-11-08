#!/usr/bin/env bash
# Run the remaining AF2 screens back to back, one process at a time.
#
# STRICTLY SEQUENTIAL, and that is the point. AF2 holds ~61 GB on the 80 GB
# card, so two concurrent screens would OOM. Each screen is checkpointed per
# design, so an interruption costs one design and a re-run resumes.
#
# Deliberately NOT chained with `pgrep` waits. An earlier attempt at that
# deadlocked: pgrep -f matches full command lines, and the harness shell
# wrappers carry the heredoc text of the scripts they create -- including the
# name of the script being waited on -- so each waiter blocked on a process
# containing its own search string (results/RUNLOG.md Stage 11b). This script is
# launched once, after the main screen has finished, and simply runs in order.
#
# Usage: bash scripts/45_remaining_screens.sh [beta,deep,long]
set -uo pipefail
cd "$(dirname "$0")/.."
source /home/HarjasG/miniforge3/etc/profile.d/conda.sh

WHICH="${1:-beta,deep,long}"
log() { echo -e "\n=== [$(date -u +%H:%M)] $* ==="; }

run_screen() {   # name indir scores ckpt outdir resultcsv
  local name="$1" indir="$2" work="$3" out="$4"
  local n; n=$(ls "$indir"/*.pdb 2>/dev/null | wc -l)
  if [ "$n" -eq 0 ]; then log "$name: nothing threaded, skipping"; return 0; fi
  log "$name screen: ${n} designs"
  conda activate af2ig
  python scripts/af2ig/predict_pdb.py \
    -pdbdir "$indir" -outpdbdir "$work/out" \
    -scorefilename "$work/scores.sc" -checkpoint_name "$work/check.point" \
    2>&1 | grep -vE "WARNING:tensorflow|Instructions for updating|^$" | tail -40
  conda activate boltz
  python -I scripts/42_af2_filter.py --scores "$work/scores.sc" \
    --threaded "$indir/threaded_manifest.csv" --out "$out"
}

case ",$WHICH," in *,beta,*)
  run_screen "BETA (topology diversity)" \
    /mnt/comp5data/comp5/af2_beta/in /mnt/comp5data/comp5/af2_beta \
    results/processed/af2_beta.csv ;;
esac

case ",$WHICH," in *,deep,*)
  run_screen "DEEP (24 seqs x 4 pools on 12 proven backbones, temp 0.25)" \
    /mnt/comp5data/comp5/deep/in /mnt/comp5data/comp5/deep \
    results/processed/af2_deep.csv ;;
esac

case ",$WHICH," in *,long,*)
  # The long arm (86-100 aa) may still be generating when the earlier screens
  # finish. MPNN + threading are cheap, so do them here rather than earlier, to
  # pick up every backbone that completed in the meantime.
  NL=$(ls gen/rfdiff_long/L_*.pdb 2>/dev/null | wc -l)
  if [ "$NL" -gt 0 ]; then
    log "LONG arm: ${NL} backbones -> MPNN + thread"
    conda activate mpnn
    [ -d /mnt/comp5data/comp5/mpnn_long/sol_neutral/seqs ] || \
      python -I scripts/04_mpnn.py --rfdiff gen/rfdiff_long \
        --out /mnt/comp5data/comp5/mpnn_long --num-seq 2
    conda activate boltz
    python -I scripts/40_thread_af2.py --all \
      --mpnn /mnt/comp5data/comp5/mpnn_long --backbones gen/rfdiff_long \
      --out /mnt/comp5data/comp5/af2_long/in \
      --pools "sol_neutral,sol_basic,van_neutral,van_basic"
    run_screen "LONG (86-100 aa)" \
      /mnt/comp5data/comp5/af2_long/in /mnt/comp5data/comp5/af2_long \
      results/processed/af2_long.csv
  else
    log "LONG arm: no backbones yet, skipping"
  fi
esac

log "ALL REMAINING SCREENS DONE"
conda activate boltz
python -I scripts/44_hit_backbones.py --out /mnt/comp5data/comp5/deep/hit_backbones_final.txt
