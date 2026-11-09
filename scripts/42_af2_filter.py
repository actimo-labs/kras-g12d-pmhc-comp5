#!/usr/bin/env python
"""Parse the AF2-initial-guess scorefile and apply the canonical filter.

Thresholds are the published ones (Bennett et al. 2023, `dl_binder_design`),
NOT calibrated from our negatives, and that difference is deliberate:
`pae_interaction < 10` is the only criterion in this pipeline with reported
experimental correlation. Re-deriving it from our own controls would replace a
validated number with an unvalidated one.

`target_aligned_rmsd` is required alongside it. Measured on 32 pilot designs,
median target_aligned_rmsd was 40.7 Å -- most designs fold correctly and then
relocate the binder wholesale. One design (basic_p_5_s1) scored pae_interaction
8.91, a clear pass, with target_aligned_rmsd 39.3: AF2 was confident about an
interface that is not the designed one. Confidence without pose agreement
describes a different complex that happens to score well, so both gates apply.

Emits one row per design with a `af2_pass` boolean and the reason for failure,
so rejected candidates carry a reason rather than vanishing (CLAUDE.md §4).

Run:
  conda run -n boltz python scripts/42_af2_filter.py \
      --scores /mnt/comp5data/comp5/af2_screen/scores.sc \
      --threaded /mnt/comp5data/comp5/af2_screen/in/threaded_manifest.csv
"""
from __future__ import annotations

import argparse
import csv
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

# Published canonical criteria. Do not quietly loosen these to raise yield.
PAE_MAX = 10.0
PLDDT_BINDER_MIN = 80.0
TARGET_RMSD_MAX = 5.0
BINDER_RMSD_MAX = 1.5      # reported, not gated -- see docstring in 41_af2_screen.sh

POOLS = ("sol_neutral", "sol_basic", "van_neutral", "van_basic",
         "neutral", "basic")   # legacy pilot names last


def read_scorefile(p: Path) -> list[dict]:
    rows, hdr = [], None
    with open(p) as fh:
        for ln in fh:
            parts = ln.split()
            if not parts or parts[0] != "SCORE:":
                continue
            if hdr is None:
                hdr = parts[1:]
                continue
            if len(parts[1:]) != len(hdr):
                continue
            rows.append(dict(zip(hdr, parts[1:])))
    return rows


def f(r: dict, k: str):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--threaded", type=Path, default=None)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/processed/af2_screen.csv")
    args = ap.parse_args()

    rows = read_scorefile(args.scores)
    if not rows:
        raise SystemExit(f"no SCORE lines in {args.scores}")

    seqs = {}
    if args.threaded and args.threaded.exists():
        with open(args.threaded) as fh:
            seqs = {r["candidate_id"]: r for r in csv.DictReader(fh)}

    out = []
    for r in rows:
        cid = r.get("description", "")
        pae, pl = f(r, "pae_interaction"), f(r, "plddt_binder")
        trm, brm = f(r, "target_aligned_rmsd"), f(r, "binder_aligned_rmsd")

        why = []
        if pae is None or pae >= PAE_MAX:
            why.append(f"pae_interaction>={PAE_MAX}")
        if pl is None or pl <= PLDDT_BINDER_MIN:
            why.append(f"plddt_binder<={PLDDT_BINDER_MIN}")
        if trm is None or trm >= TARGET_RMSD_MAX:
            why.append(f"target_aligned_rmsd>={TARGET_RMSD_MAX}")

        meta = seqs.get(cid, {})
        out.append({
            "candidate_id": cid,
            "sequence": meta.get("sequence", ""),
            "length": meta.get("length", ""),
            "backbone": meta.get("backbone", ""),
            # Match by FULL pool name, longest first: pool names such as
            # "sol_neutral" contain underscores, so split("_",1)[0] yields "sol"
            # and silently leaves this column empty -- which it did, hiding the
            # whole soluble-vs-vanilla comparison in the first main-screen report.
            "mpnn_pool": next((_p for _p in sorted(POOLS, key=len, reverse=True)
                               if cid.startswith(_p + "_")), ""),
            "af2_pae_interaction": pae,
            "af2_plddt_binder": pl,
            "af2_target_aligned_rmsd": trm,
            "af2_binder_aligned_rmsd": brm,
            "af2_pass": not why,
            "af2_reject_reason": ";".join(why),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)

    passed = [r for r in out if r["af2_pass"]]
    print(f"AF2 screen: {len(out)} designs, {len(passed)} pass "
          f"({100*len(passed)/len(out):.1f}%)")
    print(f"  gates: pae_interaction<{PAE_MAX}, plddt_binder>{PLDDT_BINDER_MIN}, "
          f"target_aligned_rmsd<{TARGET_RMSD_MAX}")

    for k, lbl in (("af2_pae_interaction", "pae_interaction"),
                   ("af2_plddt_binder", "plddt_binder"),
                   ("af2_target_aligned_rmsd", "target_aligned_rmsd"),
                   ("af2_binder_aligned_rmsd", "binder_aligned_rmsd")):
        v = [r[k] for r in out if r[k] is not None]
        if v:
            s = sorted(v)
            print(f"  {lbl:22s} median {st.median(v):7.2f}  best {min(v):7.2f}  "
                  f"p10 {s[int(0.10*(len(s)-1))]:7.2f}")

    from collections import Counter
    print("  failure modes (first listed reason):")
    for reason, n in Counter(r["af2_reject_reason"].split(";")[0]
                             for r in out if r["af2_reject_reason"]).most_common():
        print(f"    {n:5d}  {reason}")

    if passed:
        print("\n  passing designs:")
        passed.sort(key=lambda r: r["af2_pae_interaction"])
        for r in passed[:25]:
            print(f"    {r['candidate_id']:<28} pae {r['af2_pae_interaction']:6.2f}  "
                  f"tgt_rmsd {r['af2_target_aligned_rmsd']:6.2f}  "
                  f"plddt {r['af2_plddt_binder']:5.1f}")

    rej = ROOT / "results/rejected/af2_rejected.csv"
    rej.parent.mkdir(parents=True, exist_ok=True)
    with open(rej, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows([r for r in out if not r["af2_pass"]])

    print(f"\nwrote {rel(args.out, ROOT)}")
    print(f"wrote {rel(rej, ROOT)}  ({len(out)-len(passed)} rejected, with reasons)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
