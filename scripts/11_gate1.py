#!/usr/bin/env python
"""Evaluate the four Gate 1 criteria from docs/PIPELINE.md, from measured pilot data.

Criterion 3 is the one that needs the decoy arm: "the ipTM spread between g12d
and empty is wider for high-peptide-contact designs than for framework-dominated
ones." That is a two-sample comparison. The peptide-targeted arm supplies the
first group and the beta2m decoys supply the second.

Prints a verdict per criterion and exits non-zero if any fails, so the pilot
cannot be reported as passing by accident.

Run:
  conda run -n boltz python scripts/11_gate1.py
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load(p: Path) -> list[dict]:
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def f(row: dict, key: str):
    v = row.get(key)
    if v in (None, "", "not_run", "None"):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def summarise(vals: list[float]) -> str:
    if not vals:
        return "n=0"
    return (f"n={len(vals)} median={st.median(vals):+.3f} "
            f"mean={st.mean(vals):+.3f} min={min(vals):+.3f} max={max(vals):+.3f}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", type=Path, default=ROOT / "results/processed/manifest.csv")
    ap.add_argument("--decoy", type=Path,
                    default=ROOT / "results/processed/manifest_decoy.csv")
    ap.add_argument("--timing", type=Path, default=ROOT / "pilot/boltz/main_timing.json")
    ap.add_argument("--gate", type=float, default=0.25)
    args = ap.parse_args()

    main_rows, decoy_rows = load(args.main), load(args.decoy)
    results: dict[str, bool] = {}

    print("=" * 72)
    print("GATE 1 — docs/PIPELINE.md Stage 1")
    print("=" * 72)

    # ---- 1. end to end -----------------------------------------------------
    ok1 = bool(main_rows)
    folded = sum(1 for r in main_rows if r.get("qc_status") == "ok")
    print(f"\n[1] Pipeline runs end to end, no manual intervention")
    print(f"    manifest rows: {len(main_rows)} main, {len(decoy_rows)} decoy; "
          f"{folded} folded ok")
    results["end_to_end"] = ok1 and folded > 0

    # ---- 2. peptide contact ------------------------------------------------
    fa = [f(r, "peptide_contact_fraction_area") for r in main_rows]
    fa = [x for x in fa if x is not None]
    bsap = [f(r, "bsa_peptide") for r in main_rows]
    bsap = [x for x in bsap if x is not None]
    npc = [f(r, "n_peptide_contacts") for r in main_rows]
    npc = [x for x in npc if x is not None]
    n_pass = sum(x >= args.gate for x in fa)
    print(f"\n[2] Some designs reach peptide_contact_fraction >= {args.gate}")
    print(f"    area fraction : {summarise(fa)}")
    print(f"    bsa_peptide A^2: {summarise(bsap)}")
    print(f"    n_peptide_contacts: {summarise(npc)}")
    print(f"    designs >= {args.gate}: {n_pass}/{len(fa)}")
    print(f"    NOTE measured ceiling on this target is 0.199-0.342 depending on")
    print(f"         total interface size (peptide presents only 239.4 A^2).")
    results["peptide_contact"] = n_pass > 0

    # ---- 3. the differential discriminates ---------------------------------
    def margins(rows):
        return [x for x in (f(r, "iptm_margin_empty") for r in rows) if x is not None]

    m_main, m_decoy = margins(main_rows), margins(decoy_rows)
    print(f"\n[3] Three-state differential discriminates at all")
    print(f"    iptm_margin_empty, peptide-targeted : {summarise(m_main)}")
    print(f"    iptm_margin_empty, beta2m framework : {summarise(m_decoy)}")

    ok3 = False
    if m_main and m_decoy:
        diff = st.median(m_main) - st.median(m_decoy)
        print(f"    separation (median main - median decoy): {diff:+.3f}")
        ok3 = diff > 0
        if not ok3:
            print("    -> framework decoys are NOT separated from peptide-targeted")
            print("       designs. The empty-state gate does not discriminate as built.")
    else:
        print("    -> cannot evaluate: one arm is empty")
    results["differential_discriminates"] = ok3

    # correlation between peptide contact and the margin, within the main arm
    pairs = [(f(r, "peptide_contact_fraction_area"), f(r, "iptm_margin_empty"))
             for r in main_rows]
    pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
    if len(pairs) >= 4:
        hi = [b for a, b in pairs if a >= st.median([x for x, _ in pairs])]
        lo = [b for a, b in pairs if a < st.median([x for x, _ in pairs])]
        if hi and lo:
            print(f"    within main arm, high vs low peptide contact:")
            print(f"      high-contact half: {summarise(hi)}")
            print(f"      low-contact  half: {summarise(lo)}")
            print(f"      separation: {st.median(hi) - st.median(lo):+.3f}")

    # ---- 4. throughput -----------------------------------------------------
    print(f"\n[4] Measured throughput fits the Stage 3 matrix")
    ok4 = False
    if args.timing.exists():
        t = json.load(open(args.timing))
        tot = sum(t["timings_sec"].values())
        n = t["n_designs"] * len(t["states"])
        per = tot / max(n, 1)
        print(f"    measured: {tot:.0f}s for {n} predictions = {per:.1f}s/prediction")
        # Stage 3: 600 designs x 3 states x 2 models, Boltz is one of the two
        boltz_n = 600 * 3
        boltz_h = boltz_n * per / 3600
        print(f"    Stage 3 Boltz-2 leg (600 designs x 3 states) = "
              f"{boltz_n} predictions = {boltz_h:.1f} GPU-h")
        print(f"    plus a triage pass (g12d only, ~1500 designs) = "
              f"{1500 * per / 3600:.1f} GPU-h")
        print(f"    docs/PIPELINE.md budgets ~20 GPU-h for Stage 3 TOTAL, and that")
        print(f"    must also cover the Chai-1 leg, which is not installed.")
        ok4 = boltz_h < 20
        if not ok4:
            print(f"    -> the Boltz leg ALONE exceeds the whole Stage 3 budget.")
    else:
        print(f"    no timing file at {args.timing}")
    results["throughput"] = ok4

    # ---- verdict -----------------------------------------------------------
    print("\n" + "=" * 72)
    for k, v in results.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    allok = all(results.values())
    print(f"\nGATE 1: {'PASS' if allok else 'FAIL — see above'}")
    print("=" * 72)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
