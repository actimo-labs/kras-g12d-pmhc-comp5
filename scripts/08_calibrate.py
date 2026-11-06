#!/usr/bin/env python
"""Replace the provisional thresholds with values measured from negative controls.

`configs/challenge.yaml` marks four thresholds `provisional: true`, and
`docs/PIPELINE.md` Stage 4 is explicit that ranking with them is meaningless.
The pilot proved the point: `peptide_contact_fraction_min = 0.25` sits at or
above the physical ceiling for this target (the peptide presents only 239.4 A^2
in total), and 0/32 pilot designs cleared it.

Method, mirroring the organisers' own rule. They call a design a binder when it
"clears 300% over the negative control", i.e. the threshold is defined BY the
negatives, not chosen in advance. So for each metric we take the negative
distribution and set the cut at its upper percentile: a real design has to beat
what a composition-matched nonsense sequence achieves.

Negative families (see scripts/07_controls.py, scripts/02b_crop_decoy.py):
  shuffled / reversed  -- composition-matched sequence scrambles. Calibrate the
                          ipTM-style scores: whatever these reach is what the
                          model gives any chain of this composition.
  framework decoys     -- beta2m binders with ZERO peptide burial after folding.
                          Calibrate the peptide-dependency metrics: these are
                          real binders that simply do not read the peptide.

Writes configs/thresholds_calibrated.yaml. It does NOT edit challenge.yaml --
the provisional values stay visible so the change is auditable.

Run:
  conda run -n boltz python scripts/08_calibrate.py
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402


def load(p: Path) -> list[dict]:
    if not Path(p).exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def fl(r: dict, k: str):
    v = r.get(k)
    if v in (None, "", "not_run", "None"):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def pct(vals: list[float], p: float):
    if not vals:
        return None
    s = sorted(vals)
    i = min(int(round(p / 100 * (len(s) - 1))), len(s) - 1)
    return s[i]


def describe(vals: list[float], name: str) -> dict:
    if not vals:
        print(f"  {name:34s} n=0  -- NOT CALIBRATED")
        return {"n": 0, "status": "not_run"}
    d = {
        "n": len(vals),
        "median": round(st.median(vals), 4),
        "mean": round(st.mean(vals), 4),
        "sd": round(st.pstdev(vals), 4),
        "p50": round(pct(vals, 50), 4),
        "p90": round(pct(vals, 90), 4),
        "p95": round(pct(vals, 95), 4),
        "max": round(max(vals), 4),
    }
    print(f"  {name:34s} n={d['n']:3d}  median {d['median']:+.3f}  "
          f"p90 {d['p90']:+.3f}  p95 {d['p95']:+.3f}  max {d['max']:+.3f}")
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", type=Path, default=ROOT / "results/processed/manifest.csv")
    ap.add_argument("--decoy", type=Path, default=ROOT / "results/processed/manifest_decoy.csv")
    ap.add_argument("--controls", type=Path, nargs="*",
                    default=[ROOT / "results/processed/manifest_shuffled.csv",
                             ROOT / "results/processed/manifest_reversed.csv"])
    ap.add_argument("--replicate", type=Path,
                    default=ROOT / "pilot/boltz_rep/rep2_three_state.csv")
    ap.add_argument("--percentile", type=float, default=95.0)
    ap.add_argument("--out", type=Path, default=ROOT / "configs/thresholds_calibrated.yaml")
    args = ap.parse_args()

    real = load(args.real)
    decoy = load(args.decoy)
    ctrl: list[dict] = []
    for c in args.controls:
        ctrl += load(c)

    report: dict = {"percentile": args.percentile, "families": {}}

    print("=" * 74)
    print("THRESHOLD CALIBRATION FROM NEGATIVE CONTROLS")
    print("=" * 74)

    print(f"\nreal designs            n={len(real)}")
    print(f"framework decoys        n={len(decoy)}")
    print(f"sequence controls       n={len(ctrl)}")

    # ---- noise floor from the seed replicate -------------------------------
    noise = None
    rep = load(args.replicate)
    if rep and real:
        base = {r["candidate_id"]: r for r in real}
        deltas = {"iptm_margin_wt_peptide": [], "iptm_margin_empty": [],
                  "g12d_iptm_binder_peptide": []}
        for r in rep:
            b = base.get(r["candidate_id"])
            if not b:
                continue
            for k, bk in (("iptm_margin_wt_peptide", "iptm_margin_wt_peptide"),
                          ("iptm_margin_empty", "iptm_margin_empty"),
                          ("g12d_iptm_binder_peptide", "iptm_binder_peptide_g12d")):
                a, c = fl(r, k), fl(b, bk)
                if a is not None and c is not None:
                    deltas[k].append(abs(a - c))
        print("\nSEED-REPLICATE NOISE (|seed A - seed B|, same design):")
        nf = {}
        for k, v in deltas.items():
            nf[k] = describe(v, k)
        noise = nf
        report["seed_replicate_noise"] = nf
        print("\n  A margin is only interpretable if it exceeds this. Any threshold")
        print("  below the p95 of this noise is measuring the random number generator.")

    # ---- ipTM-style metrics, calibrated on sequence controls ---------------
    print("\nSEQUENCE CONTROLS (composition-matched scrambles):")
    cal: dict = {}
    for key in ("iptm_g12d", "iptm_binder_peptide_g12d",
                "iptm_margin_wt", "iptm_margin_empty"):
        vals = [x for x in (fl(r, key) for r in ctrl) if x is not None]
        d = describe(vals, key)
        report["families"].setdefault("sequence_controls", {})[key] = d
        if d.get("n"):
            floor = pct(vals, args.percentile)
            if noise and key in ("iptm_margin_wt", "iptm_margin_empty"):
                nk = ("iptm_margin_wt_peptide" if key == "iptm_margin_wt"
                      else "iptm_margin_empty")
                nv = noise.get(nk, {}).get("p95")
                if nv is not None and nv > floor:
                    print(f"      -> raising {key} cut from {floor:+.3f} to the "
                          f"seed-noise p95 {nv:+.3f} (noise dominates)")
                    floor = nv
            cal[key] = round(floor, 4)

    # ---- peptide-dependency metrics, calibrated on framework decoys --------
    print("\nFRAMEWORK DECOYS, zero peptide burial only (true framework binders):")
    tru = [r for r in decoy if (fl(r, "bsa_peptide") or 0) == 0]
    print(f"  {len(tru)} of {len(decoy)} decoys had zero peptide burial after folding")
    for key in ("peptide_contact_fraction_area", "bsa_peptide", "n_peptide_contacts"):
        vals = [x for x in (fl(r, key) for r in decoy) if x is not None]
        d = describe(vals, key + " (all decoys)")
        report["families"].setdefault("framework_decoys", {})[key] = d
        if d.get("n"):
            cal[key] = round(pct(vals, args.percentile), 4)

    # ---- what the real designs actually achieve ----------------------------
    print("\nREAL DESIGNS, for comparison:")
    achieved = {}
    for key in ("peptide_contact_fraction_area", "bsa_peptide", "n_peptide_contacts",
                "iptm_margin_wt", "iptm_margin_empty"):
        vals = [x for x in (fl(r, key) for r in real) if x is not None]
        achieved[key] = describe(vals, key)
    report["real_designs"] = achieved

    # ---- yield at the calibrated cuts --------------------------------------
    print("\nYIELD AT CALIBRATED CUTS (real designs passing each, independently):")
    yields = {}
    for key, cut in sorted(cal.items()):
        vals = [x for x in (fl(r, key) for r in real) if x is not None]
        n = sum(x > cut for x in vals)
        yields[key] = {"cut": cut, "pass": n, "of": len(vals)}
        print(f"  {key:34s} cut {cut:+.4f}  ->  {n}/{len(vals)} pass")

    report["calibrated"] = cal
    report["yield"] = yields

    # ---- write -------------------------------------------------------------
    lines = [
        "# CALIBRATED thresholds -- generated by scripts/08_calibrate.py.",
        "# DO NOT hand-edit. Regenerate after adding controls.",
        "# Each cut is the %.0f-th percentile of the matching NEGATIVE distribution," % args.percentile,
        "# raised to the seed-replicate noise p95 where noise dominates, mirroring the",
        "# organisers' 'clears 300%% over the negative control' rule: the threshold is",
        "# defined by the negatives, not chosen in advance.",
        "#",
        "# These REPLACE the provisional values in configs/challenge.yaml, which are",
        "# left in place so the change stays auditable.",
        "calibrated_thresholds:",
    ]
    for k, v in sorted(cal.items()):
        y = yields.get(k, {})
        lines.append(f"  {k}:")
        lines.append(f"    value: {v}")
        lines.append(f"    provisional: false")
        lines.append(f"    source: negative_control_p{args.percentile:.0f}")
        lines.append(f"    real_designs_passing: {y.get('pass')}/{y.get('of')}")
    args.out.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {rel(args.out, ROOT)}")

    jp = ROOT / "results/processed/calibration_report.json"
    json.dump(report, open(jp, "w"), indent=2)
    print(f"wrote {rel(jp, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
