#!/usr/bin/env python
"""Did averaging over diffusion samples rescue the ipTM margins?

THE BAR IS RANK CORRELATION, NOT SIGNAL-TO-NOISE.

We rank designs; we never use the absolute margin value. So the question is
whether the ORDERING survives a reseed. At one sample per design, measured on
n=32 (results/RUNLOG.md Stage 08):

    iptm_margin_wt     S/N 1.16x   Spearman rho +0.239   <- ranking ~random
    iptm_margin_empty  S/N 2.37x   Spearman rho +0.853   <- ranking stable

S/N alone would have called margin_wt borderline-survivable. The rank
correlation shows it is not. That is why rho is the criterion here.

Compares two independent averaged runs (different seeds, same settings), so the
noise of the AVERAGED estimator is measured exactly as the single-sample noise
was in Stage 08.
"""
from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# single-sample baselines, n=32, results/RUNLOG.md Stage 08
BASE = {
    "iptm_margin_wt_peptide": {"rho": 0.239, "noise": 0.136, "sd": 0.157,
                               "label": "margin_wt"},
    "iptm_margin_empty": {"rho": 0.853, "noise": 0.136, "sd": 0.322,
                          "label": "margin_empty"},
}
RHO_BAR = 0.70


def load(p: Path) -> dict[str, dict]:
    if not p.exists():
        return {}
    with open(p) as fh:
        return {r["candidate_id"]: r for r in csv.DictReader(fh)}


def fl(r: dict, k: str):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def spearman(a: list[float], b: list[float]) -> float:
    if len(a) < 3:
        return float("nan")

    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos
        return r

    ra, rb = rank(a), rank(b)
    n = len(a)
    return 1 - 6 * sum((x - y) ** 2 for x, y in zip(ra, rb)) / (n * (n * n - 1))


def main() -> int:
    a = load(ROOT / "pilot/boltz_samp/samp5_three_state.csv")
    b = load(ROOT / "pilot/boltz_samp/samp5b_three_state.csv")
    if not a or not b:
        print("sample-test outputs missing; nothing to analyse", file=sys.stderr)
        return 1

    common = sorted(set(a) & set(b))
    print(f"two independent averaged runs, {len(common)} designs in common")
    n_samp = fl(next(iter(a.values())), "g12d_n_samples_averaged")
    if n_samp:
        print(f"diffusion samples averaged per prediction: {n_samp:.0f}")

    any_rescued = False
    for key, meta in BASE.items():
        label = meta["label"]
        A, B, d = [], [], []
        for cid in common:
            x, y = fl(a[cid], key), fl(b[cid], key)
            if x is None or y is None:
                continue
            A.append(x); B.append(y); d.append(abs(x - y))
        if not d:
            print(f"\n--- {label}: no data ---")
            continue

        rho = spearman(A, B)
        noise = st.mean(d)
        print(f"\n--- {label}, averaged (n={len(d)}) ---")
        for cid, x, y in zip(common, A, B):
            print(f"  {cid:<24} {x:+8.3f} {y:+8.3f}  |d|={abs(x-y):.3f}")
        print(f"  noise mean|d| : {noise:.3f}   (single-sample {meta['noise']:.3f})")
        print(f"  Spearman rho  : {rho:+.3f}   (single-sample {meta['rho']:+.3f})")
        if noise > 0:
            print(f"  S/N           : {meta['sd']/noise:.2f}x   "
                  f"(single-sample {meta['sd']/meta['noise']:.2f}x)")
        ok = rho > RHO_BAR
        print(f"  VERDICT (rho > {RHO_BAR:.2f}): "
              f"{'USABLE for ranking' if ok else 'NOT USABLE for ranking'}")
        if ok and meta["rho"] <= RHO_BAR:
            any_rescued = True
            print(f"  -> averaging RESCUED {label}. Objective 1 regains a proxy;")
            print(f"     price the extra samples into the Stage 3 budget.")
        elif not ok:
            print(f"  -> averaging did NOT rescue {label} at this sample count.")

    print()
    if any_rescued:
        print("ACTION: re-enable iptm_margin_wt in configs/challenge.yaml ranking.order")
        print("        and re-cost Stage 3 at the higher sample count.")
    else:
        print("ACTION: leave objective 1 UNPREDICTED. Rank on folded-complex geometry")
        print("        (bsa_peptide, hotspot coverage), which does discriminate")
        print("        framework from peptide binders: 68.3 vs 0.0 A^2.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
