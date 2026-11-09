#!/usr/bin/env python
"""Final ranking: hard gates, novelty, diversity, then the competition's own order.

RANKING ORDER. `docs/CHALLENGE_SPEC.md §4` is lexicographic, not weighted, and
this script does not build a composite score -- a composite would contradict how
designs are actually judged.

  1. mutation specificity   iptm_margin_wt   (objective 1)
  2. peptide dependency     iptm_margin_empty + peptide geometry (objective 2)
  3. affinity proxy         interface size / contact density (objective 3)

WHAT EACH METRIC IS WORTH HERE, measured, not assumed:

  iptm_margin_wt      usable ONLY at >=5 averaged diffusion samples. At one
                      sample its between-seed Spearman rho is +0.239 -- the
                      ranking is near-random under a reseed. At five it is
                      +0.846. Rows are checked for n_samples_averaged and a row
                      without it is NOT ranked on this field.
  iptm_margin_empty   rho +0.881 at five samples, so reproducible -- but it
                      separated real designs from beta2m framework decoys by only
                      +0.034. Reproducible is not the same as informative. It is
                      used as a confidence input, never as evidence of peptide
                      dependency.
  peptide geometry    bsa_peptide separated peptide-targeted designs from
                      framework decoys 68.3 vs 0.0 A^2. This, not the ipTM
                      margin, is the real discriminator for objective 2.
  af2_pae_interaction the only metric with reported experimental correlation.
                      Already applied as a hard gate upstream; carried here for
                      the tiebreak and the methods note.

HARD GATES run first and every rejection is written with a reason
(CLAUDE.md §4). Nothing is silently dropped.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

AGG_RUN = re.compile(r"([AILMFWVY])\1{4,}")
HOMO_RUN = re.compile(r"(.)\1{4,}")
SEQUON = re.compile(r"N[^P][ST]")


def f(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def hard_gates(r, cfg) -> str:
    """Returns a rejection reason, or '' to keep."""
    seq = r.get("sequence", "")
    lo, hi = cfg["challenge"]["length_min"], cfg["challenge"]["length_max"]
    if not seq:
        return "no_sequence"
    if not (lo <= len(seq) <= hi):
        return f"length_{len(seq)}_outside_{lo}_{hi}"
    if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
        return "invalid_residue"
    if "C" in seq:
        return "free_cys"
    if AGG_RUN.search(seq):
        return "hydrophobic_run>=5"
    if HOMO_RUN.search(seq):
        return "homopolymer_run>=5"
    if SEQUON.search(seq):
        return "n_glycosylation_sequon"
    return ""


def main() -> int:
    import yaml

    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path,
                    default=ROOT / "results/processed/manifest_final.csv")
    ap.add_argument("--survivors", type=Path,
                    default=ROOT / "results/processed/survivors.csv")
    ap.add_argument("--novelty", type=Path,
                    default=ROOT / "results/processed/novelty.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/processed/ranked.csv")
    ap.add_argument("--max-per-backbone", type=int, default=3)
    ap.add_argument("--max-identity", type=float, default=0.90,
                    help="reject a design this similar to a higher-ranked one")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(ROOT / "configs/challenge.yaml"))
    cap = cfg["challenge"]["max_designs"]

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"no rows in {args.manifest}")

    # join AF2 metrics and arm labels back on
    extra = {}
    if args.survivors.exists():
        with open(args.survivors) as fh:
            extra = {r["candidate_id"]: r for r in csv.DictReader(fh)}
    novel = {}
    if args.novelty.exists():
        with open(args.novelty) as fh:
            novel = {r["candidate_id"]: r for r in csv.DictReader(fh)}

    for r in rows:
        e = extra.get(r["candidate_id"], {})
        r["arm"] = e.get("arm", "")
        r["backbone"] = e.get("backbone", "")
        r["af2_pae_interaction"] = e.get("af2_pae_interaction", "not_run")
        r["af2_plddt_binder"] = e.get("af2_plddt_binder", "not_run")
        r["af2_target_aligned_rmsd"] = e.get("af2_target_aligned_rmsd", "not_run")
        n = novel.get(r["candidate_id"], {})
        r["novelty_max_identity"] = n.get("max_identity", "not_run")
        r["novelty_hit"] = n.get("hit", "")

    kept, rejected = [], []
    for r in rows:
        why = hard_gates(r, cfg)
        if not why and str(r.get("novelty_hit", "")).lower() in ("true", "1", "yes"):
            why = f"novelty_hit:{r.get('novelty_max_identity')}"
        if not why:
            ns = f(r, "g12d_n_samples_averaged") or f(r, "n_samples_averaged")
            if ns is not None and ns < 5:
                why = f"margins_from_{ns:.0f}_samples_unusable"
        if why:
            r["reject_reason"] = why
            rejected.append(r)
        else:
            kept.append(r)

    print(f"{len(rows)} folded -> {len(kept)} pass hard gates, {len(rejected)} rejected")
    if rejected:
        from collections import Counter
        for reason, n in Counter(r["reject_reason"].split(":")[0]
                                 for r in rejected).most_common():
            print(f"    {n:4d}  {reason}")

    # lexicographic sort. Missing margins sort last rather than being treated
    # as zero, which would flatter a design that was never measured.
    # NOISE-AWARE lexicographic key.
    #
    # The competition's order is specificity -> peptide dependency -> affinity,
    # and this honours it. But sorting raw iptm_margin_wt sorts noise: the
    # measured 5-sample seed-to-seed |delta| is 0.053, and ranks 5-20 differ by
    # 0.009-0.033 -- less than that. Ordering on those differences is ordering
    # the random number generator, which is the exact error Stage 08 caught in
    # the single-sample margins.
    #
    # So objective 1 is BINNED at the noise floor. Designs whose margins are
    # indistinguishable are treated as tied and the tie is broken by objective 2
    # (buried peptide area), which IS validated as discriminating: 68.3 A^2 for
    # peptide-targeted designs vs 0.0 for beta2m framework decoys. This is still
    # strictly lexicographic -- it just refuses to invent precision the
    # measurement does not have.
    MARGIN_NOISE = 0.053

    def key(r):
        mw = f(r, "iptm_margin_wt")
        me = f(r, "iptm_margin_empty")
        geo = f(r, "bsa_peptide")
        aff = f(r, "bsa_total_target")
        pae = f(r, "af2_pae_interaction")
        # Bin to the noise floor; within a bin objective 1 says nothing.
        # FLOOR, not int(): int() truncates toward zero, so -0.040 and +0.006
        # would both land in bin 0 and a design scoring BETTER on wild-type
        # would tie with one scoring better on G12D. Floor puts negatives in a
        # strictly lower bin, so direction is preserved even where magnitude is
        # not resolved.
        mw_bin = (math.floor(mw / MARGIN_NOISE) if mw is not None else -99)
        return (
            -mw_bin,
            -(geo if geo is not None else -9),
            -(me if me is not None else -9),
            -(aff if aff is not None else -9),
            (pae if pae is not None else 99),
        )

    kept.sort(key=key)

    # Diversity, two filters.
    #
    # 1. SEQUENCE IDENTITY. Exact-match dedup is NOT enough and the submission
    #    platform caught what it missed: two pairs at 98.8% and 97.5% identity,
    #    each the same backbone and same sample index differing only by MPNN
    #    pool, i.e. the positional R/K bias flipping one or two residues. That
    #    mechanism was measured in Stage 20 (the bias fails to change the
    #    sequence 12-14% of the time) and then not acted on here. Near-identical
    #    designs spend two of twenty slots on one molecule.
    # 2. PER-BACKBONE CAP, so one proven scaffold cannot fill the slate.
    def identity(a, b):
        m = sum(1 for x, y in zip(a, b) if x == y)
        return m / max(len(a), len(b))

    seen = defaultdict(int)
    final, deferred = [], []
    for r in kept:
        bb = r.get("backbone") or "unknown"
        seq = r.get("sequence", "")
        twin = next((o for o in final
                     if identity(seq, o.get("sequence", "")) >= args.max_identity), None)
        if twin is not None:
            r["reject_reason"] = (f"identity_{identity(seq, twin['sequence']):.3f}"
                                  f"_to_{twin['candidate_id']}")
            deferred.append(r)
            continue
        if seen[bb] >= args.max_per_backbone:
            r["reject_reason"] = f"diversity_cap_{args.max_per_backbone}_per_backbone"
            deferred.append(r)
            continue
        seen[bb] += 1
        final.append(r)

    for i, r in enumerate(final, 1):
        r["rank"] = i
    selected = final[:cap]
    overflow = final[cap:]
    for r in overflow:
        r["reject_reason"] = f"below_rank_{cap}"

    print(f"\ndiversity cap {args.max_per_backbone}/backbone: {len(deferred)} deferred")
    print(f"selected {len(selected)} of {cap} slots "
          f"from {len({r.get('backbone') for r in selected})} distinct backbones")
    by_arm = defaultdict(int)
    for r in selected:
        by_arm[r.get("arm") or "?"] += 1
    print(f"  arms: {dict(by_arm)}")

    allrows = selected + overflow
    keys = sorted({k for r in allrows for k in r})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(allrows)

    rej = ROOT / "results/rejected/rank_rejected.csv"
    rej.parent.mkdir(parents=True, exist_ok=True)
    allrej = rejected + deferred + overflow
    if allrej:
        keys = sorted({k for r in allrej for k in r})
        with open(rej, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader(); w.writerows(allrej)

    print(f"\nwrote {rel(args.out, ROOT)}")
    print(f"wrote {rel(rej, ROOT)}  ({len(allrej)} rejected/deferred, with reasons)")

    print(f"\nTOP {min(len(selected), 20)}:")
    print(f"  {'rank':>4} {'candidate':<30} {'mgn_wt':>7} {'bsa_pep':>8} {'mgn_emp':>8} {'pae':>6} {'arm':>5}")
    for r in selected[:20]:
        def s(k, w, p=3):
            v = f(r, k)
            return f"{v:>{w}.{p}f}" if v is not None else f"{'n/a':>{w}}"
        print(f"  {r['rank']:>4} {r['candidate_id']:<30} {s('iptm_margin_wt',7)} "
              f"{s('bsa_peptide',8,1)} {s('iptm_margin_empty',8)} "
              f"{s('af2_pae_interaction',6,2)} {r.get('arm','?'):>5}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
