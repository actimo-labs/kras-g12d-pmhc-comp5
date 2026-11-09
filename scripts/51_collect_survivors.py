#!/usr/bin/env python
"""Collect AF2 survivors from every arm, deduplicate, and cap by AF2 quality.

DEDUPLICATION IS MANDATORY, not hygiene. Measured: 7,200 main-screen designs
contain 6,742 unique sequences, and the 28 main hits contain only 24 unique
sequences. The soluble and vanilla pools never collide, but within a weight set
the +1.2 positional R/K bias fails to flip ProteinMPNN's near-deterministic
argmax 12-14% of the time, so the same sequence arrives under two pool labels.
Twenty submission slots cannot afford to spend two on one molecule.

ARM BALANCING. The deep arm hit at 10.42% against 0.39% for fresh geometry, so
raw AF2 ranking would hand almost every slot to 11 backbones. That is the
monoculture failure docs/PIPELINE.md section 5 warns about, arriving by a
different route. A per-backbone cap is applied here so that no single proven
backbone monopolises the fold budget, and the cheaper diversity arms (beta
topology, long length) keep representation even though their rates are lower.

Ranking inside the cap is by `pae_interaction` ascending -- the published
criterion, the only metric here with reported experimental correlation.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

ARMS = [
    ("main", "results/processed/af2_screen.csv"),
    ("beta", "results/processed/af2_beta.csv"),
    ("deep", "results/processed/af2_deep.csv"),
    ("long", "results/processed/af2_long.csv"),
]


def f(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=90)
    ap.add_argument("--per-backbone", type=int, default=6)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    rows = []
    for arm, path in ARMS:
        p = ROOT / path
        if not p.exists():
            print(f"  {arm}: not present, skipping")
            continue
        with open(p) as fh:
            got = [r for r in csv.DictReader(fh) if r.get("af2_pass") == "True"]
        for r in got:
            r["arm"] = arm
        rows += got
        print(f"  {arm}: {len(got)} passing")

    if not rows:
        raise SystemExit("no survivors in any arm")

    # dedupe by sequence, keeping the best pae and recording the collision
    best: dict[str, dict] = {}
    dupes = 0
    for r in rows:
        s = r.get("sequence", "")
        if not s:
            continue
        cur = best.get(s)
        if cur is None:
            r["duplicate_of"] = ""
            best[s] = r
        else:
            dupes += 1
            if (f(r, "af2_pae_interaction") or 99) < (f(cur, "af2_pae_interaction") or 99):
                r["duplicate_of"] = cur["candidate_id"]
                best[s] = r
            else:
                cur["duplicate_of"] = (cur.get("duplicate_of") or "") + f"{r['candidate_id']};"
    uniq = list(best.values())
    print(f"\n{len(rows)} passing -> {len(uniq)} unique sequences ({dupes} duplicates folded)")

    uniq.sort(key=lambda r: (f(r, "af2_pae_interaction") or 99))

    # per-backbone cap, so 11 deep backbones cannot take every slot
    seen: dict[str, int] = defaultdict(int)
    chosen, capped = [], 0
    for r in uniq:
        bb = r.get("backbone", "?")
        if seen[bb] >= args.per_backbone:
            capped += 1
            continue
        seen[bb] += 1
        chosen.append(r)
        if len(chosen) >= args.top:
            break

    print(f"cap {args.per_backbone}/backbone: {capped} deferred, {len(chosen)} selected")
    by_arm = defaultdict(int)
    for r in chosen:
        by_arm[r["arm"]] += 1
    print(f"  arms: {dict(by_arm)}")
    print(f"  distinct backbones: {len({r.get('backbone') for r in chosen})}")

    args.out.mkdir(parents=True, exist_ok=True)
    for r in chosen:
        cid = r["candidate_id"]
        (args.out / f"{cid}.fa").write_text(
            f">placeholder_backbone_record\nAAAA\n"
            f">T=0.1, sample=1, score=0.0, global_score=0.0\n{r['sequence']}\n")

    man = ROOT / "results/processed/survivors.csv"
    keys = sorted({k for r in chosen for k in r})
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(chosen)
    print(f"\nwrote {len(chosen)} fastas -> {rel(args.out, ROOT)}")
    print(f"wrote {rel(man, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
