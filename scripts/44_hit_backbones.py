#!/usr/bin/env python
"""Collect backbones that produced an AF2 hit, from every completed screen.

Reads at run time rather than hardcoding, so re-running after the beta and
long-binder arms finish picks up their hits automatically.

A "hit" is the full canonical criterion set, matching scripts/42_af2_filter.py:
pae_interaction < 10, plddt_binder > 80, target_aligned_rmsd < 5.

Emits `<backbone_dir>:<backbone_stem>` lines so the caller knows where to find
each PDB -- the three arms live in different directories and their prefixes
(g_, b_, L_) are not globally unique by construction.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

PAE_MAX, PLDDT_MIN, TRMSD_MAX = 10.0, 80.0, 5.0
POOLS = ("sol_neutral", "sol_basic", "van_neutral", "van_basic", "neutral", "basic")

# screen scorefile -> directory holding that screen's backbones
SCREENS = [
    ("/mnt/comp5data/comp5/af2_screen/scores.sc", "gen/rfdiff"),
    ("/mnt/comp5data/comp5/af2_beta/scores.sc", "gen/rfdiff_beta"),
    ("/mnt/comp5data/comp5/af2_long/scores.sc", "gen/rfdiff_long"),
]


def read_scores(p: Path) -> list[dict]:
    rows, hdr = [], None
    if not p.exists():
        return rows
    with open(p) as fh:
        for ln in fh:
            q = ln.split()
            if not q or q[0] != "SCORE:":
                continue
            if hdr is None:
                hdr = q[1:]
                continue
            if len(q[1:]) != len(hdr):
                continue
            rows.append(dict(zip(hdr, q[1:])))
    return rows


def f(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def backbone_of(desc: str) -> str | None:
    for p in sorted(POOLS, key=len, reverse=True):
        if desc.startswith(p + "_"):
            core = desc[len(p) + 1:]
            return core.rsplit("_s", 1)[0] if "_s" in core else core
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--pae-max", type=float, default=PAE_MAX)
    args = ap.parse_args()

    found: dict[str, str] = {}
    total_hits = 0
    for sc, bbdir in SCREENS:
        rows = read_scores(Path(sc))
        if not rows:
            continue
        hits = [r for r in rows
                if (f(r, "pae_interaction") or 99) < args.pae_max
                and (f(r, "plddt_binder") or 0) > PLDDT_MIN
                and (f(r, "target_aligned_rmsd") or 99) < TRMSD_MAX]
        total_hits += len(hits)
        n_new = 0
        for r in hits:
            bb = backbone_of(r["description"])
            if bb and (ROOT / bbdir / f"{bb}.pdb").exists():
                if bb not in found:
                    n_new += 1
                found[bb] = bbdir
        print(f"  {sc}: {len(rows)} scored, {len(hits)} hits, {n_new} new backbones")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(f"{d}:{b}\n" for b, d in sorted(found.items())))
    print(f"\n{total_hits} hits span {len(found)} distinct backbones -> {rel(args.out, ROOT)}")
    for b, d in sorted(found.items()):
        print(f"  {d}/{b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
