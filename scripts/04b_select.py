#!/usr/bin/env python
"""Pick which MPNN sequences are worth folding, and say why each one was dropped.

Co-folding is the budget bottleneck (docs/PIPELINE.md Stage 3), so this is the
cheap gate in front of it. Selection is deliberately boring: composition
filter, then best MPNN score, then at most `--per-backbone` per backbone so one
lucky scaffold cannot fill the whole batch.

Diversity note: capping per backbone is NOT the same as the TM-score clustering
in Stage 5. It is a crude stand-in that costs nothing, and it only guarantees
scaffold spread, not binding-mode spread.

Rejections are appended to results/rejected/ with a reason, per CLAUDE.md §4 —
never silently dropped.

Run:
  conda run -n boltz python scripts/04b_select.py \
      --mpnn pilot/mpnn --pools neutral,basic --top-per-pool 16
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

AGG_RUN = re.compile(r"([AILMFWVY])\1{4,}")
HOMO_RUN = re.compile(r"(.)\1{4,}")
SEQUON = re.compile(r"N[^P][ST]")
HYD = set("AILMFWVY")
CHG = set("DEKR")


def composition_reject(seq: str) -> str:
    if "C" in seq:
        return "free_cys"
    m = AGG_RUN.search(seq)
    if m:
        return f"hydrophobic_run:{m.group(0)}"
    m = HOMO_RUN.search(seq)
    if m:
        return f"homopolymer_run:{m.group(0)}"
    if SEQUON.search(seq):
        return "n_glycosylation_sequon"
    return ""


def parse_pool(seqs_dir: Path) -> list[dict]:
    out = []
    for fa in sorted(seqs_dir.glob("*.fa")):
        recs, name, buf = [], None, []
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        for hdr, seq in recs[1:]:
            s = seq.split("/")[0]
            if not s:
                continue
            gs = re.search(r"global_score=([\d.]+)", hdr)
            sc = re.search(r"score=([\d.]+)", hdr)
            sm = re.search(r"sample=(\d+)", hdr)
            out.append({
                "backbone": fa.stem,
                "candidate_id": f"{fa.stem}_s{sm.group(1) if sm else len(out)}",
                "sequence": s,
                "length": len(s),
                "mpnn_score": float(sc.group(1)) if sc else float("inf"),
                "mpnn_global_score": float(gs.group(1)) if gs else float("inf"),
                "hydrophobic_frac": round(sum(c in HYD for c in s) / len(s), 3),
                "charged_frac": round(sum(c in CHG for c in s) / len(s), 3),
                "net_charge": sum(c in "KR" for c in s) - sum(c in "DE" for c in s),
            })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpnn", type=Path, default=ROOT / "pilot/mpnn")
    ap.add_argument("--pools", default="neutral,basic")
    ap.add_argument("--top-per-pool", type=int, default=16)
    ap.add_argument("--per-backbone", type=int, default=1)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--tag", default="selected")
    args = ap.parse_args()

    out_root = args.out or (args.mpnn / args.tag)
    pools = [p.strip() for p in args.pools.split(",") if p.strip()]
    all_rej, all_sel = [], []

    for pool in pools:
        seqs = args.mpnn / pool / "seqs"
        if not seqs.exists():
            print(f"pool {pool}: no seqs dir at {seqs}", file=sys.stderr)
            continue
        rows = parse_pool(seqs)

        kept = []
        for r in rows:
            why = composition_reject(r["sequence"])
            if why:
                all_rej.append({**r, "pool": pool, "reject_reason": why,
                                "stage": "04b_select_composition"})
            else:
                kept.append(r)

        kept.sort(key=lambda r: r["mpnn_global_score"])
        seen: dict[str, int] = {}
        chosen = []
        for r in kept:
            n = seen.get(r["backbone"], 0)
            if n >= args.per_backbone:
                all_rej.append({**r, "pool": pool,
                                "reject_reason": f"per_backbone_cap>{args.per_backbone}",
                                "stage": "04b_select_diversity"})
                continue
            seen[r["backbone"]] = n + 1
            chosen.append(r)
        overflow = chosen[args.top_per_pool:]
        chosen = chosen[:args.top_per_pool]
        for r in overflow:
            all_rej.append({**r, "pool": pool,
                            "reject_reason": f"below_top_{args.top_per_pool}_by_mpnn_score",
                            "stage": "04b_select_rank"})

        import statistics as st
        print(f"pool {pool}: {len(rows)} seqs -> {len(kept)} pass composition "
              f"-> {len(chosen)} selected (cap {args.per_backbone}/backbone)")
        if chosen:
            print(f"   mpnn_global_score {min(r['mpnn_global_score'] for r in chosen):.3f}"
                  f"-{max(r['mpnn_global_score'] for r in chosen):.3f}  "
                  f"hydrophobic {st.median([r['hydrophobic_frac'] for r in chosen]):.2f}  "
                  f"charged {st.median([r['charged_frac'] for r in chosen]):.2f}  "
                  f"net {st.median([r['net_charge'] for r in chosen]):+.0f}")

        d = out_root / pool / "seqs"
        d.mkdir(parents=True, exist_ok=True)
        for r in chosen:
            # Namespace the id by pool. Both pools design the SAME backbones, so
            # without this, neutral/p_0_s1.fa and basic/p_0_s1.fa collide on one
            # dict key in scripts/05_boltz3.py:read_fastas -- one pool silently
            # overwrites the other and every mpnn_pool label downstream is wrong.
            r["backbone_id"] = r["candidate_id"]
            r["candidate_id"] = f"{pool}_{r['candidate_id']}"
            (d / f"{r['candidate_id']}.fa").write_text(
                f">placeholder_backbone_record\nAAAA\n"
                f">T=0.1, sample=1, score={r['mpnn_score']}, "
                f"global_score={r['mpnn_global_score']}\n{r['sequence']}\n")
            r["pool"] = pool
        all_sel += chosen

    if all_rej:
        rp = ROOT / "results/rejected/select_rejected.csv"
        rp.parent.mkdir(parents=True, exist_ok=True)
        keys = sorted({k for r in all_rej for k in r})
        with open(rp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader(); w.writerows(all_rej)
        print(f"\n{len(all_rej)} rejected -> {rel(rp, ROOT)}")
        from collections import Counter
        for reason, n in Counter(r["reject_reason"].split(":")[0]
                                 for r in all_rej).most_common():
            print(f"   {n:4d}  {reason}")

    sp = ROOT / "results/processed/selected.csv"
    sp.parent.mkdir(parents=True, exist_ok=True)
    if all_sel:
        keys = sorted({k for r in all_sel for k in r})
        with open(sp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader(); w.writerows(all_sel)
    print(f"\n{len(all_sel)} selected -> {rel(out_root, ROOT)} "
          f"(manifest {rel(sp, ROOT)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
