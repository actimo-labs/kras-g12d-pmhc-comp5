#!/usr/bin/env python
"""Write the submission CSV and audit it. NEVER submits.

`CLAUDE.md §2.5` and `configs/challenge.yaml:submission.auto_submit=false`:
Harjas submits, this script only prepares and verifies. It writes
`submission/comp5_track3.csv` plus an audit report and a methods note.

COLUMN NAMES ARE UNVERIFIED. `docs/CHALLENGE_SPEC.md §6` lists the submission
CSV columns as UNKNOWN for this round and `configs/challenge.yaml` carries
`columns_verified: false`. Prior rounds used `name,sequence,molecule_class`,
which is what this writes, and the audit prints a loud warning until someone
checks the live form. Guessing a column name silently is how a valid design set
gets rejected on a technicality.

The audit re-derives every gate from the CSV itself rather than trusting the
upstream manifest, so a bug in ranking cannot produce a CSV that passes its own
check.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Verified against the live challenge page 2026-11-09. molecule_class must be
# one of these EXACT strings. We previously wrote "protein_minibinder", which is
# descriptive but not a legal value and would have been rejected. "Protein
# minibinder, 40-100 aa" is the published molecule CATEGORY used for scoring
# stratification -- it is not the molecule_class field.
VALID_MOLECULE_CLASS = {"protein", "nanobody", "scfv", "fab_kappa", "fab_lambda"}


def f(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def main() -> int:
    import yaml

    ap = argparse.ArgumentParser()
    ap.add_argument("--ranked", type=Path,
                    default=ROOT / "results/processed/ranked.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "submission/comp5_track3.csv")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(ROOT / "configs/challenge.yaml"))
    ch = cfg["challenge"]
    cap, lo, hi = ch["max_designs"], ch["length_min"], ch["length_max"]
    cols = cfg["submission"]["columns"]
    verified = cfg["submission"]["columns_verified"]

    with open(args.ranked) as fh:
        rows = [r for r in csv.DictReader(fh)
                if r.get("rank") and str(r["rank"]).isdigit()]
    rows.sort(key=lambda r: int(r["rank"]))
    rows = [r for r in rows if int(r["rank"]) <= cap]
    if not rows:
        raise SystemExit(f"no ranked rows in {args.ranked}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            name = f"comp5_t3_{int(r['rank']):02d}_{r['candidate_id'][:36]}"
            mol = "protein"   # single-chain de novo protein; see VALID_MOLECULE_CLASS
            w.writerow([name, r["sequence"], mol])

    # ---------------- audit, re-derived from the written CSV ----------------
    with open(args.out) as fh:
        wrote = list(csv.DictReader(fh))

    problems, notes = [], []
    if len(wrote) > cap:
        problems.append(f"row count {len(wrote)} exceeds cap {cap}")
    if len(wrote) < cap:
        notes.append(f"row count {len(wrote)} is under the cap of {cap} "
                     f"(allowed, but slots are unused)")

    seqs = [r["sequence"] for r in wrote]
    dup = [s for s, n in Counter(seqs).items() if n > 1]
    if dup:
        problems.append(f"{len(dup)} duplicate sequence(s) in the submission")
    names = [r["name"] for r in wrote]
    if len(set(names)) != len(names):
        problems.append("duplicate design names")

    for r in wrote:
        s = r["sequence"]
        if not (lo <= len(s) <= hi):
            problems.append(f"{r['name']}: length {len(s)} outside [{lo},{hi}]")
        bad = sorted(set(s) - VALID_AA)
        if bad:
            problems.append(f"{r['name']}: invalid residues {bad}")
        if "C" in s:
            notes.append(f"{r['name']}: contains Cys (expression system undisclosed)")
    for r in wrote:
        if r["molecule_class"] not in VALID_MOLECULE_CLASS:
            problems.append(f"{r['name']}: molecule_class {r['molecule_class']!r} not in "
                            f"{sorted(VALID_MOLECULE_CLASS)}")

    if not verified:
        problems.append("submission.columns_verified is FALSE -- the column names "
                        f"{cols} are inherited from prior rounds and NOT checked "
                        "against this round's live form")
    if str(ch.get("track")) != "3":
        problems.append(f"track is {ch.get('track')}, expected 3")

    lens = [len(s) for s in seqs]
    cat1 = sum(1 for n in lens if 40 <= n <= 100)

    report = [
        "# comp5 Track 3 — pre-submission audit",
        "",
        f"- file: `{rel(args.out, ROOT)}`",
        f"- rows: **{len(wrote)}** (cap {cap}, Track {ch.get('track')})",
        f"- columns: `{','.join(cols)}` — **verified: {verified}**",
        f"- length range: {min(lens)}–{max(lens)} aa (bounds {lo}–{hi})",
        f"- in published molecule category 1 (40–100 aa): {cat1}/{len(wrote)}",
        f"- unique sequences: {len(set(seqs))}/{len(wrote)}",
        f"- sha256 of file: `{hashlib.sha256(args.out.read_bytes()).hexdigest()[:16]}`",
        "",
        "## Blocking problems",
        "",
    ]
    report += [f"- **{p}**" for p in problems] or ["- none"]
    report += ["", "## Notes", ""]
    report += [f"- {n}" for n in notes] or ["- none"]
    report += [
        "",
        "## Verdict",
        "",
        ("**NOT READY** — resolve the blocking problems above." if problems
         else "**READY for Harjas to submit.** This script never submits."),
    ]

    ap_path = ROOT / "submission/AUDIT.md"
    ap_path.write_text("\n".join(report) + "\n")

    print("\n".join(report[:14]))
    print()
    for p in problems:
        print(f"  PROBLEM: {p}")
    for n in notes[:5]:
        print(f"  note: {n}")
    print(f"\nwrote {rel(args.out, ROOT)}")
    print(f"wrote {rel(ap_path, ROOT)}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
