#!/usr/bin/env python
"""Negative controls for threshold calibration (comp2 04b/04c pattern).

`configs/challenge.yaml` marks every threshold `provisional: true` and
`docs/PIPELINE.md` Stage 4 says they MUST be replaced by calibrated values
before ranking, or the ranking is meaningless. The pilot showed why this is not
optional: the `peptide_contact_fraction_min` of 0.25 sits above the measured
physical ceiling for this target, so as written it rejects everything.

Two control families, matching PIPELINE.md:

1. `shuffled` -- each real design's sequence randomly permuted, composition
   held EXACTLY constant. This asks whether the ipTM signal is sequence
   specific or whether any 75-residue chain of the same composition scores the
   same against this groove. Composition is preserved because an ipTM
   difference driven by amino-acid content rather than by sequence would be a
   confound, not a control.

2. `reversed` -- each sequence reversed. Preserves composition AND local
   run-length structure, so it is a stricter control than shuffling: it keeps
   whatever compositional texture the shuffle destroys.

The framework-binder decoy, the third control PIPELINE.md asks for, is built
separately by scripts/02b_crop_decoy.py because it needs its own generation run.

Writes fasta dirs that scripts/05_boltz3.py consumes unchanged.

Run:
  conda run -n boltz python scripts/07_controls.py --from results/processed/manifest.csv
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", type=Path,
                    default=ROOT / "results/processed/manifest.csv")
    ap.add_argument("--out", type=Path, default=ROOT / "pilot/controls")
    ap.add_argument("--kinds", default="shuffled,reversed")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.src)))
    rows = [r for r in rows if r.get("sequence")]
    if args.limit:
        rows = rows[:args.limit]
    if not rows:
        raise SystemExit(f"no sequences in {args.src}")

    rng = random.Random(args.seed)
    kinds = [k.strip() for k in args.kinds.split(",") if k.strip()]

    for kind in kinds:
        d = args.out / kind / "seqs"
        d.mkdir(parents=True, exist_ok=True)
        n = 0
        for r in rows:
            seq = r["sequence"]
            if kind == "shuffled":
                chars = list(seq)
                rng.shuffle(chars)
                new = "".join(chars)
            elif kind == "reversed":
                new = seq[::-1]
            else:
                raise SystemExit(f"unknown control kind {kind}")

            assert sorted(new) == sorted(seq), "composition not preserved"
            cid = f"{kind}_{r['candidate_id']}_s1"
            (d / f"{cid}.fa").write_text(
                f">placeholder_backbone_record\nAAAA\n"
                f">control={kind} source={r['candidate_id']}, sample=1, "
                f"score=0.0, global_score=0.0\n{new}\n")
            n += 1
        print(f"{kind}: {n} control sequences -> {rel(d, ROOT)}")

    print("\nComposition is identical to the parent design in every control, so any")
    print("ipTM difference is attributable to sequence order, not to amino-acid content.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
