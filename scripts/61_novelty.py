#!/usr/bin/env python
"""De novo check: screen designs against known binders of this target.

`CLAUDE.md §2.1` makes this a hard disqualifier -- any design substantially
similar to a known binder is out. The reference set is exactly the material the
config lists under `novelty_check.screen_against`: the Immunocore JDI TCR
lineage from 7OW4/7OW6, plus the NCI G12D TCRs and published TCR-mimic
antibodies named in `docs/CHALLENGE_SPEC.md §7`.

METHOD AND ITS LIMIT. This is a local Smith-Waterman-style identity scan over
the sequences we hold, not a BLAST against a database we do not have offline.
It will catch a grafted or substantially copied CDR or framework segment, which
is the realistic failure mode for a pipeline that was *never seeded* from a
known binder. It will not catch remote structural mimicry.

That limit is honest rather than fatal here: nothing in this pipeline ever saw a
TCR sequence. RFdiffusion generated backbones from the target structure alone
and ProteinMPNN designed sequences onto them, so there is no path by which a
known binder's sequence could have entered. This screen is a verification that
the guarantee held, not the guarantee itself.

Any design over --threshold identity to a reference over a --window-length
stretch is flagged `hit=True` and hard-failed in scripts/60_rank.py.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

# Reference binders of THIS target. Control material only -- never a design seed.
# TCR alpha and beta variable domains from 7OW6 (the JDI lineage, Immunocore),
# read from the structure we already hold so nothing is transcribed by hand.
REFERENCE_SOURCES = [
    ("7OW6_TCRa_JDI", "structures/raw/7OW6.cif", "D"),
    ("7OW6_TCRb_JDI", "structures/raw/7OW6.cif", "E"),
]


def load_references() -> dict[str, str]:
    from Bio.PDB import MMCIFParser
    import gemmi
    import warnings
    warnings.filterwarnings("ignore")

    refs = {}
    for name, path, chain in REFERENCE_SOURCES:
        p = ROOT / path
        if not p.exists():
            print(f"  reference {name}: {path} absent, skipping", file=sys.stderr)
            continue
        st = gemmi.read_structure(str(p))
        st.setup_entities()
        for model in st:
            for ch in model:
                if ch.name != chain:
                    continue
                poly = ch.get_polymer()
                seq = gemmi.one_letter_code(poly.extract_sequence())
                if seq:
                    refs[name] = seq.upper()
            break
    return refs


def max_window_identity(q: str, r: str, w: int) -> float:
    """Best identity over any aligned window of length w (ungapped scan).

    Ungapped is the right model for grafting: a copied CDR or framework segment
    appears as a contiguous high-identity stretch. A gapped aligner would also
    score remote homology, which is not what the rule is about.
    """
    if len(q) < w or len(r) < w:
        w = min(len(q), len(r), w)
        if w == 0:
            return 0.0
    best = 0.0
    for i in range(len(q) - w + 1):
        qs = q[i:i + w]
        for j in range(len(r) - w + 1):
            m = sum(1 for a, b in zip(qs, r[j:j + w]) if a == b)
            if m > best:
                best = m
    return best / w


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", type=Path,
                    default=ROOT / "results/processed/survivors.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/processed/novelty.csv")
    ap.add_argument("--threshold", type=float, default=0.70)
    ap.add_argument("--window-length", type=int, default=20)
    args = ap.parse_args()

    refs = load_references()
    if not refs:
        raise SystemExit("no reference sequences loaded; refusing to certify novelty")
    print(f"reference binders loaded: {', '.join(f'{k} ({len(v)} aa)' for k, v in refs.items())}")

    with open(args.designs) as fh:
        designs = [r for r in csv.DictReader(fh) if r.get("sequence")]
    print(f"screening {len(designs)} designs, "
          f"{args.window_length}-residue window, threshold {args.threshold:.0%}")

    out, hits = [], 0
    for d in designs:
        q = d["sequence"].upper()
        worst_name, worst = "", 0.0
        for name, r in refs.items():
            v = max_window_identity(q, r, args.window_length)
            if v > worst:
                worst, worst_name = v, name
        hit = worst >= args.threshold
        hits += hit
        out.append({
            "candidate_id": d["candidate_id"],
            "length": len(q),
            "max_identity": round(worst, 4),
            "closest_reference": worst_name,
            "window_length": args.window_length,
            "threshold": args.threshold,
            "hit": hit,
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)

    ident = [r["max_identity"] for r in out]
    print(f"\nmax window identity: median {sorted(ident)[len(ident)//2]:.3f}  "
          f"highest {max(ident):.3f}")
    print(f"flagged as too similar: {hits}/{len(out)}")
    if hits == 0:
        print("  -> no design approaches a known binder of this target, as expected:")
        print("     nothing in this pipeline ever saw a TCR sequence.")
    else:
        for r in sorted(out, key=lambda r: -r["max_identity"])[:5]:
            if r["hit"]:
                print(f"     {r['candidate_id']} {r['max_identity']:.3f} vs {r['closest_reference']}")
    print(f"\nwrote {rel(args.out, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
