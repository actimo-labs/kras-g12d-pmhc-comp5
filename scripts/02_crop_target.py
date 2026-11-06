#!/usr/bin/env python
"""Crop the pMHC to the alpha1/alpha2 platform + peptide, for RFdiffusion ONLY.

Why crop: RFdiffusion cost scales ~quadratically in total denoised length. The
full state is 385 residues (A 276 + B 99 + C 10); with a 70 aa binder that is
455, which does not fit 3,500 backbones in the Stage 2 block. The peptide-
binding platform is only the alpha1/alpha2 superdomain of the heavy chain.

Why a DOMAIN crop and not the configured 15 A distance shell: measured, the
15 A shell around chain C selects 158 chain-A residues in EIGHT discontiguous
segments plus 8 stray b2m residues. Feeding that to RFdiffusion means an
8-segment contigmap, and any mis-specified gap invents peptide bonds across
segment boundaries -- the exact failure tfr1-comp4/scripts/02_crop_target.py
warns about. The alpha1/alpha2 platform (heavy chain 1-180) is ONE contiguous
segment, is only 22 residues larger than the shell, and strictly contains it:
the 15 A shell spans A3-A179 and the 12 A shell spans A4-A179. Asserted below.

b2m (chain B) is dropped from the GENERATION target. It sits under the
platform; only 8 of its residues come within 15 A of the peptide and none are
groove-facing, so it cannot be part of a peptide-directed interface. It is
retained in full for VALIDATION -- structures/{g12d,wt,empty}.pdb are
unchanged and remain the scoring targets.

Chain C is kept WHOLE and its author numbering 1..10 is preserved, so hotspot
p4/p5/p7/p8 == C4/C5/C7/C8 in the cropped file exactly as in the full file.
That is the comp3 round-2 bug this script asserts away.

Writes:
  structures/crop/g12d_crop.pdb     cropped coords, original chain IDs/numbering
  structures/crop/crop_contigs.txt  contigmap.contigs for the cropped input
  structures/crop/crop_report.json  segments, counts, hotspot verification
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

PLATFORM = (1, 180)      # alpha1 + alpha2 superdomain of the HLA heavy chain
HEAVY, PEP = "A", "C"


def main() -> int:
    from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select

    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, default=ROOT / "structures/g12d.pdb")
    ap.add_argument("--outdir", type=Path, default=ROOT / "structures/crop")
    ap.add_argument("--verify-shell", type=float, default=15.0,
                    help="assert the distance shell at this radius fits inside PLATFORM")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(ROOT / "configs/challenge.yaml"))
    hs = cfg["hotspots"]
    pep_seq = cfg["states"]["g12d"]["peptide"]

    model = PDBParser(QUIET=True).get_structure("x", str(args.pdb))[0]

    # --- assert the domain crop contains the configured distance shell --------
    pep_atoms = list(model[PEP].get_atoms())
    ns = NeighborSearch([a for ch in model for a in ch.get_atoms() if ch.id != PEP])
    shell: dict[str, set[int]] = {}
    for a in pep_atoms:
        for b in ns.search(a.coord, args.verify_shell):
            pr = b.get_parent()
            shell.setdefault(pr.get_parent().id, set()).add(pr.id[1])

    heavy_shell = sorted(shell.get(HEAVY, set()))
    outside = [n for n in heavy_shell if not (PLATFORM[0] <= n <= PLATFORM[1])]
    if outside:
        raise SystemExit(
            f"CROP ASSERTION FAILED: chain {HEAVY} residues {outside} are within "
            f"{args.verify_shell} A of the peptide but fall outside the platform "
            f"crop {PLATFORM}. Widen PLATFORM or the groove is being truncated.")
    print(f"[assert] {args.verify_shell:.0f} A shell on chain {HEAVY} spans "
          f"{heavy_shell[0]}-{heavy_shell[-1]} ({len(heavy_shell)} res), "
          f"inside platform {PLATFORM[0]}-{PLATFORM[1]} OK")
    dropped_b2m = sorted(shell.get("B", set()))
    if dropped_b2m:
        print(f"[note]   dropping {len(dropped_b2m)} b2m residues that fall inside the "
              f"shell ({dropped_b2m[0]}-{dropped_b2m[-1]}); not groove-facing, "
              f"retained in the validation target")

    # --- write the crop -------------------------------------------------------
    class Crop(Select):
        def accept_chain(self, chain):
            return chain.id in (HEAVY, PEP)

        def accept_residue(self, residue):
            ch = residue.get_parent().id
            if ch == PEP:
                return 1
            return 1 if PLATFORM[0] <= residue.id[1] <= PLATFORM[1] else 0

        def accept_atom(self, atom):
            return 1 if atom.get_altloc() in (" ", "A") else 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_pdb = args.outdir / "g12d_crop.pdb"
    io = PDBIO()
    io.set_structure(model)
    io.save(str(out_pdb), Crop())

    # --- verify what landed on disk ------------------------------------------
    crop = PDBParser(QUIET=True).get_structure("c", str(out_pdb))[0]
    kept = {ch.id: sorted(r.id[1] for r in ch) for ch in crop}
    for cid, nums in kept.items():
        runs, st, prev = [], nums[0], nums[0]
        for n in nums[1:]:
            if n != prev + 1:
                runs.append((st, prev)); st = n
            prev = n
        runs.append((st, prev))
        print(f"  chain {cid}: {len(nums)} residues, {len(runs)} segment(s) {runs}")

    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    got = "".join(index_to_one(three_to_index(r.get_resname())) for r in crop[PEP])
    if got != pep_seq:
        raise SystemExit(f"cropped peptide is {got}, expected {pep_seq}")
    if kept[PEP] != list(range(1, len(pep_seq) + 1)):
        raise SystemExit(f"cropped peptide renumbered to {kept[PEP]}; hotspots would shift")
    print(f"[assert] cropped chain {PEP} == {got}, numbering {kept[PEP][0]}..{kept[PEP][-1]} "
          f"preserved -> hotspots unshifted OK")

    hotspot_tokens = [f"{PEP}{p}" for p in hs["primary"]]
    for p in hs["primary"]:
        res = crop[PEP][p]
        aa = index_to_one(three_to_index(res.get_resname()))
        if aa != pep_seq[p - 1]:
            raise SystemExit(f"hotspot {PEP}{p} is {aa}, peptide says {pep_seq[p-1]}")
    print(f"[assert] RFdiffusion hotspot_res {hotspot_tokens} verified against crop OK")

    # contigmap for the cropped target: one contiguous heavy segment + peptide
    n_heavy = len(kept[HEAVY])
    contigs = f"{HEAVY}{kept[HEAVY][0]}-{kept[HEAVY][-1]}/0 {PEP}1-{len(pep_seq)}"
    (args.outdir / "crop_contigs.txt").write_text(contigs + "\n")
    print(f"\n  target contigs: {contigs}")
    print(f"  target length:  {n_heavy + len(pep_seq)} residues "
          f"(was {sum(1 for ch in model for _ in ch)})")

    report = {
        "source": rel(args.pdb, ROOT),
        "crop_method": "domain (alpha1/alpha2 platform) not distance shell -- see docstring",
        "platform_range": list(PLATFORM),
        "chains_kept": sorted(kept),
        "chains_dropped": ["B"],
        "residues_kept": {k: len(v) for k, v in kept.items()},
        "verify_shell_angstrom": args.verify_shell,
        "shell_heavy_span": [heavy_shell[0], heavy_shell[-1]],
        "shell_fits_inside_platform": True,
        "b2m_shell_residues_dropped": dropped_b2m,
        "peptide_numbering_preserved": True,
        "hotspot_res": hotspot_tokens,
        "contigs": contigs,
        "generation_only": True,
        "validation_target": "structures/{g12d,wt,empty}.pdb (full A+B+C, uncropped)",
    }
    json.dump(report, open(args.outdir / "crop_report.json", "w"), indent=2)
    print(f"  wrote {rel(args.outdir / 'crop_report.json', ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
