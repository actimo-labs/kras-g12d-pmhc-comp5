#!/usr/bin/env python
"""Crop a FRAMEWORK-BINDER DECOY generation target: beta2m + its heavy-chain context.

This builds a deliberate negative control, not a candidate. `docs/PIPELINE.md`
Stage 4 requires a "design deliberately placed on the alpha1/alpha2 helix [that]
must fail the peptide-dependency gate", and Gate 1 criterion 3 cannot be
evaluated without one: it asks whether the g12d-minus-empty ipTM spread is
*wider for high-peptide-contact designs than for framework-dominated ones*,
which is a comparison and needs both arms.

Why beta2m rather than the alpha1/alpha2 helix. On the cropped platform only six
residues sit >18 A from the peptide with SASA >60 A^2, and no four of them form
a patch tight enough to be a credible epitope. beta2m is a better control
anyway: the organisers name it explicitly -- "a binder that reads the HLA or the
beta-2-microglobulin instead of the peptide binds all three targets and ranks
there" -- so a beta2m binder is the exact failure mode the gate must catch, and
it is peptide-independent by construction rather than by argument.

Hotspots B44/B45/B47/B48: max Calpha-Calpha spread 12.5 A, mean 32.5 A from the
nearest peptide atom. Measured, not chosen by eye.

Generation target is beta2m ALONE (99 residues, one contiguous segment). The
heavy chain is deliberately NOT included: retaining heavy-chain residues within
15 A of beta2m shatters it into 32 discontiguous segments, and a 32-break
contigmap is exactly the invented-peptide-bond hazard the main crop was
designed to avoid. Dropping it is safe because the four hotspots sit 12.7-19.6 A
clear of the nearest heavy-chain atom (asserted below), so a binder on that
patch has room without the alpha3 domain in the way. Validation still uses the
full three states, where any clash with the heavy chain would show up.

Writes structures/crop/decoy_crop.pdb + decoy_crop_report.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

B2M, HEAVY, PEP = "B", "A", "C"
DECOY_HOTSPOTS = [44, 45, 47, 48]


def main() -> int:
    import numpy as np
    from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select

    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, default=ROOT / "structures/g12d.pdb")
    ap.add_argument("--out", type=Path, default=ROOT / "structures/crop/decoy_crop.pdb")
    ap.add_argument("--min-heavy-clearance", type=float, default=10.0)
    args = ap.parse_args()

    model = PDBParser(QUIET=True).get_structure("x", str(args.pdb))[0]

    # verify the hotspots really are far from the peptide
    pep = np.array([a.coord for a in model[PEP].get_atoms()])
    dists = {}
    for h in DECOY_HOTSPOTS:
        res = model[B2M][h]
        dists[h] = round(min(float(np.linalg.norm(a.coord - pep, axis=-1).min())
                             for a in res), 1)
    if min(dists.values()) < 20.0:
        raise SystemExit(f"decoy hotspot too close to peptide: {dists}; "
                         f"it would not be a clean framework control")
    print(f"[assert] decoy hotspot distances to peptide: {dists} -- all >20 A OK")

    # beta2m alone is only safe if the patch is clear of the heavy chain
    heavy_xyz = np.array([a.coord for a in model[HEAVY].get_atoms()])
    clear = {h: round(min(float(np.linalg.norm(a.coord - heavy_xyz, axis=-1).min())
                          for a in model[B2M][h]), 1) for h in DECOY_HOTSPOTS}
    if min(clear.values()) < args.min_heavy_clearance:
        raise SystemExit(f"decoy hotspot too close to the heavy chain: {clear}; "
                         f"beta2m-alone generation would place binders into the "
                         f"alpha3 domain")
    print(f"[assert] clearance from heavy chain: {clear} -- all >"
          f"{args.min_heavy_clearance} A OK")

    class Crop(Select):
        def accept_chain(self, chain):
            return chain.id == B2M

        def accept_residue(self, residue):
            return 1

        def accept_atom(self, atom):
            return 1 if atom.get_altloc() in (" ", "A") else 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(model)
    io.save(str(args.out), Crop())

    crop = PDBParser(QUIET=True).get_structure("c", str(args.out))[0]
    kept = {ch.id: sorted(r.id[1] for r in ch) for ch in crop}
    segs = {}
    for cid, nums in kept.items():
        runs, st, prev = [], nums[0], nums[0]
        for n in nums[1:]:
            if n != prev + 1:
                runs.append((st, prev)); st = n
            prev = n
        runs.append((st, prev))
        segs[cid] = runs
        print(f"  chain {cid}: {len(nums)} residues, {len(runs)} segment(s) {runs}")

    if kept[B2M] != list(range(1, 100)):
        raise SystemExit(f"beta2m renumbered to {kept[B2M][0]}..{kept[B2M][-1]}; "
                         f"decoy hotspots would shift")
    print(f"[assert] beta2m numbering 1..99 preserved -> hotspots "
          f"{[f'{B2M}{h}' for h in DECOY_HOTSPOTS]} unshifted OK")

    contigs = f"{B2M}1-99"
    (args.out.parent / "decoy_contigs.txt").write_text(contigs + "\n")
    total = sum(len(v) for v in kept.values())
    print(f"\n  contigs: {contigs}")
    print(f"  target length: {total} residues")

    json.dump({
        "purpose": "framework-binder NEGATIVE CONTROL, not a candidate",
        "epitope": "beta2-microglobulin",
        "hotspot_res": [f"{B2M}{h}" for h in DECOY_HOTSPOTS],
        "hotspot_dist_to_peptide_angstrom": dists,
        "max_ca_ca_spread_angstrom": 12.5,
        "heavy_chain_included": False,
        "heavy_chain_clearance_angstrom": clear,
        "heavy_chain_omitted_because": "15 A context fragments chain A into 32 segments; "
                                      "patch is 12.7-19.6 A clear so context is unnecessary",
        "chains_kept": sorted(kept),
        "residues_kept": {k: len(v) for k, v in kept.items()},
        "segments": {k: [list(s) for s in v] for k, v in segs.items()},
        "contigs": contigs,
        "expected_behaviour": "must FAIL the peptide-dependency gate: "
                              "iptm_margin_empty ~ 0 and no peptide contacts",
    }, open(args.out.parent / "decoy_crop_report.json", "w"), indent=2)
    print(f"  wrote {rel(args.out.parent / 'decoy_crop_report.json', ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
