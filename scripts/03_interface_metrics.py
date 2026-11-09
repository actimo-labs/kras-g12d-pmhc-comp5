#!/usr/bin/env python
"""Interface metrics for a binder against the pMHC, peptide-vs-MHC resolved.

THE RENUMBERING HAZARD, and how this script refuses to fall for it
------------------------------------------------------------------
RFdiffusion does not preserve the input chain IDs or numbering. Measured on
this target, input `A1-180/0 C1-10` comes back as:

    chain A = the designed binder (1..N)
    chain B = heavy chain AND peptide MERGED, renumbered N+1 .. N+190

So peptide p1..p10 land at B259..B268 for an N=78 binder, and at B251..B260 for
an N=70 binder -- the offset MOVES WITH BINDER LENGTH. Hardcoding it is the
comp3 round-2 bug with a new coat of paint.

This script therefore never guesses. It reads `complex_con_ref_pdb_idx` and
`complex_con_hal_pdb_idx` out of the sibling .trb, which RFdiffusion writes as
an explicit (original -> output) residue map, and inverts it. If the .trb is
missing or the map does not cover the peptide, it refuses to score rather than
emitting a plausible-looking number.

WHY CONTACTS ARE COUNTED ON ALL HEAVY ATOMS, BACKBONE INCLUDED
--------------------------------------------------------------
Peptide p4/p7/p9 are Gly (zero side-chain atoms) and p5 is Ala (one). Measured
side-chain SASA at those positions is 0.0 A^2 by construction. A side-chain
contact metric scores them as untouched no matter how well a binder packs the
groove, which would silently reject exactly the designs we want. See
configs/challenge.yaml:hotspots.low_side_chain_note. Contacts here are
heavy-atom, backbone included, and the backbone/side-chain split is reported
separately so the bias is visible rather than baked in.

Two fractions are emitted because they are NOT interchangeable:
  peptide_contact_fraction_area  -- buried-surface-area fraction. This is the
      definition challenge.yaml's threshold references (it cites Baker's
      "~1200 A^2 buried, ~31% from peptide"). Measured ceiling on this target
      is ~0.20-0.34 depending on total interface size, because the peptide
      only presents 239.4 A^2 in total (results/processed/target_prep.json).
  peptide_contact_fraction_count -- fraction of contacted target RESIDUES that
      are peptide. Different scale entirely; a 10-residue peptide against a
      180-residue platform cannot reach 0.3 here either.

Neither is a prediction of anything. They are geometric descriptors used to
rank and to gate, and the gate value itself must be calibrated from the
negative controls before it means anything.

Run:
  conda run -n boltz python scripts/03_interface_metrics.py --glob 'pilot/rfd/*.pdb'
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402
from comp5.rfmap import trb_map  # noqa: E402  (shared; see module docstring)

BACKBONE = {"N", "CA", "C", "O", "OXT"}
CONTACT_CUTOFF = 4.5           # side-chain-complete models
BACKBONE_CONTACT_CUTOFF = 8.0  # backbone-only models (raw RFdiffusion)


def load_model(pdb: Path):
    from Bio.PDB import PDBParser
    return PDBParser(QUIET=True).get_structure(pdb.stem, str(pdb))[0]


def heavy_atoms(res):
    return [a for a in res if a.element != "H" and a.get_id() != "OXT"]


def has_sidechains(model) -> bool:
    """RFdiffusion writes N/CA/C/O only -- no CB, no side chains.

    BSA and side-chain contact counts computed on such a file are not
    comparable to the Baker ~1200 A^2 / ~31% reference, which is measured on a
    side-chain-complete crystal structure. Detecting this is what stops the
    script emitting a number that looks like the gate metric but is not.
    """
    for ch in model:
        for r in ch:
            for a in r:
                if a.get_id() not in BACKBONE and a.get_id() != "CB":
                    return True
    return False


def per_residue_bsa(model, binder_chain: str, target_chains: list[str]) -> dict:
    """SASA(isolated) - SASA(complex) per target residue. Needs 3 SASA passes."""
    import copy
    from Bio.PDB.SASA import ShrakeRupley

    def subset(chains):
        m = copy.deepcopy(model)
        for ch in list(m):
            if ch.id not in chains:
                m.detach_child(ch.id)
        ShrakeRupley().compute(m, level="R")
        return m

    cplx = subset([binder_chain] + target_chains)
    tgt = subset(target_chains)
    bnd = subset([binder_chain])

    bsa = {}
    for ch in tgt:
        for r in ch:
            d = r.sasa - cplx[ch.id][r.id[1]].sasa
            if d > 0.1:
                bsa[(ch.id, r.id[1])] = d
    binder_bsa = sum(
        r.sasa - cplx[binder_chain][r.id[1]].sasa for r in bnd[binder_chain]
    )
    return bsa, binder_bsa


def score_one(pdb: Path, peptide_chain: str, hotspots: list[int],
              explicit: tuple[str, str] | None = None) -> dict | None:
    """explicit=(binder_chain, peptide_chain) bypasses the .trb for Boltz output."""
    from Bio.PDB import NeighborSearch

    trb = pdb.with_suffix(".trb")
    if explicit is None and not trb.exists():
        return {"design": re.sub(r"_model_\\d+$", "", pdb.stem), "qc_status": "failed",
                "reject_reason": "no .trb; cannot map output numbering to peptide"}

    model = load_model(pdb)
    chains = {ch.id: ch for ch in model}

    if explicit is None:
        # --- RFdiffusion output: resolve everything through the .trb map ---
        omap = trb_map(trb)
        mapped = {c for c, _ in omap}
        binder = [c for c in chains if c not in mapped]
        if len(binder) != 1:
            return {"design": re.sub(r"_model_\\d+$", "", pdb.stem), "qc_status": "failed",
                    "reject_reason": f"expected 1 unmapped binder chain, got {binder}"}
        binder = binder[0]
        target_chains = sorted(mapped & set(chains))
        pep_out = {k for k, v in omap.items() if v[0] == peptide_chain}
        if not pep_out:
            return {"design": re.sub(r"_model_\\d+$", "", pdb.stem), "qc_status": "failed",
                    "reject_reason": f"trb map has no chain {peptide_chain}; refusing to score"}
    else:
        # --- Boltz output: chain ids come from the YAML we wrote, and each
        # chain is numbered 1..n, so there is no renumbering to undo and no
        # .trb to read. The mapping is the identity on the peptide chain.
        binder, pep_chain_id = explicit
        if binder not in chains:
            return {"design": re.sub(r"_model_\\d+$", "", pdb.stem), "qc_status": "failed",
                    "reject_reason": f"binder chain {binder} absent; have {sorted(chains)}"}
        target_chains = sorted(c for c in chains if c != binder)
        pep_out = {(pep_chain_id, r.id[1]) for r in chains[pep_chain_id]} \
            if pep_chain_id in chains else set()
        omap = {(c, r.id[1]): (c, r.id[1]) for c in chains for r in chains[c]}
        if not pep_out:
            # the `empty` state legitimately has no peptide; score the rest
            pass

    sc_present = has_sidechains(model)
    # Backbone-only files (raw RFdiffusion output) get a WIDER cutoff. At 4.5 A
    # on N/CA/C/O alone a genuinely docked binder registers almost nothing --
    # measured 6 atom pairs < 4.5 A but 606 < 8 A on a design sitting 3.55 A off
    # the peptide. 8.0 A between backbone atoms is the standard proxy for
    # "these residues would be in side-chain contact once side chains exist".
    cutoff = CONTACT_CUTOFF if sc_present else BACKBONE_CONTACT_CUTOFF

    binder_atoms = [a for r in chains[binder] for a in heavy_atoms(r)]
    ns = NeighborSearch(binder_atoms)

    contacts: dict[tuple[str, int], dict] = {}
    for cid in target_chains:
        for res in chains[cid]:
            key = (cid, res.id[1])
            n_bb = n_sc = 0
            for a in heavy_atoms(res):
                if ns.search(a.coord, cutoff):
                    if a.get_id() in BACKBONE:
                        n_bb += 1
                    else:
                        n_sc += 1
            if n_bb or n_sc:
                contacts[key] = {"bb": n_bb, "sc": n_sc, "orig": omap[key]}

    pep_contacts = {k: v for k, v in contacts.items() if k in pep_out}

    if sc_present:
        bsa, binder_bsa = per_residue_bsa(model, binder, target_chains)
        pep_bsa = sum(v for k, v in bsa.items() if k in pep_out)
        tot_bsa = sum(bsa.values())
    else:
        # Refuse to emit a BSA fraction from a side-chain-free model. It would
        # look exactly like the gate metric and be off by a large, unknown
        # factor -- CLAUDE.md section 2: never fill a column with a
        # plausible-looking number.
        bsa, binder_bsa, pep_bsa, tot_bsa = {}, None, None, None

    # hotspot satisfaction, in ORIGINAL peptide numbering
    hit = sorted(omap[k][1] for k in pep_contacts)
    hs_hit = sorted(set(hit) & set(hotspots))

    pep_atom_contacts = sum(v["bb"] + v["sc"] for v in pep_contacts.values())
    pep_bb = sum(v["bb"] for v in pep_contacts.values())

    return {
        # strip Boltz's _model_N suffix so this joins to candidate_id
        "design": re.sub(r"_model_\d+$", "", pdb.stem),
        "binder_chain": binder,
        "binder_length": len(chains[binder]),
        "n_target_contact_res": len(contacts),
        "n_peptide_contacts": len(pep_contacts),
        "peptide_positions_contacted": hit,
        "hotspots_contacted": hs_hit,
        "n_hotspots_contacted": len(hs_hit),
        "peptide_atom_contacts": pep_atom_contacts,
        "peptide_contacts_backbone": pep_bb,
        "peptide_contacts_sidechain": pep_atom_contacts - pep_bb,
        "sidechains_present": sc_present,
        "contact_cutoff_angstrom": cutoff,
        "bsa_total_target": round(tot_bsa, 1) if tot_bsa is not None else "not_run",
        "bsa_peptide": round(pep_bsa, 1) if pep_bsa is not None else "not_run",
        "bsa_binder": round(binder_bsa, 1) if binder_bsa is not None else "not_run",
        "peptide_contact_fraction_area": (
            round(pep_bsa / tot_bsa, 4) if tot_bsa else "not_run"),
        "peptide_contact_fraction_count": (
            round(len(pep_contacts) / len(contacts), 4) if contacts else 0.0),
        "qc_status": "scored" if sc_present else "scored_backbone_only",
        "reject_reason": "",
    }


def main() -> int:
    import glob as _glob
    import yaml

    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", required=True)
    ap.add_argument("--binder-chain", default=None,
                    help="set for Boltz output (e.g. X) to bypass the .trb map")
    ap.add_argument("--peptide-chain", default=None,
                    help="peptide chain id in the scored file; defaults to config")
    ap.add_argument("--out", type=Path, default=ROOT / "results/processed/interface_metrics.csv")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(ROOT / "configs/challenge.yaml"))
    pep_chain = cfg["hotspots"]["chain"]
    hotspots = cfg["hotspots"]["primary"]

    files = sorted(Path(p) for p in _glob.glob(args.glob, recursive=True) if not p.endswith(".trb"))
    if not files:
        raise SystemExit(f"no pdbs matched {args.glob}")

    explicit = None
    if args.binder_chain:
        explicit = (args.binder_chain, args.peptide_chain or pep_chain)
        print(f"explicit chain mode: binder={explicit[0]} peptide={explicit[1]} "
              f"(no .trb needed)")
    rows = [r for r in (score_one(f, pep_chain, hotspots, explicit) for f in files) if r]
    scored = [r for r in rows if r["qc_status"].startswith("scored")]

    hdr = (f"{'design':<14} {'len':>4} {'ctc':>4} {'pep':>4} {'hs':>3} "
           f"{'BSAtot':>7} {'BSApep':>7} {'fracA':>6} {'fracN':>6}  positions")
    print(hdr)
    print("-" * len(hdr))
    def fmt(v, w, p=1):
        return f"{v:>{w}.{p}f}" if isinstance(v, (int, float)) else f"{str(v):>{w}}"

    for r in scored:
        print(f"{r['design']:<14} {r['binder_length']:>4} {r['n_target_contact_res']:>4} "
              f"{r['n_peptide_contacts']:>4} {r['n_hotspots_contacted']:>3} "
              f"{fmt(r['bsa_total_target'],7)} {fmt(r['bsa_peptide'],7)} "
              f"{fmt(r['peptide_contact_fraction_area'],7,3)} "
              f"{fmt(r['peptide_contact_fraction_count'],6,3)}  "
              f"{r['peptide_positions_contacted']}")

    for r in rows:
        if not r["qc_status"].startswith("scored"):
            print(f"{r['design']:<14} FAILED: {r['reject_reason']}")

    if scored:
        import statistics as st
        fa = [r["peptide_contact_fraction_area"] for r in scored
              if isinstance(r["peptide_contact_fraction_area"], float)]
        if fa:
            print(f"\nn={len(fa)}  peptide_contact_fraction_area: "
                  f"median {st.median(fa):.3f}  max {max(fa):.3f}  "
                  f"n>=0.25 {sum(x >= 0.25 for x in fa)}  n>0 {sum(x > 0 for x in fa)}")
        else:
            print(f"\nn={len(scored)}  peptide_contact_fraction_area: not_run "
                  f"(backbone-only input; BSA needs side chains -- score after MPNN+fold)")
        nh = [r["n_hotspots_contacted"] for r in scored]
        print(f"  hotspots contacted (of {len(hotspots)}): "
              f"median {st.median(nh):.1f}  max {max(nh)}  "
              f"designs touching peptide at all: {sum(x > 0 for x in [r['n_peptide_contacts'] for r in scored])}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for r in rows for k in r})
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    try:
        shown = rel(args.out.resolve(), ROOT)
    except ValueError:
        shown = args.out
    print(f"\nwrote {shown}  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
