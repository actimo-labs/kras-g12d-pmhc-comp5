#!/usr/bin/env python
"""Stage 0 -- build the three assay states from 7OW6 with the TCR removed.

States (configs/challenge.yaml:states):
  g12d   chains A (HLA heavy) / B (b2m) / C (VVVGADGVGK)   -- the bind target
  wt     as g12d with peptide D6->G                        -- counter-target 1
  empty  chains A / B only, groove empty                   -- counter-target 2

Three deliberate deviations from docs/PIPELINE.md Stage 0, all forced by the
environment audit (no PyRosetta, no OpenMM, no PDBFixer on this box):

1. `wt` gets NO relax. PIPELINE.md says "local relax only". Asp->Gly is a pure
   atom-subset deletion (drop CB/CG/OD1/OD2, keep N/CA/C/O), so it cannot
   introduce a steric clash -- there is nothing to relieve. It does leave a
   small unfilled cavity where the carboxylate sat, which is physically what
   the real WT complex has to cope with too. Adding atoms would require
   repacking; removing them does not. Recorded as a known approximation.

2. `empty` gets NO groove relax, same reason. PIPELINE.md already flags that
   peptide-free HLA-A*11:01 is physically unstable, so a rigid empty groove is
   an approximation of a target that is itself ill-defined. Do not read the
   `empty` differential as a quantitative number.

3. The PIPELINE.md "bulge-up" secondary state is NOT built here. It needs a
   modelling tool we do not have, and the audit shows the measured peptide
   surface is small enough (see --report) that the question it answers is not
   the bottleneck. Dropped per its own >1h rule.

Run:
  conda run -n boltz python scripts/01_prep_targets.py
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

HEAVY, B2M, PEP = "A", "B", "C"
TCR_CHAINS = {"D", "E"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
# Asp -> Gly: these side-chain atoms are deleted. Gly keeps only backbone.
ASP_SIDECHAIN = {"CB", "CG", "OD1", "OD2"}


def load_clean(cif: Path):
    """7OW6 with waters, ligands and the TCR stripped. Returns a Bio.PDB model."""
    from Bio.PDB import MMCIFParser

    st = MMCIFParser(QUIET=True).get_structure("7ow6", str(cif))
    model = st[0]
    for ch in list(model):
        if ch.id in TCR_CHAINS:
            model.detach_child(ch.id)
            continue
        for res in list(ch):
            if res.id[0] != " ":          # HOH, SO4, anything non-standard
                ch.detach_child(res.id)
    return model


def assert_peptide(model, expected_seq: str, hotspots: dict) -> None:
    """The comp3 round-2 bug: hotspots must index the TARGET chain's residues.

    Chain C in 7OW6 is already author-numbered 1..10, which is exactly peptide
    position p1..p10. Assert that, so a later crop that renumbers chain C fails
    loudly here instead of silently shifting every hotspot by an offset.
    """
    from Bio.PDB.Polypeptide import three_to_index, index_to_one

    ch = model[PEP]
    nums = [r.id[1] for r in ch]
    seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in ch)

    errs = []
    if nums != list(range(1, len(expected_seq) + 1)):
        errs.append(f"chain {PEP} numbering is {nums}, expected 1..{len(expected_seq)}")
    if seq != expected_seq:
        errs.append(f"chain {PEP} sequence is {seq}, expected {expected_seq}")

    # every hotspot/tracked index must resolve to the residue the config claims
    for group in ("primary", "tracked_not_targeted", "anchors_excluded"):
        for p in hotspots.get(group, []):
            if p not in nums:
                errs.append(f"{group} p{p} absent from chain {PEP}")
                continue
            got = index_to_one(three_to_index(ch[p].get_resname()))
            want = expected_seq[p - 1]
            if got != want:
                errs.append(f"{group} p{p}: structure has {got}, peptide has {want}")

    if errs:
        raise SystemExit("HOTSPOT ASSERTION FAILED:\n  " + "\n  ".join(errs))
    print(f"  [assert] chain {PEP} numbering 1..{len(seq)} == p1..p{len(seq)}, seq {seq} OK")
    print(f"  [assert] hotspots {hotspots['primary']} resolve to "
          f"{''.join(expected_seq[p-1] for p in hotspots['primary'])} OK")


def mutate_d6_to_g(model, pos: int = 6) -> None:
    """In-place Asp->Gly by side-chain deletion. No atoms added, so no clash."""
    res = model[PEP][pos]
    if res.get_resname() != "ASP":
        raise SystemExit(f"expected ASP at p{pos}, found {res.get_resname()}")
    for name in sorted(ASP_SIDECHAIN):
        if name in res:
            res.detach_child(name)
    res.resname = "GLY"
    left = sorted(a.get_id() for a in res)
    if set(left) - BACKBONE:
        raise SystemExit(f"p{pos} still has non-backbone atoms: {left}")
    print(f"  [wt] p{pos} ASP->GLY, atoms now {left} (deletion only, no relax needed)")


def write_state(model, chains: list[str], out_pdb: Path) -> None:
    import copy
    from Bio.PDB import PDBIO

    m = copy.deepcopy(model)
    for ch in list(m):
        if ch.id not in chains:
            m.detach_child(ch.id)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(m)
    io.save(str(out_pdb))
    n = sum(1 for ch in m for _ in ch)
    print(f"  wrote {rel(out_pdb, ROOT)}  chains={chains} residues={n}")


def to_cif(pdb: Path, cif: Path) -> None:
    """mmCIF for Boltz-2 templates. Mirrors tfr1-comp4/scripts/make_template_cif.py:
    gemmi must assign entities + label_seq_id or boltz's parse_polymer dies with a
    bare IndexError, and label_asym_id must equal the author chain name so the
    YAML can say template_id: [A, B, C] and mean it."""
    import gemmi

    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    st.assign_label_seq_id()

    # (a) Boltz matches `template_id` against label_asym_id, and gemmi invents
    #     ids like "Axp"/"Bxp"/"Cxp". Force them to the author chain names so the
    #     YAML can say template_id: [A, B, C] and mean it.
    for model in st:
        for ch in model:
            for res in ch:
                res.subchain = ch.name
    chain_names = {c.name for m in st for c in m}
    for ent in st.entities:
        ent.subchains = sorted(set(ent.subchains) & chain_names) or \
            sorted(n for n in chain_names if n == ent.name)

    # (b) gemmi leaves full_sequence empty for a file with no SEQRES records, and
    #     boltz's parse_polymer then indexes sequence[j] off the end and dies with
    #     a bare IndexError. Fill it from the observed residues of each entity.
    for ent in st.entities:
        if ent.entity_type != gemmi.EntityType.Polymer or ent.full_sequence:
            continue
        names: list[str] = []
        for model in st:
            for ch in model:
                if ch.name in ent.subchains:
                    names.extend(r.name for r in ch)
            if names:
                break
        ent.full_sequence = names

    bad = [e.name for e in st.entities
           if e.entity_type == gemmi.EntityType.Polymer and not e.full_sequence]
    if bad:
        raise SystemExit(f"{cif.name}: entities still without a sequence: {bad}")

    st.make_mmcif_document().write_file(str(cif))
    ents = ", ".join(f"{e.name}:{len(e.full_sequence)}" for e in st.entities
                     if e.entity_type == gemmi.EntityType.Polymer)
    print(f"  wrote {rel(cif, ROOT)}  entities[{ents}]")


def sasa_report(model, hotspots: dict, seq: str) -> dict:
    """Per-peptide-position SASA in the TCR-removed state.

    Reported BOTH whole-residue and side-chain-only, because p4/p7/p9 are Gly
    (zero side-chain atoms) and p5 is Ala. A side-chain-only number is
    meaningless at those positions -- their exposed surface is backbone. This is
    the configs/challenge.yaml:hotspots.low_side_chain_note point, in numbers.
    """
    import copy
    from Bio.PDB.SASA import ShrakeRupley

    def compute(chains, level):
        m = copy.deepcopy(model)
        for ch in list(m):
            if ch.id not in chains:
                m.detach_child(ch.id)
        ShrakeRupley().compute(m, level=level)
        return m

    whole_ctx = compute([HEAVY, B2M, PEP], "R")
    whole_iso = compute([PEP], "R")
    atom_ctx = compute([HEAVY, B2M, PEP], "A")
    atom_iso = compute([PEP], "A")

    def sc(m, p):
        return sum(a.sasa for a in m[PEP][p] if a.get_id() not in BACKBONE)

    rows = {}
    for p in range(1, len(seq) + 1):
        w_ctx, w_iso = whole_ctx[PEP][p].sasa, whole_iso[PEP][p].sasa
        s_ctx, s_iso = sc(atom_ctx, p), sc(atom_iso, p)
        n_sc = sum(1 for a in whole_ctx[PEP][p] if a.get_id() not in BACKBONE)
        rows[p] = {
            "aa": seq[p - 1],
            "n_sidechain_atoms": n_sc,
            "sasa_whole_pmhc": round(w_ctx, 1),
            "sasa_whole_isolated": round(w_iso, 1),
            "rel_exposure_whole": round(100 * w_ctx / w_iso, 1) if w_iso else None,
            "sasa_sidechain_pmhc": round(s_ctx, 1),
            "rel_exposure_sidechain": round(100 * s_ctx / s_iso, 1) if s_iso else None,
        }
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", type=Path, default=ROOT / "structures/raw/7OW6.cif")
    ap.add_argument("--outdir", type=Path, default=ROOT / "structures")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(ROOT / "configs/challenge.yaml"))
    g12d_seq = cfg["states"]["g12d"]["peptide"]
    wt_seq = cfg["states"]["wt"]["peptide"]
    hotspots = cfg["hotspots"]

    if hotspots["chain"] != PEP:
        raise SystemExit(f"config hotspots.chain={hotspots['chain']}, this script assumes {PEP}")

    print(f"source: {rel(args.cif, ROOT)}")
    model = load_clean(args.cif)
    print("  chains kept:", [ch.id for ch in model],
          "(TCR", sorted(TCR_CHAINS), "removed, waters + SO4 stripped)")

    print("\nASSERTIONS (comp3 round-2 hotspot-offset bug)")
    assert_peptide(model, g12d_seq, hotspots)

    print("\nSASA -- TCR-removed state, i.e. what a binder actually sees")
    rows = sasa_report(model, hotspots, g12d_seq)
    print(f"  {'pos':>3} {'aa':>3} {'nSC':>3} {'whole':>7} {'rel%':>6} {'sc':>6} {'scRel%':>7}  class")
    for p, r in rows.items():
        grp = ("primary" if p in hotspots["primary"] else
               "tracked" if p in hotspots["tracked_not_targeted"] else
               "anchor" if p in hotspots["anchors_excluded"] else "unclassified")
        scr = "n/a" if r["rel_exposure_sidechain"] is None else f"{r['rel_exposure_sidechain']:.1f}"
        print(f"  {p:>3} {r['aa']:>3} {r['n_sidechain_atoms']:>3} "
              f"{r['sasa_whole_pmhc']:7.1f} {r['rel_exposure_whole']:5.1f}% "
              f"{r['sasa_sidechain_pmhc']:6.1f} {scr:>6}%  {grp}")

    total = sum(r["sasa_whole_pmhc"] for r in rows.values())
    print(f"\n  TOTAL peptide SASA available to a binder: {total:.1f} A^2")
    print("  Ceiling on an AREA-based peptide_contact_fraction, if a binder")
    print("  buried every square angstrom of peptide on offer:")
    for bsa in (700, 900, 1100, 1200):
        print(f"    total interface {bsa:5d} A^2 -> max fraction {total/bsa:.3f}")
    gate = cfg["thresholds"]["peptide_contact_fraction_min"]["value"]
    print(f"  configured gate: {gate}  <-- compare against the ceilings above")

    print("\nSTATES")
    out = args.outdir
    write_state(model, [HEAVY, B2M, PEP], out / "g12d.pdb")
    write_state(model, [HEAVY, B2M], out / "empty.pdb")
    import copy
    wt_model = copy.deepcopy(model)
    mutate_d6_to_g(wt_model)
    write_state(wt_model, [HEAVY, B2M, PEP], out / "wt.pdb")

    # verify the wt peptide really reads as the WT sequence
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    got = "".join(index_to_one(three_to_index(r.get_resname())) for r in wt_model[PEP])
    if got != wt_seq:
        raise SystemExit(f"wt peptide is {got}, config says {wt_seq}")
    print(f"  [assert] wt peptide == {got} == config wt OK")

    print("\nTEMPLATES (mmCIF for Boltz-2)")
    for name in ("g12d", "wt", "empty"):
        to_cif(out / f"{name}.pdb", out / f"{name}.cif")

    meta = {
        "source_pdb": "7OW6",
        "source_resolution_angstrom": 2.64,
        "tcr_chains_removed": sorted(TCR_CHAINS),
        "chains": {"heavy": HEAVY, "b2m": B2M, "peptide": PEP},
        "peptide_sasa_tcr_removed": rows,
        "total_peptide_sasa_angstrom2": round(total, 1),
        "wt_build": "D6->G by side-chain atom deletion; no relax (no PyRosetta/OpenMM on box)",
        "empty_build": "chain C deleted; no relax; rigid approximation, see docstring",
        "bulge_up_state": "not_run",
        "relax_applied": False,
    }
    mp = ROOT / "results/processed/target_prep.json"
    mp.parent.mkdir(parents=True, exist_ok=True)
    json.dump(meta, open(mp, "w"), indent=2)
    print(f"\nwrote {rel(mp, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
