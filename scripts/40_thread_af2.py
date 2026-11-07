#!/usr/bin/env python
"""Thread MPNN sequences onto their RFdiffusion backbones as AF2-initial-guess inputs.

Threading logic adapted from glp1r-comp3/scripts/11_mpnn_thread.py:thread().
That script ran MPNN itself; this one consumes the MPNN output we already have,
so sequence design is not repeated.

WHY AF2-INITIAL-GUESS IS WORTH ADDING NOW
-----------------------------------------
It is the field-standard filter for de novo binders (Bennett et al. 2023,
`dl_binder_design`), and `pae_interaction < 10` is its canonical hit criterion.
Nothing in the pipeline so far asks the question it asks: *does this sequence
actually fold into this backbone and hold this interface?* Boltz re-docks the
binder and reports confidence in whatever pose it finds, which is a different
and weaker question.

It also supplies the ORTHOGONAL SECOND MODEL that `docs/PIPELINE.md` wanted from
Chai-1. AF2 and Boltz-2 are genuinely different architectures with different
training sets, so agreement between them means more than two Boltz seeds
agreeing. Dropping Chai-1 therefore costs less than it appeared to.

NO FASTRELAX. Upstream `dl_binder_design` runs Rosetta FastRelax between MPNN
and AF2. PyRosetta is not installed (Stage 00 audit), so that step is omitted.
This is a real deviation and belongs in the methods "not run" list -- it is
recorded here rather than glossed.

GEOMETRY. RFdiffusion emits the binder as poly-GLY with N/CA/C/O and no CB, so
only those four atoms are written for the binder and the residue is renamed to
the designed identity. AF2's parser maps atoms by name into its 37-atom
representation and masks what it cannot find, which is exactly what "initial
guess" wants for the chain being predicted. The target keeps its full side
chains because it is templated, not predicted.

Binder MUST be chain A and target chain B -- the RFdiffusion output layout.
AF2-ig supports at most 2 chains, which is why the CROPPED target is used here
and the full three-chain pMHC is left to Boltz.

Run:
  conda run -n boltz python scripts/40_thread_af2.py \
      --mpnn pilot/mpnn --backbones pilot/rfd --out pilot/af2_in
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

BACKBONE_ATOMS = {"N", "CA", "C", "O"}
ONE2THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}


def kabsch(P, Q):
    """Rigid transform (R, t) mapping P onto Q, both (N,3)."""
    import numpy as np
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def target_with_sidechains(backbone_pdb: Path, crop_pdb: Path):
    """Crop target atom lines moved into the RFdiffusion output frame.

    RFdiffusion RECENTRES the complex -- measured 46.2 Å mean CA displacement
    between its output target and the crop it was given -- and it writes only
    N/CA/C/O for every chain, target included. AF2-initial-guess wants the
    target templated WITH side chains, so the backbone-only output target is
    not good enough.

    Rather than move the binder back (which would need the inverse transform and
    risks drift), the side-chain-bearing crop is superposed ONTO the output
    frame using the CA correspondence recorded in the .trb. The binder therefore
    stays exactly where RFdiffusion put it, and the target gains its real side
    chains. Fails loudly if the fit is poor.
    """
    import numpy as np
    sys.path.insert(0, str(ROOT / "src"))
    from comp5.rfmap import trb_map

    trb = backbone_pdb.with_suffix(".trb")
    omap = trb_map(trb)

    def ca_index(pdb):
        out = {}
        for ln in open(pdb):
            if ln[:4] == "ATOM" and ln[12:16].strip() == "CA" and ln[16] in " A":
                out[(ln[21], int(ln[22:26]))] = np.array(
                    [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        return out

    out_ca, crop_ca = ca_index(backbone_pdb), ca_index(crop_pdb)
    P, Q = [], []
    for out_key, orig in omap.items():
        if out_key in out_ca and orig in crop_ca:
            P.append(crop_ca[orig]); Q.append(out_ca[out_key])
    if len(P) < 20:
        raise ValueError(f"{backbone_pdb.name}: only {len(P)} CA pairs for superposition")
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean()))
    if rmsd > 0.5:
        raise ValueError(f"{backbone_pdb.name}: target superposition RMSD {rmsd:.2f} A")

    moved = []
    for ln in open(crop_pdb):
        if ln[:4] != "ATOM" or ln[16] not in " A":
            continue
        xyz = R @ np.array([float(ln[30:38]), float(ln[38:46]),
                            float(ln[46:54])]) + t
        moved.append(ln[:21] + "B" + ln[22:30]
                     + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + ln[54:])
    return moved, rmsd


def thread(backbone_pdb: Path, seq: str, out_pdb: Path,
           crop_pdb: Path | None = None) -> float:
    """Binder chain A carrying `seq`; target chain B with side chains if crop given."""
    binder_res, lines = [], []
    for ln in open(backbone_pdb):
        if ln[:4] != "ATOM":
            continue
        lines.append(ln)
        if ln[21] == "A":
            n = int(ln[22:26])
            if not binder_res or binder_res[-1] != n:
                binder_res.append(n)

    if len(binder_res) != len(seq):
        raise ValueError(f"{backbone_pdb.name}: backbone has {len(binder_res)} "
                         f"binder residues, sequence has {len(seq)}")
    aamap = {n: ONE2THREE[a] for n, a in zip(binder_res, seq)}

    tgt_lines, rmsd = (None, float("nan"))
    if crop_pdb is not None:
        tgt_lines, rmsd = target_with_sidechains(backbone_pdb, crop_pdb)

    serial = 0
    with open(out_pdb, "w") as fh:
        for ln in lines:
            if ln[21] != "A":
                if tgt_lines is not None:
                    continue           # replaced by the side-chain version below
                serial += 1
                fh.write(ln[:6] + f"{serial:5d}" + ln[11:])
                continue
            if ln[12:16].strip() not in BACKBONE_ATOMS:
                continue
            serial += 1
            fh.write(ln[:6] + f"{serial:5d}" + ln[11:17]
                     + aamap[int(ln[22:26])] + ln[20:])
        fh.write("TER\n")
        if tgt_lines is not None:
            for ln in tgt_lines:
                serial += 1
                fh.write(ln[:6] + f"{serial:5d}" + ln[11:])
            fh.write("TER\n")
        fh.write("END\n")
    return rmsd


def read_selected(mpnn_root: Path, pools: list[str],
                  use_all: bool = False) -> dict[str, tuple[str, str]]:
    """candidate_id -> (backbone_stem, sequence).

    use_all=False reads the pre-selected fastas (one per backbone per pool).
    use_all=True reads EVERY MPNN sequence.

    Screening all of them is the correct protocol and the pilot showed why.
    Pre-selection ranked by MPNN global_score, and MPNN score does not predict
    AF2 success: of 32 pre-selected designs only 1 cleared the canonical AF2
    criteria, and the designs our calibrated Boltz gates liked best scored
    pae_interaction 25-28. Picking one sequence per backbone on a metric that
    does not correlate with the filter throws away most of the diversity MPNN
    produced, for no gain. AF2-ig at ~10 s/design makes screening all of them
    affordable.
    """
    out: dict[str, tuple[str, str]] = {}
    for pool in pools:
        d = (mpnn_root / pool / "seqs") if use_all else (mpnn_root / "selected" / pool / "seqs")
        if not d.exists():
            continue
        for fa in sorted(d.glob("*.fa")):
            recs, name, buf = [], None, []
            for ln in fa.read_text().splitlines():
                if ln.startswith(">"):
                    if name is not None:
                        recs.append((name, "".join(buf)))
                    name, buf = ln[1:], []
                elif ln.strip():
                    buf.append(ln.strip())
            if name is not None:
                recs.append((name, "".join(buf)))
            for _hdr, seq in recs[1:]:
                s = seq.split("/")[0]
                if not s:
                    continue
                if use_all:
                    m = re.search(r"sample=(\d+)", _hdr)
                    cid = f"{pool}_{fa.stem}_s{m.group(1) if m else len(out)}"
                else:
                    cid = fa.stem
                # <pool>_<backbone>_s<N> -> <backbone>
                core = cid.split("_", 1)[1] if cid.split("_", 1)[0] in pools else cid
                bb = re.sub(r"_s\d+$", "", core)
                out[cid] = (bb, s)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpnn", type=Path, default=ROOT / "pilot/mpnn")
    ap.add_argument("--backbones", type=Path, default=ROOT / "pilot/rfd")
    ap.add_argument("--out", type=Path, default=ROOT / "pilot/af2_in")
    ap.add_argument("--pools", default="neutral,basic")
    ap.add_argument("--crop", type=Path, default=ROOT / "structures/crop/g12d_crop.pdb",
                    help="side-chain-bearing target, superposed into the output frame")
    ap.add_argument("--no-sidechain-target", action="store_true")
    ap.add_argument("--all", action="store_true",
                    help="thread EVERY MPNN sequence, not the pre-selected subset")
    args = ap.parse_args()

    pools = [p.strip() for p in args.pools.split(",") if p.strip()]
    sel = read_selected(args.mpnn, pools, use_all=args.all)
    if not sel:
        where = "" if args.all else "/selected"
        raise SystemExit(f"no sequences under {args.mpnn}{where}")

    args.out.mkdir(parents=True, exist_ok=True)
    rows, missing, bad = [], [], []
    for cid, (bb, seq) in sorted(sel.items()):
        pdb = args.backbones / f"{bb}.pdb"
        if not pdb.exists():
            missing.append((cid, bb))
            continue
        try:
            r = thread(pdb, seq, args.out / f"{cid}.pdb",
                       None if args.no_sidechain_target else args.crop)
        except ValueError as e:
            bad.append(str(e))
            continue
        rows.append({"candidate_id": cid, "backbone": bb, "sequence": seq,
                     "length": len(seq), "target_superpose_rmsd": round(r, 3)})

    for m in missing[:5]:
        print(f"  missing backbone for {m[0]}: {m[1]}.pdb", file=sys.stderr)
    for b in bad[:5]:
        print(f"  {b}", file=sys.stderr)

    mf = args.out / "threaded_manifest.csv"
    with open(mf, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["candidate_id", "backbone", "sequence",
                                           "length", "target_superpose_rmsd"])
        w.writeheader(); w.writerows(rows)

    print(f"threaded {len(rows)} designs -> {rel(args.out, ROOT)}")
    if missing:
        print(f"  {len(missing)} skipped: backbone not found")
    if bad:
        print(f"  {len(bad)} skipped: length mismatch")
    print(f"  manifest {rel(mf, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
