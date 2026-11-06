#!/usr/bin/env python
"""Stage 1/2 -- ProteinMPNN sequence design on RFdiffusion backbones.

Ported from tfr1-comp4/scripts/03_mpnn.py. The design logic is target-agnostic
and the chain layout is identical (RFdiffusion writes binder=A, target=B), so
only the bias pool changes. Target held as fixed context; binder chain designed.

TWO POOLS, kept separate so we can tell afterwards whether the bias did anything
--------------------------------------------------------------------------------
  neutral -- SolubleMPNN, no bias. The control.
  basic   -- SolubleMPNN with a mild Arg/Lys bias applied ONLY at binder
             positions that project over peptide p6.

WHY POSITIONAL AND NOT GLOBAL BIAS
----------------------------------
configs/challenge.yaml asks for "mild Arg/Lys enrichment at positions
projecting over peptide p6". ProteinMPNN's --bias_AA_jsonl is GLOBAL: it
enriches an amino acid at every position in the chain. Applied globally, a
+1.5 Arg/Lys bias on a 70 aa binder makes the whole protein basic, which
wrecks pI and expression and would confound the very comparison the two pools
exist to make. --bias_by_res_jsonl takes a per-position 21-vector instead, so
the bias lands only where the hypothesis says it should.

The hypothesis itself (R50a/K70a electrostatics, Poole & Karuppiah 2022) is
labelled HYPOTHESIS in the config and stays labelled here. It is not a
validated design rule, and the p6 carboxylate it is aimed at is BURIED --
measured 10.4% side-chain exposure, salt-bridged to HLA R114 (2.95/3.04 A) and
H-bonded to Q70 (3.08 A), see results/processed/target_prep.json. So the bias
is aimed at long-range electrostatic complementarity with a partly occluded
charge, not at a direct salt bridge. Expect a small effect or none.

p6 is located through the .trb residue map, never through a hardcoded offset
(src/comp5/rfmap.py explains why).

Runs under the `mpnn` conda env:
  conda run -n mpnn python scripts/04_mpnn.py --rfdiff pilot/rfd --out pilot/mpnn
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.rfmap import MPNN_ALPHABET, trb_map  # noqa: E402

MPNN = Path(os.environ.get("PROTEINMPNN_DIR",
                           Path.home() / "tnf-track3" / "ProteinMPNN"))
BINDER_CHAIN = "A"
TARGET_CHAIN = "B"
PEPTIDE_CHAIN = "C"          # original chain id of the peptide, pre-renumbering


def chains_of(pdb: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[line[21]] = out.get(line[21], 0) + 1
    return out


def stage_inputs(src: Path, dest: Path, limit: int | None) -> list[Path]:
    """Copy backbones into a clean dir; MPNN's parser globs a whole directory.

    Skips anything whose layout isn't binder=A/target=B -- a silent mismatch
    would design the target instead of the binder.
    """
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)

    pdbs = sorted(p for p in src.glob("*.pdb"))
    if limit:
        pdbs = pdbs[:limit]

    staged, skipped = [], []
    for p in pdbs:
        ch = chains_of(p)
        if set(ch) != {BINDER_CHAIN, TARGET_CHAIN} or ch[BINDER_CHAIN] < 10:
            skipped.append((p.name, ch))
            continue
        shutil.copy(p, dest / p.name)
        if p.with_suffix(".trb").exists():
            shutil.copy(p.with_suffix(".trb"), dest / (p.stem + ".trb"))
        staged.append(p)

    if skipped:
        print(f"skipped {len(skipped)} backbone(s) with unexpected chains:", file=sys.stderr)
        for name, ch in skipped[:5]:
            print(f"  {name}: {ch}", file=sys.stderr)
    return staged


def basic_bias_jsonl(staged_dir: Path, out: Path, radius: float, weight: float) -> dict:
    """Per-position Arg/Lys bias at binder residues projecting over peptide p6.

    ProteinMPNN --bias_by_res_jsonl format:
        {pdb_name: {chain: [[21 floats] * n_residues]}}

    It must be ONE JSON object covering every design, on a single line. The
    loader at ProteinMPNN/protein_mpnn_run.py:145 is

        for json_str in json_list:
            bias_by_res_dict = json.loads(json_str)

    which REASSIGNS on each line, so a one-object-per-line file silently keeps
    only the last design and every other design dies with
    `KeyError: <name>` inside tied_featurize. Caught by the pilot smoke test.
    """
    from Bio.PDB import PDBParser
    import numpy as np

    iR, iK = MPNN_ALPHABET.index("R"), MPNN_ALPHABET.index("K")
    parser = PDBParser(QUIET=True)
    merged: dict = {}
    stats: dict = {}

    for pdb in sorted(staged_dir.glob("*.pdb")):
        trb = pdb.with_suffix(".trb")
        if not trb.exists():
            raise SystemExit(f"{pdb.name}: no .trb, cannot locate p6 without guessing")
        omap = trb_map(trb)
        p6 = [k for k, v in omap.items() if v == (PEPTIDE_CHAIN, 6)]
        if len(p6) != 1:
            raise SystemExit(f"{pdb.name}: p6 maps to {p6}; refusing to guess")
        p6_out = p6[0]

        model = parser.get_structure(pdb.stem, str(pdb))[0]
        p6_res = model[p6_out[0]][p6_out[1]]
        anchor = p6_res["CA"].coord

        binder = list(model[BINDER_CHAIN])
        vec = np.zeros((len(binder), len(MPNN_ALPHABET)), dtype=float)
        near = 0
        for i, res in enumerate(binder):
            if "CA" not in res:
                continue
            if float(np.linalg.norm(res["CA"].coord - anchor)) <= radius:
                vec[i, iR] = weight
                vec[i, iK] = weight
                near += 1
        stats[pdb.stem] = near
        merged[pdb.stem] = {BINDER_CHAIN: vec.tolist()}

    # single object, single line -- see docstring
    out.write_text(json.dumps(merged, separators=(",", ":")) + "\n")
    vals = list(stats.values())
    print(f"  basic bias: +{weight} R/K within {radius} A of p6 CA; "
          f"biased positions per design min={min(vals)} median="
          f"{sorted(vals)[len(vals)//2]} max={max(vals)}")
    if max(vals) == 0:
        print("  WARNING: no binder position is within range of p6 in any design. "
              "The basic pool is then identical to neutral -- report it as such, "
              "do not claim the bias was applied.", file=sys.stderr)
    return stats


def run(cmd: list[str]) -> None:
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rfdiff", type=Path, default=ROOT / "pilot/rfd")
    ap.add_argument("--out", type=Path, default=ROOT / "pilot/mpnn")
    ap.add_argument("--num-seq", type=int, default=4)
    ap.add_argument("--temp", default="0.1")
    ap.add_argument("--ala-bias", type=float, default=-0.6,
                    help="global Ala log-bias; counters ProteinMPNN's poly-Ala habit")
    ap.add_argument("--bias-weight", type=float, default=1.2)
    ap.add_argument("--bias-radius", type=float, default=12.0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--pools", default="neutral,basic")
    ap.add_argument("--seed", type=int, default=37)
    args = ap.parse_args()

    if not MPNN.exists():
        raise SystemExit(f"ProteinMPNN not found at {MPNN}; set PROTEINMPNN_DIR")

    args.out.mkdir(parents=True, exist_ok=True)
    staged_dir = args.out / "_backbones"
    staged = stage_inputs(args.rfdiff, staged_dir, args.limit)
    if not staged:
        print(f"no usable backbones in {args.rfdiff}", file=sys.stderr)
        return 1
    print(f"staged {len(staged)} backbones")

    parsed = args.out / "parsed.jsonl"
    assigned = args.out / "assigned.jsonl"
    run([sys.executable, MPNN / "helper_scripts/parse_multiple_chains.py",
         f"--input_path={staged_dir}", f"--output_path={parsed}"])
    run([sys.executable, MPNN / "helper_scripts/assign_fixed_chains.py",
         f"--input_path={parsed}", f"--output_path={assigned}",
         f"--chain_list={BINDER_CHAIN}"])

    pools = [p.strip() for p in args.pools.split(",") if p.strip()]

    # GLOBAL Ala de-bias, applied to BOTH pools.
    #
    # Measured in the pilot smoke test: at temp 0.1 with no Ala bias, 6 of 8
    # designs carried a poly-Ala run of 5-13 residues and were rejected by the
    # stage-05 composition filter -- a 75% loss before folding. Runs that long
    # are not merely a filter artefact; they are low-complexity surface with no
    # core packing. This is ProteinMPNN's documented poly-Ala habit on
    # backbone-only helical input.
    #
    # Applied to both pools so it does NOT confound the neutral-vs-basic
    # comparison, which is the only thing those pools exist to measure.
    aa_bias_path = args.out / "bias_aa.jsonl"
    aa_bias_path.write_text(json.dumps({"A": args.ala_bias}) + "\n")
    print(f"global Ala bias {args.ala_bias} -> {aa_bias_path.name}, sampling_temp {args.temp}")

    bias_path = args.out / "bias_by_res.jsonl"
    if "basic" in pools:
        print("building positional R/K bias")
        basic_bias_jsonl(staged_dir, bias_path, args.bias_radius, args.bias_weight)

    for pool in pools:
        folder = args.out / pool
        # MPNN appends into seqs/; stale fastas from an earlier, smaller
        # backbone set would silently ride along downstream.
        if folder.exists():
            shutil.rmtree(folder)
        folder.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, MPNN / "protein_mpnn_run.py",
               "--jsonl_path", parsed,
               "--chain_id_jsonl", assigned,
               "--out_folder", folder,
               "--num_seq_per_target", args.num_seq,
               "--sampling_temp", args.temp,
               "--seed", args.seed,
               "--batch_size", 1,
               # Free Cys is an automatic developability reject later and the
               # expression system here is undisclosed (CHALLENGE_SPEC section 6),
               # so disulfide-dependent designs are a risk we decline. Cheaper to
               # forbid now than to discard finished designs at ranking.
               "--omit_AAs", "CX",
               "--bias_AA_jsonl", aa_bias_path,
               "--use_soluble_model"]
        if pool == "basic":
            cmd += ["--bias_by_res_jsonl", bias_path]
        run(cmd)
        n = len(list((folder / "seqs").glob("*.fa")))
        print(f"pool {pool}: {n} fasta files -> {folder/'seqs'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
