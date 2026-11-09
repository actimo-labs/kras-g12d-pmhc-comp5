#!/usr/bin/env python3
"""
AF2-initial-guess for binder/target complexes, without PyRosetta.

Drop-in replacement for dl_binder_design/af2_initial_guess/predict.py in its
-pdbdir mode. Silent-file mode is NOT supported (it is the only part that
genuinely needed Rosetta). See runlog/00-environment-audit.md finding F3 and
make_af2_util_nopyro.py for why this fork exists.

What is preserved exactly from upstream:
  * model_1_ptm, num_ensemble 1, max_extra_msa 5, max_msa_clusters 5
  * initial_guess enabled, recycles configurable (default 3)
  * single-sequence MSA (the binder is de novo; there is no MSA to use)
  * target residues templated, binder residues predicted
  * chainbreak detection via amide distance, residue_index truncation
  * score definitions: pae_interaction is the mean of the two off-diagonal
    binder/target blocks; plddt_binder is the mean over binder residues only

What differs:
  * input PDBs are normalised here (see normalise_pdb) rather than by a
    pose round-trip. This reproduces what pose_from_pdb -> dump_pdb did:
    standard residues only, no hydrogens/altlocs/HETATM, residues renumbered
    contiguously across chains, TER between chains.
  * output PDBs are written directly with per-residue pLDDT in the B-factor
    column and binder/target as separate chains.

Chain convention, same as upstream: THE BINDER MUST BE THE FIRST CHAIN.

Usage:
  python predict_pdb.py -pdbdir in/ -outpdbdir out/ -scorefilename out.sc
"""

import argparse
import os
import sys
import glob
from timeit import default_timer as timer

import numpy as np

# --- locate the bundled alphafold package and our de-pyrosetta'd util --------
HERE = os.path.dirname(os.path.abspath(__file__))
DBD = os.path.expanduser("~/tnf-track3/dl_binder_design/af2_initial_guess")
sys.path.insert(0, DBD)   # alphafold/
sys.path.insert(0, HERE)  # af2_util_nopyro

import jax
import jax.numpy as jnp
from jax.lib import xla_bridge

# --- biopython shim ----------------------------------------------------------
# alphafold.data.pipeline transitively imports alphafold.data.templates, which
# does `from Bio.Data import SCOPData`. Biopython removed SCOPData in 1.80; the
# table moved to Bio.Data.PDBData.protein_letters_3to1_extended. We never run
# template search (initial-guess templating is done from the input structure in
# af2_util.generate_template_features), so rather than pin biopython back to
# 1.79 -- which has no wheel for python 3.11 -- we satisfy the import.
import types as _types
try:
    from Bio.Data import SCOPData as _  # noqa: F401
except ImportError:
    import Bio.Data as _BioData
    from Bio.Data.PDBData import protein_letters_3to1_extended as _p3to1
    _shim = _types.ModuleType("Bio.Data.SCOPData")
    _shim.protein_letters_3to1 = dict(_p3to1)
    sys.modules["Bio.Data.SCOPData"] = _shim
    _BioData.SCOPData = _shim

from alphafold.common import residue_constants, protein, confidence
from alphafold.data import pipeline
from alphafold.model import data, config, model

import af2_util_nopyro as af2_util


# ============================ PDB normalisation ==============================

AA3 = set(residue_constants.restype_3to1.keys())


def normalise_pdb(in_fn, out_fn):
    """
    Reproduce the pyrosetta pose_from_pdb -> dump_pdb normalisation.

    Keeps ATOM records for the 20 standard residues (MSE promoted to MET,
    its SE handled downstream as SD), drops hydrogens, keeps only the first
    altloc, renumbers residues contiguously from 1 across all chains, and
    emits TER between chains.

    Returns (n_chains, [chain_lengths]).
    """
    chains = []          # [(chain_id, [ (resname, [atom_lines...]) ])]
    cur_chain = None
    cur_res_key = None

    for line in open(in_fn):
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        resname = line[17:20].strip()
        if resname == "MSE":
            resname = "MET"
        elif rec == "HETATM" or resname not in AA3:
            continue

        altloc = line[16]
        if altloc not in (" ", "A"):
            continue

        element = line[76:78].strip() or line[12:16].strip()[:1]
        if element == "H":
            continue

        chain_id = line[21]
        res_key = (chain_id, line[22:27])  # resseq + icode

        if cur_chain is None or chain_id != cur_chain[0]:
            cur_chain = (chain_id, [])
            chains.append(cur_chain)
            cur_res_key = None
        if res_key != cur_res_key:
            cur_chain[1].append((resname, []))
            cur_res_key = res_key
        cur_chain[1][-1][1].append(line)

    if not chains:
        raise ValueError(f"{in_fn}: no standard protein residues found")

    out_chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    serial = 1
    resseq = 1
    lengths = []
    with open(out_fn, "w") as f:
        for ci, (_orig_id, residues) in enumerate(chains):
            cid = out_chain_ids[ci]
            lengths.append(len(residues))
            for resname, atom_lines in residues:
                for al in atom_lines:
                    name = al[12:16]
                    # MSE selenium -> sulphur column, matching Rosetta
                    if name.strip() == "SE" and al[17:20].strip() == "MSE":
                        name = " SD "
                    f.write(
                        "ATOM  %5d %4s %3s %s%4d    %s%s\n"
                        % (serial, name, resname, cid, resseq,
                           al[30:54], "  1.00  0.00")
                    )
                    serial += 1
                resseq += 1
            f.write("TER\n")
        f.write("END\n")

    return len(chains), lengths


def write_complex_pdb(out_fn, aatype, atom_positions, atom_mask, plddt, binderlen):
    """Write the prediction with per-residue pLDDT as B-factor, binder = chain A."""
    atom_types = residue_constants.atom_types
    serial = 1
    with open(out_fn, "w") as f:
        for i in range(len(aatype)):
            cid = "A" if (binderlen is not None and i < binderlen) else "B"
            resnum = i + 1 if cid == "A" else i + 1 - (binderlen or 0)
            resname = residue_constants.restype_1to3[
                residue_constants.restypes[aatype[i]]
            ]
            for j, aname in enumerate(atom_types):
                if atom_mask[i][j] < 0.5:
                    continue
                x, y, z = atom_positions[i][j]
                name = (" %-3s" % aname) if len(aname) < 4 else aname
                f.write(
                    "ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\n"
                    % (serial, name, resname, cid, resnum, x, y, z,
                       plddt[i], aname[0])
                )
                serial += 1
            if binderlen is not None and i + 1 == binderlen:
                f.write("TER\n")
        f.write("TER\nEND\n")


# ================================ runner =====================================

class AF2Runner:
    def __init__(self, args):
        self.args = args
        self.model_name = "model_1_ptm"

        cfg = config.model_config(self.model_name)
        cfg.data.eval.num_ensemble = 1
        cfg.data.common.num_recycle = args.recycle
        cfg.model.num_recycle = args.recycle
        cfg.model.embeddings_and_evoformer.initial_guess = not args.no_initial_guess
        cfg.data.common.max_extra_msa = 5
        cfg.data.eval.max_msa_clusters = 5

        params_dir = os.path.join(DBD, "model_weights")
        params = data.get_model_haiku_params(model_name=self.model_name,
                                             data_dir=params_dir)
        self.runner = model.RunModel(cfg, params)

    def featurize(self, pdb_fn, seq, binderlen, monomer):
        all_pos, all_mask = af2_util.af2_get_atom_positions(pdb_fn)
        initial_guess = af2_util.parse_initial_guess(all_pos)

        if monomer:
            residue_mask = [False for _ in range(len(seq))]
        else:
            residue_mask = [int(i) > binderlen for i in range(len(seq))]

        template_dict = af2_util.generate_template_features(
            seq, all_pos, all_mask, residue_mask
        )

        feature_dict = {
            **pipeline.make_sequence_features(sequence=seq, description="none",
                                              num_res=len(seq)),
            **pipeline.make_msa_features(msas=[[seq]],
                                         deletion_matrices=[[[0] * len(seq)]]),
            **template_dict,
        }

        breaks = [] if monomer else af2_util.check_residue_distances(
            all_pos, all_mask, self.args.max_amide_dist
        )
        feature_dict["residue_index"] = af2_util.insert_truncations(
            feature_dict["residue_index"], breaks
        )
        feature_dict = self.runner.process_features(feature_dict, random_seed=0)
        return feature_dict, initial_guess, all_pos

    def run(self, tag, src_pdb, tmpdir):
        t0 = timer()
        norm_fn = os.path.join(tmpdir, f"{tag}.norm.pdb")
        n_chains, lengths = normalise_pdb(src_pdb, norm_fn)

        seq_chains = af2_util.get_seq_from_pdb(norm_fn)
        seq_chains = [s for s in seq_chains if s]

        if self.args.force_monomer or n_chains == 1:
            monomer, binderlen = True, len(seq_chains[0])
        elif n_chains == 2:
            monomer, binderlen = False, len(seq_chains[0])
        else:
            raise ValueError(
                f"{tag}: {n_chains} chains. Only monomer or binder+target "
                f"(2 chains) are supported; crop the target first."
            )
        seq = "".join(seq_chains)

        feature_dict, initial_guess, init_pos = self.featurize(
            norm_fn, seq, binderlen, monomer
        )

        t_pred = timer()
        result = self.runner.apply(self.runner.params, jax.random.PRNGKey(0),
                                   feature_dict, initial_guess)
        pred_s = timer() - t_pred

        sm = result["structure_module"]
        atom_positions = np.array(sm["final_atom_positions"])
        atom_mask = np.array(sm["final_atom_mask"])
        aatype = np.array(feature_dict["aatype"][0])

        plddt = confidence.compute_plddt(
            np.array(result["predicted_lddt"]["logits"])
        )
        pae = confidence.compute_predicted_aligned_error(
            np.array(result["predicted_aligned_error"]["logits"]),
            np.array(result["predicted_aligned_error"]["breaks"]),
        )["predicted_aligned_error"]

        target_mask = np.zeros(len(seq), dtype=bool)
        target_mask[binderlen:] = True
        rmsds = af2_util.calculate_rmsds(init_pos, atom_positions, target_mask)

        if monomer:
            scores = {
                "plddt_total": float(np.mean(plddt)),
                "plddt_binder": float(np.mean(plddt)),
                "plddt_target": float("nan"),
                "pae_binder": float(np.mean(pae)),
                "pae_target": float("nan"),
                "pae_interaction": float("nan"),
            }
        else:
            scores = {
                "plddt_total": float(np.mean(plddt)),
                "plddt_binder": float(np.mean(plddt[:binderlen])),
                "plddt_target": float(np.mean(plddt[binderlen:])),
                "pae_binder": float(np.mean(pae[:binderlen, :binderlen])),
                "pae_target": float(np.mean(pae[binderlen:, binderlen:])),
                "pae_interaction": float(
                    (np.mean(pae[:binderlen, binderlen:])
                     + np.mean(pae[binderlen:, :binderlen])) / 2
                ),
            }
        scores["binder_aligned_rmsd"] = float(rmsds["binder_aligned_rmsd"])
        scores["target_aligned_rmsd"] = float(rmsds["target_aligned_rmsd"])
        scores["binderlen"] = int(binderlen)
        scores["af2_time"] = round(pred_s, 2)
        scores["total_time"] = round(timer() - t0, 2)

        os.makedirs(self.args.outpdbdir, exist_ok=True)
        write_complex_pdb(
            os.path.join(self.args.outpdbdir, f"{tag}_af2pred.pdb"),
            aatype, atom_positions, atom_mask, plddt,
            None if monomer else binderlen,
        )
        os.remove(norm_fn)
        return scores


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-pdbdir", required=True)
    p.add_argument("-outpdbdir", default="outputs")
    p.add_argument("-scorefilename", default="out.sc")
    p.add_argument("-checkpoint_name", default="check.point")
    p.add_argument("-runlist", default="")
    p.add_argument("-recycle", type=int, default=3)
    p.add_argument("-max_amide_dist", type=float, default=3.0)
    p.add_argument("-no_initial_guess", action="store_true")
    p.add_argument("-force_monomer", action="store_true")
    p.add_argument("-debug", action="store_true")
    args = p.parse_args()

    dev = xla_bridge.get_backend().platform
    print(f"[af2ig] backend={dev} devices={jax.devices()}")
    if dev != "gpu":
        print("[af2ig] WARNING: no GPU, this will be very slow")

    pdbs = sorted(glob.glob(os.path.join(args.pdbdir, "*.pdb")))
    if args.runlist:
        keep = {l.strip() for l in open(args.runlist) if l.strip()}
        pdbs = [p_ for p_ in pdbs
                if os.path.basename(p_).rsplit(".", 1)[0] in keep]

    # --- sort by complex length: this is worth hours ------------------------
    # JAX retraces the whole model for every new input shape. Measured here:
    # 122 s on a new length, 10.5 s on a repeat. Processing in filename order
    # groups designs by BACKBONE, and each backbone has its own length, so
    # every backbone pays the 122 s compile -- the worst possible ordering.
    # Sorting by total residue count means each distinct length compiles once.
    # On 1800 designs spanning ~31 lengths that is ~63 min of compile instead
    # of ~15 h.
    def n_res(p_):
        seen = set()
        for l in open(p_):
            if l[:4] == "ATOM" and l[12:16].strip() == "CA":
                seen.add((l[21], l[22:27]))
        return len(seen)

    sizes = {p_: n_res(p_) for p_ in pdbs}
    pdbs.sort(key=lambda p_: (sizes[p_], p_))
    n_shapes = len(set(sizes.values()))
    print(f"[af2ig] {len(pdbs)} designs, {n_shapes} distinct lengths "
          f"-> {n_shapes} compiles (~{n_shapes * 122 / 60:.0f} min) "
          f"+ ~{(len(pdbs) - n_shapes) * 10.5 / 60:.0f} min inference")

    finished = set()
    if os.path.isfile(args.checkpoint_name):
        finished = {l.strip() for l in open(args.checkpoint_name) if l.strip()}

    os.makedirs(args.outpdbdir, exist_ok=True)
    tmpdir = os.path.join(args.outpdbdir, ".tmp")
    os.makedirs(tmpdir, exist_ok=True)

    runner = AF2Runner(args)
    n_ok = n_fail = 0
    wall0 = timer()

    for src in pdbs:
        tag = os.path.basename(src).rsplit(".", 1)[0]
        if tag in finished:
            print(f"[af2ig] {tag} already done, skipping")
            continue
        try:
            scores = runner.run(tag, src, tmpdir)
            write_header = not os.path.isfile(args.scorefilename)
            af2_util.add2scorefile(tag, args.scorefilename, write_header,
                                   scores, None)
            print(f"[af2ig] {tag}  pae_int={scores['pae_interaction']:.2f}  "
                  f"plddt_binder={scores['plddt_binder']:.1f}  "
                  f"{scores['total_time']}s")
            n_ok += 1
        except KeyboardInterrupt:
            sys.exit("interrupted")
        except Exception as e:
            n_fail += 1
            print(f"[af2ig] {tag} FAILED: {type(e).__name__}: {e}")
            if args.debug:
                raise
        with open(args.checkpoint_name, "a") as f:
            f.write(tag + "\n")

    wall = timer() - wall0
    done = n_ok + n_fail
    print(f"\n[af2ig] {n_ok} ok, {n_fail} failed, {wall:.1f}s total"
          + (f", {wall/done:.1f}s/design" if done else ""))


if __name__ == "__main__":
    main()
