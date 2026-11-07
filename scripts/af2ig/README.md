# AF2-initial-guess without PyRosetta

Upstream `dl_binder_design/af2_initial_guess/predict.py` (commit `cafa3853`)
imports PyRosetta at module scope in both `predict.py` and `af2_util.py`, even
for `-pdbdir` input. We have no PyRosetta licence, so this directory holds a
`-pdbdir`-only reimplementation.

## Why this is a small change, not a rewrite

PyRosetta was doing I/O, not science. The whole feature path is PDB **text**:

```
af2_util.af2_get_atom_positions(pose, tmp_fn):
    pose.dump_pdb(tmp_fn)                  # pose  -> PDB text
    lines = open(tmp_fn).readlines()       # PDB text -> coords
    os.remove(tmp_fn)
```

With `-pdbdir` the input is already a PDB, so `pose_from_pdb` → `dump_pdb` →
re-parse is a round-trip. Feed the same PDB text and the features are identical.
Scoring (`pae_interaction`, `plddt_binder`, RMSDs) is pure numpy on AF2 outputs
and never touched Rosetta at all.

The one thing the round-trip *did* do was normalise the PDB. `normalise_pdb()`
in `predict_pdb.py` reproduces that explicitly: standard residues only, MSE→MET,
no hydrogens, first altloc only, residues renumbered contiguously across chains,
`TER` between chains. That also enforces the unique-residue-index invariant
upstream checked for in `input_check`.

## Files

| File | Role |
|---|---|
| `make_af2_util_nopyro.py` | Generator. Applies 4 asserted edits to upstream `af2_util.py`. Fails loudly if upstream changed. |
| `af2_util_nopyro.py` | **Generated — do not hand-edit.** |
| `predict_pdb.py` | Driver. `-pdbdir` only. |

Regenerate with `python make_af2_util_nopyro.py`.

## Edits applied to af2_util.py

1. Drop the `pyrosetta` / `rosetta` module-scope imports.
2. `af2_get_atom_positions(pose, tmp_fn)` → `af2_get_atom_positions(pdb_fn)`;
   read the file directly. Parsing below the round-trip is byte-for-byte upstream.
3. Fix a latent upstream `NameError`: the selenomethionine branch referenced
   `res.get_resname()`, and `res` is not in scope. It never fired upstream
   because Rosetta converted MSE→MET before this code saw it; reading raw PDBs
   it would. Now tests `atom[1] == 'MSE'`.
4. Excise `insert_Rosetta_chainbreaks` (output-side pose surgery; the driver
   writes chains itself).

## Not supported

* `-silent` / `-outsilent`. Silent files are a Rosetta format and are the only
  part that genuinely required PyRosetta.
* More than 2 chains. Crop the target first.

## Environment

conda env `af2ig`. The version pins are load-bearing and were painful:

| Package | Pin | Why |
|---|---|---|
| `jax` / `jaxlib` | 0.4.20 (+cuda12.cudnn89) | AF2 2.3 code; newer jax drops `jax.core` APIs haiku 0.0.11 uses |
| `dm-haiku` | 0.0.11 | contemporary with jax 0.4.20 |
| `numpy` | 1.26.4 | jax 0.4.20 predates numpy 2.0 (`np.trapz` removed) |
| `ml-dtypes` | 0.3.1 | required by jax 0.4.20 |
| `tensorflow-cpu` | 2.16.1 | AF2's `model/tf` feature pipeline; 2.21 demands ml-dtypes ≥0.5.1 and conflicts with jax |
| `nvidia-*-cu12` | CUDA 12.3 era | pip resolved these to 12.9 by default, and `jaxlib 0.4.20` then failed with "Unable to load cuSOLVER" and silently fell back to CPU |

Biopython is 1.88 and `Bio.Data.SCOPData` is gone (removed in 1.80). It is
shimmed from `Bio.Data.PDBData` at the top of `predict_pdb.py`; it is only
reached through the template-search import chain, which we never execute.
