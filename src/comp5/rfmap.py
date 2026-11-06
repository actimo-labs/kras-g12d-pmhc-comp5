"""RFdiffusion output <-> input residue mapping.

RFdiffusion renumbers and re-chains its output. Measured on this target, input
`A1-180/0 C1-10` returns as chain A = designed binder (1..N) and chain B =
heavy chain AND peptide merged and renumbered N+1..N+190. The peptide offset
therefore MOVES WITH BINDER LENGTH, so it can never be hardcoded -- that is the
comp3 round-2 bug.

Every consumer resolves positions through `trb_map` instead, which reads the
explicit map RFdiffusion writes into the sibling .trb.
"""
from __future__ import annotations

from pathlib import Path

# ProteinMPNN's fixed alphabet order, used for bias_by_res vectors.
MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"


def trb_map(trb: Path) -> dict[tuple[str, int], tuple[str, int]]:
    """output (chain, resnum) -> original (chain, resnum).

    Raises rather than guessing if the map is absent or malformed.
    """
    import numpy as np

    d = np.load(trb, allow_pickle=True)
    for ref_k, hal_k in (("complex_con_ref_pdb_idx", "complex_con_hal_pdb_idx"),
                         ("con_ref_pdb_idx", "con_hal_pdb_idx")):
        if ref_k in d and len(d[ref_k]):
            ref, hal = d[ref_k], d[hal_k]
            if len(ref) != len(hal):
                raise ValueError(f"{trb.name}: ref/hal length mismatch "
                                 f"{len(ref)} vs {len(hal)}")
            return {(str(h[0]), int(h[1])): (str(r[0]), int(r[1]))
                    for r, h in zip(ref, hal)}
    raise ValueError(f"{trb.name}: no usable residue map; refusing to guess an offset")


def output_residue_for(trb: Path, orig_chain: str, orig_num: int) -> tuple[str, int]:
    """Where did original (chain, num) end up in the output file?"""
    for out_key, orig in trb_map(trb).items():
        if orig == (orig_chain, orig_num):
            return out_key
    raise KeyError(f"{trb.name}: original {orig_chain}{orig_num} not present in map")


def binder_chain(trb: Path, chains_present: set[str]) -> str:
    """The binder is the only chain with no entries in the map (it is de novo)."""
    mapped = {c for c, _ in trb_map(trb)}
    free = sorted(chains_present - mapped)
    if len(free) != 1:
        raise ValueError(f"{trb.name}: expected exactly 1 unmapped binder chain, got {free}")
    return free[0]
