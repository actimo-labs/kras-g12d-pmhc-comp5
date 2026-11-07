#!/usr/bin/env python3
"""
Generate af2_util_nopyro.py from the upstream dl_binder_design af2_util.py.

Why this exists
---------------
AF2-initial-guess hard-imports PyRosetta at module scope in both predict.py and
af2_util.py. We have no PyRosetta licence (see runlog/00-environment-audit.md,
finding F3), so the dependency has to go.

The dependency turns out to be almost entirely I/O boilerplate. Only two
functions in af2_util.py touch a pose:

  af2_get_atom_positions(pose, tmp_fn)
      dumps the pose to a PDB, re-reads the text, parses coordinates out of it,
      deletes the temp file. When the input is ALREADY a PDB (-pdbdir mode) the
      pose round-trip is a no-op, so we read the file directly. The parsing code
      below the round-trip is untouched, so features are identical.

  insert_Rosetta_chainbreaks(pose, binderlen)
      output-side only: splits the predicted structure into binder/target
      chains. predict_pdb.py writes chains itself, so this is dropped.

Everything else in af2_util.py is pure numpy/jax and is copied verbatim.

This is a GENERATOR rather than a hand-edited copy so the diff against upstream
stays auditable for the methods report. Each edit asserts that it matched; if
upstream ever changes, this fails loudly instead of silently producing a file
that differs from what we think it does.

Usage:  python make_af2_util_nopyro.py
"""

import os
import sys

UPSTREAM = os.path.expanduser(
    "~/tnf-track3/dl_binder_design/af2_initial_guess/af2_util.py"
)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "af2_util_nopyro.py")

HEADER = '''#!/usr/bin/env python3
# ===================================================================
# GENERATED FILE -- do not edit by hand.
# Produced by scripts/af2ig/make_af2_util_nopyro.py from upstream
# dl_binder_design/af2_initial_guess/af2_util.py (commit cafa3853).
# Re-run the generator to regenerate. See that script for the rationale.
# ===================================================================
'''

# (description, exact text to find, replacement)
EDITS = [
    (
        "drop the PyRosetta module-scope import",
        """from pyrosetta import *
from rosetta import *
# I don't think we need to call init() here but I'm not sure
""",
        """# PyRosetta imports removed -- see make_af2_util_nopyro.py.
""",
    ),
    (
        "af2_get_atom_positions: take a PDB path instead of a pose",
        """def af2_get_atom_positions(pose, tmp_fn) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Given a pose, return the AF2 atom positions array and atom mask array for the protein.
    '''

    # write pose to pdb file
    pose.dump_pdb(tmp_fn)

    with open(tmp_fn, 'r') as pdb_file:
        lines = pdb_file.readlines()

    # Delete the temporary file
    os.remove(tmp_fn)
""",
        """def af2_get_atom_positions(pdb_fn) -> Tuple[np.ndarray, np.ndarray]:
    '''
    Given a PDB file, return the AF2 atom positions array and atom mask array.

    Upstream signature was (pose, tmp_fn): it dumped the pose to tmp_fn and read
    the text straight back. We are handed a normalised PDB already, so we read
    it directly. The parsing below is byte-for-byte upstream.
    '''

    with open(pdb_fn, 'r') as pdb_file:
        lines = pdb_file.readlines()
""",
    ),
    (
        "fix upstream NameError on selenomethionine ('res' is not in scope)",
        """            elif atom_name.upper() == 'SE' and res.get_resname() == 'MSE':""",
        """            elif atom_name.upper() == 'SE' and atom[1] == 'MSE':""",
    ),
]


def main():
    if not os.path.isfile(UPSTREAM):
        sys.exit(f"upstream not found: {UPSTREAM}")

    src = open(UPSTREAM).read()

    for desc, find, repl in EDITS:
        n = src.count(find)
        if n != 1:
            sys.exit(
                f"FAILED edit ({desc}): expected exactly 1 match, found {n}.\n"
                f"Upstream af2_util.py has changed; review before regenerating."
            )
        src = src.replace(find, repl)
        print(f"  ok  {desc}")

    # Excise insert_Rosetta_chainbreaks entirely (output-side pose surgery).
    start = src.index("def insert_Rosetta_chainbreaks(")
    end = src.index("def check_residue_distances(")
    src = src[:start] + src[end:]
    print("  ok  excise insert_Rosetta_chainbreaks")

    # Strip the now-unused shebang from upstream, prepend our header.
    src = src.replace("#!/usr/bin/env python3\n", "", 1)
    src = HEADER + src

    # Belt and braces: nothing rosetta-shaped may survive.
    for bad in ("pyrosetta", "from rosetta", "pose.", "core.pose"):
        if bad in src:
            sys.exit(f"FAILED: '{bad}' still present in generated file")

    open(OUT, "w").write(src)
    print(f"wrote {OUT} ({len(src.splitlines())} lines)")


if __name__ == "__main__":
    main()
