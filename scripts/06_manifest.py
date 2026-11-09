#!/usr/bin/env python
"""Build results/processed/manifest.parquet (+ CSV mirror) to the CLAUDE.md §5 contract.

Every column in the contract is emitted. Columns with no tool behind them on
this box are written as the literal string `not_run`, per CLAUDE.md §2 — never
null-filled with a plausible-looking number. As of the Stage 00 audit that is:

  rosetta_ddg_g12d / rosetta_ddg_wt   PyRosetta is NOT installed. docs/PIPELINE.md §0
                                      lists it under "Keep (proven in comp2/3/4)", but
                                      comp3 and comp4 actually ran
                                      scripts/af2ig/af2_util_nopyro.py — the explicit
                                      no-PyRosetta variant. It was never in the stack.
  prodigy_g12d                        PRODIGY not installed.
  thermompnn_score                    ThermoMPNN not installed.
  anarci_class                        ANARCI not installed; Track B is gated at H+14
                                      anyway, and Track A designs are not antibodies.
  electrostatic_delta                 APBS not installed. Flagged as a hypothesis metric
                                      in docs/PIPELINE.md, not a gate.
  mpnn_peptide_recovery_gap           Stage 2 prefilter, not run in the pilot.

`peptide_contact_fraction` is deliberately NOT collapsed to one number here.
The contract names a single column, but the measured ceiling on this target
makes the area-based and count-based definitions non-interchangeable (see
results/RUNLOG.md Stage 01). Both are carried, plus the absolute buried area,
and the contract column is set to the area definition so downstream code
reading the contract name still works.

Run:
  conda run -n boltz python scripts/06_manifest.py \
      --boltz pilot/boltz/pilot_three_state.csv \
      --interface results/processed/interface_metrics.csv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

NOT_RUN = "not_run"

COLUMNS = [
    "candidate_id", "sequence", "length", "molecule_class", "track",
    "design_method", "model_version", "params_hash", "seed", "target_state",
    "input_structure", "peptide_contact_fraction", "n_peptide_contacts",
    "iptm_g12d", "iptm_wt", "iptm_empty", "iptm_margin_wt", "iptm_margin_empty",
    "pae_interface_g12d", "rosetta_ddg_g12d", "rosetta_ddg_wt", "prodigy_g12d",
    "mpnn_peptide_recovery_gap", "electrostatic_delta", "thermompnn_score",
    "anarci_class", "qc_status", "reject_reason", "diversity_cluster", "rank",
    "artifact_path",
]

# extra columns beyond the contract, kept because the contract's single
# peptide_contact_fraction is ambiguous on this target
EXTRA = [
    "peptide_contact_fraction_area", "peptide_contact_fraction_count",
    "bsa_peptide", "bsa_total_target", "hotspots_contacted",
    "n_hotspots_contacted", "peptide_positions_contacted",
    "iptm_margin_wt_peptide", "iptm_margin_wt_mhc",
    "iptm_binder_peptide_g12d", "iptm_binder_peptide_wt",
    "mpnn_pool", "sidechains_present",
    "pae_interface_mhc_g12d", "pae_interface_peptide_g12d",
]


def git_rev() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


def params_hash(cfg_path: Path) -> str:
    return hashlib.sha256(cfg_path.read_bytes()).hexdigest()[:12]


def read_csv(path: Path) -> dict[str, dict]:
    if not path or not path.exists():
        return {}
    with open(path) as fh:
        return {r.get("candidate_id") or r.get("design"): r
                for r in csv.DictReader(fh)}


def num(v, default=NOT_RUN):
    if v in (None, "", NOT_RUN, "None"):
        return default
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boltz", type=Path, default=ROOT / "pilot/boltz/pilot_three_state.csv")
    ap.add_argument("--interface", type=Path, default=None,
                    help="interface metrics CSV keyed by design/candidate_id")
    ap.add_argument("--out", type=Path, default=ROOT / "results/processed/manifest")
    ap.add_argument("--track", default="pilot")
    ap.add_argument("--design-method", default="rfdiffusion+solublempnn")
    ap.add_argument("--seed", type=int, default=37)
    args = ap.parse_args()

    cfg_path = ROOT / "configs/challenge.yaml"
    cfg = yaml.safe_load(open(cfg_path))
    ph = params_hash(cfg_path)
    rev = git_rev()

    boltz = read_csv(args.boltz)
    iface = read_csv(args.interface) if args.interface else {}
    if not boltz:
        raise SystemExit(f"no boltz rows in {args.boltz}")

    rows = []
    for cid, b in sorted(boltz.items()):
        # the interface CSV is keyed by backbone stem; a candidate is
        # <backbone>_s<N>, so strip the sample suffix to find its backbone
        # ids are <pool>_<backbone>_s<N>; strip both ends to reach the backbone
        _core = cid.split("_", 1)[1] if cid.split("_", 1)[0] in ("neutral", "basic") else cid
        bb = _core.rsplit("_s", 1)[0]
        i = iface.get(cid) or iface.get(_core) or iface.get(bb) or {}

        seq = b.get("sequence", "")
        length = int(b.get("length") or len(seq) or 0)
        mol = ("protein_minibinder_cat1"
               if cfg["challenge"]["minibinder_range"][0] <= length
               <= cfg["challenge"]["minibinder_range"][1]
               else "single_chain_protein")

        frac_area = i.get("peptide_contact_fraction_area", NOT_RUN)
        frac_count = i.get("peptide_contact_fraction_count", NOT_RUN)

        row = {
            "candidate_id": cid,
            "sequence": seq,
            "length": length,
            "molecule_class": mol,
            "track": args.track,
            "design_method": args.design_method,
            "model_version": f"rfdiffusion:docker@latest|proteinmpnn:soluble_v_48_020|"
                             f"boltz:2.2.1|repo:{rev}",
            "params_hash": ph,
            "seed": args.seed,
            "target_state": "g12d|wt|empty",
            "input_structure": "structures/crop/g12d_crop.pdb (gen) / "
                               "structures/{g12d,wt,empty}.pdb (score)",
            "peptide_contact_fraction": frac_area,       # contract column = area defn
            "n_peptide_contacts": i.get("n_peptide_contacts", NOT_RUN),
            "iptm_g12d": num(b.get("g12d_iptm_binder_mhc")),
            "iptm_wt": num(b.get("wt_iptm_binder_mhc")),
            "iptm_empty": num(b.get("empty_iptm_binder_mhc")),
            "iptm_margin_wt": num(b.get("iptm_margin_wt_peptide")),
            "iptm_margin_empty": num(b.get("iptm_margin_empty")),
            "pae_interface_g12d": num(b.get("g12d_pae_interface_peptide"),
                                 num(b.get("g12d_pae_interface_mhc"))),
            "rosetta_ddg_g12d": NOT_RUN,   # PyRosetta absent - see docstring
            "rosetta_ddg_wt": NOT_RUN,
            "prodigy_g12d": NOT_RUN,
            "mpnn_peptide_recovery_gap": NOT_RUN,
            "electrostatic_delta": NOT_RUN,
            "thermompnn_score": NOT_RUN,
            "anarci_class": NOT_RUN,
            "qc_status": b.get("g12d_status", "unknown"),
            "reject_reason": i.get("reject_reason", ""),
            "diversity_cluster": NOT_RUN,
            "rank": NOT_RUN,
            "artifact_path": f"pilot/boltz/{args.track}/g12d/pred",
            # extras
            "peptide_contact_fraction_area": frac_area,
            "peptide_contact_fraction_count": frac_count,
            "bsa_peptide": i.get("bsa_peptide", NOT_RUN),
            "bsa_total_target": i.get("bsa_total_target", NOT_RUN),
            "hotspots_contacted": i.get("hotspots_contacted", NOT_RUN),
            "n_hotspots_contacted": i.get("n_hotspots_contacted", NOT_RUN),
            "peptide_positions_contacted": i.get("peptide_positions_contacted", NOT_RUN),
            "iptm_margin_wt_peptide": num(b.get("iptm_margin_wt_peptide")),
            "iptm_margin_wt_mhc": num(b.get("iptm_margin_wt_mhc")),
            "iptm_binder_peptide_g12d": num(b.get("g12d_iptm_binder_peptide")),
            "iptm_binder_peptide_wt": num(b.get("wt_iptm_binder_peptide")),
            # from the id prefix, not the file path: one CSV holds both pools
            "mpnn_pool": (cid.split("_", 1)[0]
                          if cid.split("_", 1)[0] in ("neutral", "basic") else NOT_RUN),
            "sidechains_present": i.get("sidechains_present", NOT_RUN),
            "pae_interface_mhc_g12d": num(b.get("g12d_pae_interface_mhc")),
            "pae_interface_peptide_g12d": num(b.get("g12d_pae_interface_peptide")),
        }
        rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = COLUMNS + EXTRA
    csv_path = args.out.with_suffix(".csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {rel(csv_path, ROOT)}  ({len(rows)} rows, {len(cols)} cols)")

    try:
        import pandas as pd
        pq = args.out.with_suffix(".parquet")
        pd.DataFrame(rows, columns=cols).to_parquet(pq, index=False)
        print(f"wrote {rel(pq, ROOT)}")
    except Exception as e:
        print(f"parquet mirror skipped: {type(e).__name__}: {e}", file=sys.stderr)
        print("  (CSV is authoritative; install pyarrow/fastparquet for the parquet mirror)",
              file=sys.stderr)

    missing = [c for c in COLUMNS if all(r[c] == NOT_RUN for r in rows)]
    print(f"\ncontract columns that are not_run for every row ({len(missing)}): {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
