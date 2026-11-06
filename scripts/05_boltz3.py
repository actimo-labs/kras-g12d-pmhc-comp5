#!/usr/bin/env python
"""Three-state Boltz-2 co-folding: g12d / wt / empty.

Ported from tfr1-comp4/scripts/04_boltz.py. Same `msa: empty` single-sequence
mode that comp3 and comp4 both used -- the self-hosted MSA server at
10.148.0.19:8080 is DOWN (Stage 00 audit) and was never on the proven path.

WHICH NUMBER WE READ, AND WHY NOT GLOBAL ipTM
---------------------------------------------
This target is FOUR chains: binder + HLA heavy + b2m + peptide. Global ipTM
averages every pair, including heavy<->b2m, which is a real, tight, native
interface that scores high no matter what the binder does. Reading global ipTM
here measures the MHC assembling with itself. We read `pair_chains_iptm` for
binder<->target-chain pairs only.

Two pair scores, because they answer different questions:

  iptm_binder_mhc      binder <-> HLA heavy chain. Present in ALL THREE states,
                       so this is the only score comparable across the empty
                       differential (the `empty` state has no peptide chain to
                       pair against).
  iptm_binder_peptide  binder <-> peptide. Exists for g12d and wt only. The
                       most mutation-sensitive quantity, so it is the primary
                       signal for the WT differential.

Margins follow from that:
  iptm_margin_wt    = iptm_binder_peptide(g12d) - iptm_binder_peptide(wt)
                      and the same on iptm_binder_mhc as a cross-check
  iptm_margin_empty = iptm_binder_mhc(g12d) - iptm_binder_mhc(empty)

A framework binder sitting on the a1/a2 helices should show margin_empty ~ 0,
because removing the peptide does not disturb its epitope. That is the whole
point of the gate.

TEMPLATES
---------
Each state is templated from its own prepared structure
(`structures/{state}.cif`) so the groove geometry is the crystallographic one
and is IDENTICAL across states apart from the intended difference. Without
that the three states differ by prediction noise as much as by construction,
and the differential measures nothing. Mechanics of the CIF (entities +
label_seq_id + author chain naming) are in scripts/01_prep_targets.py; without
them boltz's parse_polymer dies with a bare IndexError.

Boltz-2's affinity head is NOT used. docs/PIPELINE.md and CLAUDE.md both say it
is weak for protein-protein; nothing here is reported as a K_D.

Runs under the `boltz` conda env:
  conda run -n boltz python scripts/05_boltz3.py --fasta-dir pilot/mpnn/neutral/seqs --limit 16
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import yaml as pyyaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from comp5.paths import rel  # noqa: E402

BINDER_ID = "X"
HEAVY_ID, B2M_ID, PEP_ID = "A", "B", "C"
STATES = ("g12d", "wt", "empty")

# Composition pre-filter: these are automatic developability rejects later, and
# folding them is wasted GPU. Same rule as comp4 stage 04.
AGG_RUN = re.compile(r"([AILMFWVY])\1{4,}")
SEQUON = re.compile(r"N[^P][ST]")


def read_fastas(dirs: list[Path], per_backbone: int | None) -> dict[str, str]:
    """ProteinMPNN fasta: record 0 is the input backbone, 1.. are designs."""
    out: dict[str, str] = {}
    for d in dirs:
        for fa in sorted(d.glob("*.fa")):
            recs, name, buf = [], None, []
            for line in fa.read_text().splitlines():
                if line.startswith(">"):
                    if name is not None:
                        recs.append((name, "".join(buf)))
                    name, buf = line[1:], []
                elif line.strip():
                    buf.append(line.strip())
            if name is not None:
                recs.append((name, "".join(buf)))
            designs = recs[1:]
            if per_backbone:
                designs = designs[:per_backbone]
            for hdr, seq in designs:
                # scripts/04b_select.py writes one design per file already named
                # <backbone>_s<N>.fa. Appending another _s<N> here would produce
                # ids like p_0_s1_s1, which then break the backbone join in
                # scripts/06_manifest.py (rsplit("_s",1) would yield p_0_s1, not
                # p_0). If the stem already carries a sample index, trust it.
                if re.search(r"_s\d+$", fa.stem):
                    key = fa.stem
                else:
                    m = re.search(r"sample=(\d+)", hdr)
                    key = f"{fa.stem}_s{m.group(1) if m else len(out)}"
                # chain A only; MPNN writes designed chains slash-joined
                out[key] = seq.split("/")[0]
    return out


def composition_reject(seq: str) -> str:
    if "C" in seq:
        return "free_cys"
    if AGG_RUN.search(seq):
        return "hydrophobic_run>=5"
    if SEQUON.search(seq):
        return "n_glycosylation_sequon"
    return ""


def write_yaml(path: Path, binder: str, target: dict[str, str], template: Path) -> list[str]:
    """Returns the chain order, which is what pair_chains_iptm is indexed by."""
    order = [BINDER_ID] + [c for c in (HEAVY_ID, B2M_ID, PEP_ID) if c in target]
    seqs = {BINDER_ID: binder, **target}
    lines = ["version: 1", "sequences:"]
    for cid in order:
        lines += ["  - protein:", f"      id: {cid}",
                  f"      sequence: {seqs[cid]}", "      msa: empty"]
    tchains = [c for c in order if c != BINDER_ID]
    lines += ["templates:", f"  - cif: {template}",
              f"    chain_id: [{', '.join(tchains)}]",
              f"    template_id: [{', '.join(tchains)}]"]
    path.write_text("\n".join(lines) + "\n")
    return order


def parse_confidence(pred_root: Path, name: str, order: list[str],
                     average_samples: bool = False) -> dict | None:
    """Read binder<->target pair ipTM.

    With `average_samples`, average across ALL diffusion samples (model_0,
    model_1, ...) instead of taking model_0 alone. This exists because the
    pilot measured single-sample `iptm_margin_wt` to have signal-to-noise of
    0.84 -- the seed-to-seed variance of one design exceeded the spread across
    all 32 designs, so the single-sample ranking was measuring the RNG. If
    averaging n samples cuts the noise by ~sqrt(n), it is far cheaper than n
    separate runs because the model loads once.
    """
    hits = sorted(pred_root.rglob(f"confidence_{name}_model_*.json"))
    if not hits:
        return None
    if not average_samples:
        hits = [h for h in hits if h.name.endswith("_model_0.json")] or hits[:1]

    datas = [json.loads(h.read_text()) for h in hits]

    def avg(key, default=None):
        vals = [d.get(key) for d in datas]
        vals = [v for v in vals if isinstance(v, (int, float))]
        return sum(vals) / len(vals) if vals else default

    data = datas[0]
    pair_list = [d.get("pair_chains_iptm") or {} for d in datas]
    pair = pair_list[0]

    def pair_val(i: int, j: int):
        got = []
        for pc in pair_list:                      # average over samples too
            for a, b in ((i, j), (j, i)):
                try:
                    got.append(float(pc[str(a)][str(b)]))
                except (KeyError, TypeError, ValueError):
                    continue
        return sum(got) / len(got) if got else None

    bi = order.index(BINDER_ID)
    out = {
        "n_samples_averaged": len(datas),
        "complex_plddt": avg("complex_plddt"),
        "complex_iplddt": avg("complex_iplddt"),
        "global_iptm": avg("iptm"),        # recorded, NEVER used to rank
        "ptm": avg("ptm"),
    }
    for cid, label in ((HEAVY_ID, "mhc"), (B2M_ID, "b2m"), (PEP_ID, "peptide")):
        out[f"iptm_binder_{label}"] = (
            pair_val(bi, order.index(cid)) if cid in order else None)
    return out


def parse_pae(pred_root: Path, name: str, order: list[str],
              lengths: dict[str, int]) -> dict:
    """Mean interface PAE for binder<->each target chain, from the pae npz.

    Chain blocks are laid out in YAML order, so offsets come from `order` and
    the sequence lengths actually written into the YAML -- not from the
    structure file, which may model fewer residues than the query.
    """
    import numpy as np

    hits = sorted(pred_root.rglob(f"pae_{name}_model_0.npz"))
    if not hits:
        return {}
    try:
        arr = np.load(hits[0])
        pae = arr[arr.files[0]]
    except Exception:
        return {}
    if pae.ndim != 2:
        return {}

    spans, off = {}, 0
    for cid in order:
        n = lengths[cid]
        spans[cid] = (off, off + n)
        off += n
    if off != pae.shape[0]:
        # layout does not match the query; refuse rather than mis-slice
        return {"pae_interface_note": f"shape {pae.shape[0]} != expected {off}"}

    bs, be = spans[BINDER_ID]
    out = {}
    for cid, label in ((HEAVY_ID, "mhc"), (B2M_ID, "b2m"), (PEP_ID, "peptide")):
        if cid not in spans:
            continue
        ts, te = spans[cid]
        block = np.concatenate([pae[bs:be, ts:te].ravel(), pae[ts:te, bs:be].ravel()])
        out[f"pae_interface_{label}"] = round(float(block.mean()), 3)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", type=Path, nargs="+", required=True)
    ap.add_argument("--out", type=Path, default=ROOT / "pilot/boltz")
    ap.add_argument("--tag", default="pilot")
    ap.add_argument("--states", default=",".join(STATES))
    ap.add_argument("--per-backbone", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--recycling-steps", type=int, default=3)
    ap.add_argument("--diffusion-samples", type=int, default=1)
    ap.add_argument("--devices", type=int, default=1)
    ap.add_argument("--no-composition-filter", action="store_true")
    ap.add_argument("--seed", type=int, default=None,
                    help="Boltz RNG seed. Set it to measure run-to-run variance: "
                         "without a noise estimate a margin spread cannot be "
                         "distinguished from sampling jitter.")
    ap.add_argument("--average-samples", action="store_true",
                    help="average over all diffusion samples, not just model_0")
    ap.add_argument("--only", default=None,
                    help="comma-separated candidate_ids to restrict to")
    args = ap.parse_args()

    cfg = pyyaml.safe_load(open(ROOT / "configs/challenge.yaml"))
    heavy = cfg["target"]["heavy_chain"]["seq"]
    b2m = cfg["target"]["b2m"]["seq"]
    states = [s.strip() for s in args.states.split(",") if s.strip()]

    state_target = {
        "g12d": {HEAVY_ID: heavy, B2M_ID: b2m, PEP_ID: cfg["states"]["g12d"]["peptide"]},
        "wt": {HEAVY_ID: heavy, B2M_ID: b2m, PEP_ID: cfg["states"]["wt"]["peptide"]},
        "empty": {HEAVY_ID: heavy, B2M_ID: b2m},
    }

    designs = read_fastas(list(args.fasta_dir), args.per_backbone)
    if args.only:
        want = {x.strip() for x in args.only.split(",") if x.strip()}
        missing = want - set(designs)
        if missing:
            raise SystemExit(f"--only ids not found: {sorted(missing)}")
        designs = {k: v for k, v in designs.items() if k in want}
    if args.limit:
        designs = dict(sorted(designs.items())[:args.limit])
    if not designs:
        raise SystemExit(f"no sequences found in {args.fasta_dir}")

    rejected = {}
    if not args.no_composition_filter:
        for k, v in list(designs.items()):
            why = composition_reject(v)
            if why:
                rejected[k] = why
                designs.pop(k)
    print(f"{len(designs)} designs to fold, {len(rejected)} composition-rejected")
    if rejected:
        rp = ROOT / "results/rejected/composition_prefilter.csv"
        rp.parent.mkdir(parents=True, exist_ok=True)
        new = not rp.exists()
        with open(rp, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["candidate_id", "stage", "reject_reason"])
            for k, v in rejected.items():
                w.writerow([k, "05_boltz3_composition_prefilter", v])
        print(f"  appended to {rel(rp, ROOT)}")

    args.out.mkdir(parents=True, exist_ok=True)
    orders: dict[tuple[str, str], list[str]] = {}
    timings = {}

    for state in states:
        tmpl = ROOT / f"structures/{state}.cif"
        if not tmpl.exists():
            raise SystemExit(f"missing template {tmpl}; run scripts/01_prep_targets.py")
        ydir = args.out / args.tag / state / "yaml"
        ydir.mkdir(parents=True, exist_ok=True)
        for name, seq in designs.items():
            orders[(state, name)] = write_yaml(
                ydir / f"{name}.yaml", seq, state_target[state], tmpl)

        pdir = args.out / args.tag / state / "pred"
        pdir.mkdir(parents=True, exist_ok=True)
        cmd = ["boltz", "predict", str(ydir),
               "--out_dir", str(pdir),
               "--cache", str(Path.home() / ".boltz"),
               "--recycling_steps", str(args.recycling_steps),
               "--diffusion_samples", str(args.diffusion_samples),
               "--devices", str(args.devices),
               "--output_format", "pdb",
               "--override"]
        if args.seed is not None:
            cmd += ["--seed", str(args.seed)]
        print(f"\n=== state {state}: {len(designs)} predictions ===")
        print("+", " ".join(cmd), flush=True)
        t0 = time.time()
        rc = subprocess.run(cmd).returncode
        dt = time.time() - t0
        timings[state] = dt
        print(f"state {state}: rc={rc} wall={dt:.0f}s "
              f"({dt/max(len(designs),1):.1f}s/prediction)")
        if rc != 0:
            print(f"WARNING: boltz returned {rc} for state {state}", file=sys.stderr)

    # ---- collect ----------------------------------------------------------
    rows = []
    for name in designs:
        row = {"candidate_id": name, "sequence": designs[name],
               "length": len(designs[name])}
        per = {}
        for state in states:
            pdir = args.out / args.tag / state / "pred"
            order = orders[(state, name)]
            c = parse_confidence(pdir, name, order, args.average_samples)
            if c is None:
                row[f"{state}_status"] = "failed"
                continue
            lengths = {BINDER_ID: len(designs[name]),
                       **{k: len(v) for k, v in state_target[state].items()}}
            c.update(parse_pae(pdir, name, order, lengths))
            row[f"{state}_status"] = "ok"
            per[state] = c
            for k, v in c.items():
                row[f"{state}_{k}"] = v

        if "g12d" in per and "wt" in per:
            a, b = per["g12d"].get("iptm_binder_peptide"), per["wt"].get("iptm_binder_peptide")
            row["iptm_margin_wt_peptide"] = round(a - b, 4) if None not in (a, b) else "not_run"
            a, b = per["g12d"].get("iptm_binder_mhc"), per["wt"].get("iptm_binder_mhc")
            row["iptm_margin_wt_mhc"] = round(a - b, 4) if None not in (a, b) else "not_run"
        if "g12d" in per and "empty" in per:
            a, b = per["g12d"].get("iptm_binder_mhc"), per["empty"].get("iptm_binder_mhc")
            row["iptm_margin_empty"] = round(a - b, 4) if None not in (a, b) else "not_run"
        rows.append(row)

    outcsv = args.out / f"{args.tag}_three_state.csv"
    keys = sorted({k for r in rows for k in r})
    with open(outcsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    try:
        shown = rel(outcsv.resolve(), ROOT)
    except ValueError:
        shown = outcsv
    print(f"\nwrote {shown}  ({len(rows)} rows)")
    tot = sum(timings.values())
    n = len(designs) * len(states)
    print(f"throughput: {tot:.0f}s for {n} predictions = {tot/max(n,1):.1f}s each")
    for s, d in timings.items():
        print(f"  {s}: {d:.0f}s")
    json.dump({"timings_sec": timings, "n_designs": len(designs), "states": states},
              open(args.out / f"{args.tag}_timing.json", "w"), indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
