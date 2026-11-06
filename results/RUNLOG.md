# RUNLOG — comp5 / KRAS G12D pMHC

Append-only. One row per stage: stage, start, end, GPU-h, inputs, outputs, exit status.
Times are Australia/Melbourne. Kickoff = Sat 7 Nov 2026 00:00 (H+0).

---

## Stage 00 — environment audit

| field | value |
|---|---|
| start | 2026-11-07 00:36 (H+0.6) |
| end | 2026-11-07 00:48 |
| GPU-h | 0 |
| inputs | `CLAUDE.md §3` checklist |
| outputs | this section |
| status | **complete, 2 blockers found** |

Measured, not assumed:

- **GPU** — 1× A100-SXM4-80GB, driver 550.90.07, CUDA 12.4, 1 MiB used, 0 processes, no other tenant. Free.
- **CPU/RAM** — 12 cores, 167 GB RAM, no swap.
- **Disk — BLOCKER.** `/dev/sda1` 2.0 T, **1.9 T used, 15 G available, 100%**. Cause is `/home/jaberjaber` = **1.6 T, another user's data**. Our own footprint is 63 G. Not ours to delete.
- **MSA server — DOWN.** `10.148.0.19:8080` unreachable: `/health`, `/`, `/status` all HTTP 000, 100% packet loss on ping. We are on `10.148.0.21/32`. Host is gone, not just the service.
  - *Not actually on the critical path*: comp3 and comp4 both ran Boltz-2 with `msa: empty` (grep `tfr1-comp4/scripts/04_boltz.py:173`, `glp1r-comp3/scripts/13_boltz.py:48`). The proven pattern never used this server.
- **Internet** — RCSB 200, PyPI 200. ColabFold API reachable but rejected a malformed probe (400, not a timeout).
- **Docker** — `rosettacommons/rfdiffusion:latest` (12.2 GB) present. Image ships **9** checkpoints incl. `Complex_base_ckpt.pt`, `ActiveSite_ckpt.pt`, `Complex_Fold_base_ckpt.pt`.
- **conda envs** — `base`, `af2ig`, `boltz`, `mpnn`, `viz`.
  - `boltz` 2.2.1 + Bio, torch, numpy, pandas, scipy, sklearn, gemmi
  - `mpnn` torch 2.5.1+cu121, CUDA available, Bio
  - `af2ig` jax + Bio
- **Weights on disk** — Boltz-2 `boltz2_conf.ckpt` (2.29 GB) + `boltz2_aff.ckpt` (2.06 GB) + `mols/` (45,227 entries). ProteinMPNN vanilla/soluble/ca `v_48_002/010/020/030`. RFdiffusion 3 local ckpts (superseded by the image's 9).

### Missing — contradicts `docs/PIPELINE.md §0` "Keep (proven, zero setup risk)"

| component | PIPELINE.md claim | reality |
|---|---|---|
| PyRosetta FastRelax + InterfaceAnalyzer | "proven in comp2/3/4" | **never installed.** comp3/comp4 ran `scripts/af2ig/af2_util_nopyro.py` — explicitly the no-PyRosetta variant. Interfaces were scored with Bio.PDB SASA. |
| PRODIGY | "proven" | not installed |
| ThermoMPNN | manifest column | not installed |
| ANARCI | "proven" | not installed (Track B only) |
| Chai-1 | "add" | not installed |
| FreeBindCraft | "add" | not installed |
| OpenDDE | "add" | not installed |

Consequence: manifest columns `rosetta_ddg_g12d`, `rosetta_ddg_wt`, `prodigy_g12d`,
`thermompnn_score` have no tool behind them. Per CLAUDE.md §2 they are `not_run`
until either a tool is installed or a documented substitute is agreed.

---

## Stage 00b — config corrections

| field | value |
|---|---|
| start | 2026-11-07 00:48 | | end | 2026-11-07 00:52 |
| status | **complete** |

- `scripts/verify_targets.py`: `CFG` repointed from `/mnt/user-data/uploads/...` to
  `ROOT/configs/challenge.yaml`; FASTA output repointed from a dead scratchpad path to
  `configs/target_states.fasta`. **Exit 0, ALL SEQUENCE CHECKS PASSED.**
- `hotspots.tracked_not_targeted`: `[6]` → `[6, 9]`. p9 was unclassified; Gly9 is one of
  the four critical positions in the Poole & Karuppiah 2022 alanine scan (D6, G4, A5, G9),
  so contact accounting was silently ignoring a position the reference TCR depends on.
  Tracked, **not** promoted to a primary hotspot — see Stage 01 SASA.
- Added `hotspots.low_side_chain_positions: [4,5,7,9]` + note: p4/p7/p9 are Gly and p5 is
  Ala, so this interface is carried by shape and electrostatics, not side-chain packing.
  Weak side-chain contact at those positions must not be scored as design failure.
- Fixed a label bug in `verify_targets.py` that printed p9 as "buried in bulge_down",
  reusing p6's text. p9's burial was never established. It now prints the recorded rationale.

---

## Stage 01 — target prep

| field | value |
|---|---|
| start | 2026-11-07 00:52 | | end | 2026-11-07 01:05 |
| GPU-h | 0 (CPU only) |
| inputs | `structures/raw/7OW6.cif` (RCSB, 2.64 Å) |
| outputs | `structures/{g12d,wt,empty}.{pdb,cif}`, `structures/crop/`, `results/processed/target_prep.json` |
| status | **complete** |

7OW6 fetched and confirmed: *"Crystal structure of a TCR in complex with HLA-A\*11:01 bound
to KRAS G12D peptide (VVVGADGVGK)"*, 2.64 Å. Chains A (heavy, 276 modelled), B (β2m, 99 —
exact match to config), C (peptide `VVVGADGVGK` — exact match), D/E (TCR, removed).

States built, waters + SO4 stripped:
- `g12d` A+B+C, 385 residues
- `wt` D6→G by side-chain atom deletion, 385 residues — verified to read `VVVGAGGVGK`
- `empty` A+B, 375 residues

**Deviations from PIPELINE.md Stage 0, all forced by the audit** (no PyRosetta/OpenMM/PDBFixer):
1. `wt` **not relaxed**. Asp→Gly deletes atoms only (CB/CG/OD1/OD2), so it cannot create a
   clash — there is nothing to relieve. Leaves an unfilled cavity, which is what the real WT
   complex also has to cope with. Recorded as an approximation.
2. `empty` **not relaxed**, same reason. PIPELINE.md already flags peptide-free HLA-A\*11:01
   as physically unstable, so this is a rigid model of an ill-defined target. The `empty`
   differential is directional, not quantitative.
3. **`bulge_up` secondary state: `not_run`.** Needs a modelling tool we do not have. Dropped
   under its own >1 h rule.

### Hotspot assertions (the comp3 round-2 bug)

```
[assert] chain C numbering 1..10 == p1..p10, seq VVVGADGVGK OK
[assert] hotspots [4, 5, 7, 8] resolve to GAGV OK
[assert] wt peptide == VVVGAGGVGK == config wt OK
```

### Measured peptide SASA, TCR-removed state — what a binder actually sees

| pos | aa | n_SC | whole Å² | rel% | side-chain Å² | SC rel% | class |
|---|---|---|---|---|---|---|---|
| 1 | V | 3 | 1.2 | 0.6% | 1.2 | 1.1% | anchor |
| 2 | V | 3 | 2.1 | 1.5% | 0.0 | 0.0% | **unclassified** |
| 3 | V | 3 | 5.9 | 4.6% | 4.8 | 4.9% | **unclassified** |
| 4 | G | 0 | 37.0 | 58.3% | 0.0 | n/a | primary |
| 5 | A | 1 | 26.9 | 27.2% | 9.7 | 16.0% | primary |
| 6 | D | 4 | 21.8 | 18.5% | 8.8 | 10.4% | tracked |
| 7 | G | 0 | 75.5 | 83.8% | 0.0 | n/a | primary |
| 8 | V | 3 | 46.6 | 38.1% | 37.4 | 37.8% | primary |
| 9 | G | 0 | 15.7 | 28.0% | 0.0 | n/a | tracked |
| 10 | K | 5 | 6.7 | 2.6% | 2.4 | 1.6% | anchor |

Findings:
- **The buried-Asp guardrail is CONFIRMED, not assumed.** D6 side chain is 10.4% exposed
  (8.8 of 85.2 Å²), bidentate salt bridge to HLA **Arg114** (OD1–NH1 3.04 Å, OD2–NH2 2.95 Å)
  and H-bond to **Gln70** (OD2–NE2 3.08 Å). A binder cannot read this side chain directly.
  The earlier whole-residue figure of 21.8 Å² is mostly backbone and must not be read as
  side-chain accessibility.
- Primary hotspots [4,5,7,8] are the four most exposed non-anchor positions, and p8 is the
  only one offering real side-chain surface (37.4 Å²). Config is correct.
- p9 measures 15.7 Å² / 28% — genuinely marginal. Correctly tracked, not targeted.
- p2 and p3 are buried (1.5%, 4.6%) and should be reclassified from `unclassified` to
  anchor/buried. **Not edited — flagged for Harjas**, since the brief only authorised the p9 edit.

### The primary hard gate is set above its own physical ceiling

Total peptide surface available to a binder: **239.4 Å²** (223.5 Å² excluding anchors).

| total interface buried | max possible `peptide_contact_fraction_area` |
|---|---|
| 700 Å² | 0.342 |
| 900 Å² | 0.266 |
| 1100 Å² | 0.218 |
| **1200 Å² (Baker mart1-3 reference)** | **0.199** |

`configs/challenge.yaml` sets `peptide_contact_fraction_min = 0.25`. At the Baker reference
interface size the ceiling is **0.199 — below the gate**. A design reproducing the Baker
geometry would be rejected by our own filter. The threshold is marked `provisional` and must
be recalibrated; it cannot be used as written. **Needs a decision — see report.**

### Crop for generation

Domain crop, **not** the configured 15 Å distance shell: the shell selects 158 chain-A
residues in **8 discontiguous segments** plus 8 stray β2m residues, and an 8-segment
contigmap risks inventing peptide bonds across boundaries. The α1/α2 platform (A1–180) is
**one contiguous segment**, only 22 residues larger, and strictly contains the shell
(15 Å shell spans A3–A179 — asserted in code).

```
chain A: 180 residues, 1 segment [(1, 180)]
chain C:  10 residues, 1 segment [(1, 10)]
[assert] cropped chain C == VVVGADGVGK, numbering preserved -> hotspots unshifted OK
[assert] RFdiffusion hotspot_res ['C4','C5','C7','C8'] verified against crop OK
contigs: A1-180/0 C1-10      190 residues (was 385)
```

β2m dropped from the **generation** target only (sits under the platform, not groove-facing).
Retained in full for **validation** — `structures/{g12d,wt,empty}.pdb` are the scoring targets.

---

## Stage 02 — RFdiffusion smoke test

| field | value |
|---|---|
| start | 2026-11-07 01:19 | | end | 2026-11-07 01:22 |
| GPU-h | 0.045 |
| inputs | `structures/crop/g12d_crop.pdb`, hotspots C4/C5/C7/C8, 65–85 aa |
| outputs | `pilot/rfd/smoke_{0,1}.{pdb,trb}`, `logs/rfd_smoke.log` |
| status | **pass** |

2 designs, 162 s wall including container start, IGSO3 cache build and model load.
Steady state **~62 s/backbone** ("Finished design in 1.04 minutes").

**Both designs contact all 4 primary hotspots.** Binder sits 3.55 Å off the peptide.

| design | binder len | target res contacted | peptide positions | hotspots hit |
|---|---|---|---|---|
| smoke_0 | 78 | 22 | 4,5,6,7,8 | **4/4** |
| smoke_1 | 70 | 20 | 4,5,6,7,8,9 | **4/4** |

Two things this surfaced:

1. **RFdiffusion re-chains and renumbers its output.** Input `A1-180/0 C1-10` returns as
   chain A = binder (1..N), chain B = heavy chain **and peptide merged**, renumbered
   N+1..N+190. So peptide p1–p10 land at B259–B268 for a 78 aa binder but B251–B260 for a
   70 aa one — **the offset moves with binder length.** All downstream code resolves
   positions through the `.trb` map (`src/comp5/rfmap.py`) and raises rather than guessing.
   This is the comp3 round-2 bug and it is live on this target.
2. **RFdiffusion output is backbone-only — N/CA/C/O, not even CB.** BSA and side-chain
   contact counts computed on it are not comparable to the Baker 1200 Å² / 31% reference,
   which comes from a side-chain-complete crystal structure. At a 4.5 Å side-chain cutoff a
   correctly docked design registered only 6 atom pairs; at 8 Å between backbone atoms it
   registered 606. `scripts/03_interface_metrics.py` now auto-detects side chains, widens the
   cutoff for backbone-only input, and emits `peptide_contact_fraction_area = not_run` rather
   than a number that would look like the gate metric but be wrong by an unknown factor.
   **Area metrics are only valid after MPNN + folding.**

### Throughput — the Stage 2 target does not fit

62 s/backbone × 3,500 backbones = **~60 GPU-h**. `docs/PIPELINE.md` budgets **~16 GPU-h**
for the H+7→H+26 block. The configured `n_backbones: 3500` is **~3.8× over budget**.
16 GPU-h at measured rate = **~930 backbones**. **Needs a decision — see report.**

---

## Stage 03 — pilot generation (32 designs)

| field | value |
|---|---|
| start | 2026-11-07 01:28 |
| inputs | as Stage 02, `inference.num_designs=32` |
| outputs | `pilot/rfd/p_*.{pdb,trb}`, `logs/rfd_pilot32.log` |
| status | *running at time of writing — see report for outcome* |

---

## Stage 03b — component smoke tests (run concurrently with pilot generation)

| field | value |
|---|---|
| start | 2026-11-07 01:40 | | end | 2026-11-07 02:10 |
| GPU-h | ~0.1 (shared with generation; Boltz timing below is CONTENDED) |
| status | **pass, 3 bugs found and fixed** |

### Template mmCIF — two defects, both fixed

`scripts/01_prep_targets.py:to_cif()` initially produced CIFs that Boltz would have
rejected. Both were already documented in `tfr1-comp4/scripts/make_template_cif.py` and
both were present here until ported:

- `full_sequence` was **empty** for every entity (gemmi leaves it empty with no SEQRES).
  Boltz's `parse_polymer` then indexes `sequence[j]` off the end and dies with a bare
  `IndexError`. Now filled from observed residues → `A:276, B:99, C:10`.
- `label_asym_id` was `Axp/Bxp/Cxp`, not `A/B/C`. Boltz matches `template_id` against
  `label_asym_id`, so `template_id: [A, B, C]` would silently not have matched. Forced to
  the author chain names.

Also recorded: structure chain A is 276 residues numbered **1–276 with no internal gaps**,
matching the config's 281-aa construct except a single terminal residue (276 is Pro in
7OW6, Leu in the construct). That residue sits in the α3 domain, 96 residues from the
groove. No design impact; noted rather than silently ignored.

### ProteinMPNN `--bias_by_res_jsonl` — silent data loss, fixed

Writing one JSON object per line loses every design but the last. The loader at
`ProteinMPNN/protein_mpnn_run.py:145` is

```python
for json_str in json_list:
    bias_by_res_dict = json.loads(json_str)   # REASSIGNS, does not merge
```

so designs other than the last die inside `tied_featurize` with `KeyError: <name>`.
Now written as **one JSON object on one line**. Caught only because the smoke test ran
two backbones rather than one.

### ProteinMPNN poly-Alanine — 75% of designs lost before folding

At `sampling_temp 0.1` with no Ala bias, **6 of 8** designs carried a poly-Ala run of
5–13 residues (`AAAAAAAAAAAAA` at position 41 in one case) and were rejected by the
stage-05 composition filter. Runs that long are not a filter artefact — they are
low-complexity surface with no core packing. This is ProteinMPNN's documented poly-Ala
habit on backbone-only helical input.

Measured grid (n=8 per row, neutral pool, 2 backbones — small sample):

| ala_bias / temp | hydrophobic | charged | max run | rejected |
|---|---|---|---|---|
| 0.0 / 0.1 | 0.60 | 0.31 | 13 | **75%** |
| −0.3 / 0.1 | 0.53 | 0.38 | 8 | 38% |
| −0.5 / 0.1 | 0.39 | 0.55 | 5 | 12% |
| −0.8 / 0.1 | 0.32 | 0.58 | 4 | 0% |
| −0.8 / 0.2 | 0.32 | 0.57 | 4 | 0% |

The 0.60 hydrophobic fraction in row 1 is **misleading** — it is mostly Ala, which is
nominally hydrophobic but too small to pack a core.

Charged fraction sits at 0.55–0.58 in every setting that controls poly-Ala, so the charge
is **inherent to SolubleMPNN on these backbones** (it is trained to raise surface charge to
avoid aggregation), not caused by the Ala bias. Against a ~0.25–0.30 reference for general
de novo minibinders this is high, and it is a developability risk to carry into Stage 5.

**Operating point set to `--ala-bias -0.6 --temp 0.1`**, applied to BOTH pools so it cannot
confound the neutral-vs-basic comparison. Deliberately not tuned further on n=8; the
distribution gets re-measured on the 32-backbone pilot.

### Positional R/K bias works and is not a no-op

9–10 binder positions per design fall within 12 Å of the p6 Cα, so the hypothesis is
actually applicable. Measured effect, neutral → basic: R+K **9–11 → 12–17** per design, for
a small MPNN score penalty (1.05–1.21 → 1.09–1.13).

Implemented as `--bias_by_res_jsonl` (per-position), **not** `--bias_AA_jsonl` (global).
The config asks for enrichment "at positions projecting over p6"; a global bias would make
the entire binder basic, wrecking pI and confounding the only comparison the two pools
exist to make.

### Boltz-2 three-state — runs, and the guardrail is load-bearing

One design, `g12d` state, 460 residues (78 binder + 281 + 99 + 10): **31 s inference,
83 s wall including model load**, while contended with RFdiffusion.

The per-pair extraction is **not optional**:

| quantity | value |
|---|---|
| `global_iptm` | **0.801** |
| binder ↔ peptide | 0.284 |
| binder ↔ HLA heavy | 0.256 |
| binder ↔ β2m | 0.237 |
| heavy ↔ β2m *(target-internal)* | **0.957** |
| heavy ↔ peptide *(target-internal)* | **0.924** |
| β2m ↔ peptide *(target-internal)* | **0.909** |

Global ipTM of 0.80 is almost entirely the MHC assembling with itself. Ranking on it would
have promoted every design equally and measured nothing about binding. This is the
four-chain version of the comp2 homotrimer problem and it is live here.

This design's binder pairs (0.24–0.28) are weak — expected for an un-triaged MPNN sequence
with no self-consistency filtering yet. It is a plumbing validation, not a hit.

---

## Stage 03 — pilot generation (32 designs) — COMPLETE

| field | value |
|---|---|
| start | 2026-11-07 01:28 | | end | 2026-11-07 02:03 |
| GPU-h | **0.58** (2076 s) |
| inputs | `structures/crop/g12d_crop.pdb`, hotspots C4/C5/C7/C8, contigs `A1-180/0 C1-10/0 65-85` |
| outputs | `pilot/rfd/p_*.{pdb,trb}` (32), `results/processed/backbone_metrics.csv` |
| status | **complete, rc=0, 32/32** |

**Measured throughput: 64.9 s/backbone** (2076 s / 32), matching the 62 s smoke estimate.

### Generation yield is far better than docs/PIPELINE.md predicts

PIPELINE.md Stage 2 says *"Expect only a small fraction to make extensive peptide contact —
that is the known behaviour of this protocol, not a bug."* Measured on the cropped target:

| primary hotspots contacted | backbones |
|---|---|
| 4 of 4 | **15 (47%)** |
| 3 of 4 | 10 (31%) |
| 2 of 4 | 7 (22%) |
| fewer than 2 | **0** |

**32/32 (100%) contact the peptide at all.** Binder lengths 66–85 aa, median 74 — all inside
the published 40–100 aa category 1. Median 24 target residues contacted.

This materially weakens the case for over-generating 3,500 backbones: the protocol is not
throwing away most of its output here, because the 190-residue domain crop leaves the binder
little else to bind.

### Per-position contact frequency tracks measured SASA

| pos | contacted | measured SASA Å² | class |
|---|---|---|---|
| p1 | **0%** | 1.2 | anchor (buried) |
| p2 | 9% | 2.1 | unclassified (buried) |
| p3 | 16% | 5.9 | unclassified (buried) |
| p4 | **47%** | 37.0 | primary |
| p5 | 84% | 26.9 | primary |
| p6 | 94% | 21.8 | tracked |
| p7 | **100%** | 75.5 | primary |
| p8 | 94% | 46.6 | primary |
| p9 | 78% | 15.7 | tracked |
| p10 | 56% | 6.7 | anchor (buried) |

Two things worth acting on:

1. **p4 is the outlier.** It is the third-most-exposed position (37.0 Å², 58% relative) but is
   contacted by only 47% of backbones, while p6 — deliberately *not* a hotspot — is contacted
   by 94%. The generated binders centre on p5–p9, not p4–p8. If Stage 2 wants the bulge
   N-terminus engaged, the hotspot set may be better as `[5,6,7,8]`; if not, p4 is costing a
   conditioning slot for little return. Flagged, not changed — the hotspot set is a config
   decision.
2. **p6 at 94% is the desired behaviour, not a violation.** Binders sit over the mutant Asp
   without being able to read its side chain (10.4% exposed, salt-bridged to Arg114). That is
   precisely the indirect-specificity mechanism `docs/CHALLENGE_SPEC.md §7` describes, and it
   is what the WT differential will have to detect.
3. p1 at 0% and p2/p3 at 9/16% confirm the anchor exclusions are correct and that nothing is
   wasting interface on buried positions. p10 at 56% is the groove rim near the F-pocket, not
   the buried Lys side chain itself.

---

## Stage 04 — ProteinMPNN, two pools, 34 backbones

| field | value |
|---|---|
| start | 2026-11-07 02:09 | | end | 2026-11-07 02:28 |
| inputs | `pilot/rfd/*.pdb` — **34** backbones (32 `p_*` + the 2 `smoke_*`, identical settings) |
| outputs | `pilot/mpnn/{neutral,basic}/seqs` (136 sequences each) |
| status | complete |

Operating point `--ala-bias -0.6 --temp 0.1`, applied to both pools.

### Composition, measured on n=136 per pool (not the n=8 tuning grid)

| | neutral | basic | reference |
|---|---|---|---|
| hydrophobic fraction (median) | 0.304 | 0.286 | 0.30–0.40 |
| charged fraction (median) | **0.600** | **0.609** | 0.25–0.30 |
| R+K count (median) | 18 | **20** | — |
| net charge (median) | −7 | −6 | — |
| composition rejects | 19.1% | 21.3% | — |

- **Poly-Ala is gone**: zero hydrophobic runs in 272 sequences, against 75% rejection before
  the Ala de-bias. The fix worked.
- **It relocated rather than vanished.** All remaining rejects are homopolymer runs of
  **Glu (43 occurrences)** and **Arg (18)** — longest `EEEEEE` and `RRRRRR`. Suppressing Ala
  pushed ProteinMPNN into charge runs instead.
- **Charged fraction 0.60 is the real quality concern**, roughly double the reference for de
  novo minibinders. Combined with net charge −6/−7 on a ~74 aa binder, these are polar,
  weakly-packed sequences. SolubleMPNN is trained to raise surface charge, so some of this is
  by design, but not this much.
- **The basic pool's bias is measurable but weak**: R+K median 18 → 20, net charge −7 → −6.
  A +2 shift on a net −6 scaffold does not plausibly test the R50α/K70α electrostatic
  hypothesis. If that hypothesis is worth testing, `bias_weight` (1.2) and/or `bias_radius`
  (12 Å) need raising — **flagged, not changed**, since it is a config decision and the pools
  must stay comparable.

**Caveat on net charge direction.** The discrimination mechanism in
`docs/CHALLENGE_SPEC.md §7` is that the buried Asp makes the HLA surface *more negative*, and
the reference TCR reads that with Arg50α/Lys70α. A binder at net −6 is electrostatically
mismatched for reading a negative surface. This does not invalidate the designs — shape
complementarity may carry the interface — but it is the opposite of what the stated mechanism
suggests, and it should be checked against the differential rather than assumed either way.

---

## Stage 05 — selection

| field | value |
|---|---|
| status | complete, after one crash (below) |

| pool | sequences | pass composition | selected |
|---|---|---|---|
| neutral | 136 | 110 | 16 |
| basic | 136 | 107 | 16 |
| decoy (neutral) | 32 | 22 | 8 |

40 designs total. Rejections recorded in `results/rejected/select_rejected.csv`
(per-backbone cap and below-top-N are *not* quality failures and are labelled separately).

### Self-inflicted failure, recorded because it cost a pipeline abort

`Path.relative_to(ROOT)` raises `ValueError` when the path is relative and ROOT absolute.
I hit this **four separate times** (`03`, `05`, `06`, then `04b`). The fourth one fired
*after* selection had successfully written all 32 fastas, and `set -e` killed the pilot at
step 3 of 7 — losing nothing but requiring a restart.

Fixes applied:
- `src/comp5/paths.py:rel()` now wraps it; **all** `relative_to(ROOT)` call sites swept.
  A regex sweep mangled three nested call sites (`str(x.relative_to(ROOT))` →
  `rel(str(x, ROOT))`), which a per-script re-run caught and I repaired by hand. Every script
  in `scripts/` was then re-executed to confirm it still runs.
- `scripts/10_pilot.sh` now **skips any stage whose output already exists**, so a late crash
  no longer re-burns GPU on generation or MPNN. This is what CLAUDE.md §1 asks for and it
  should have been there from the start.

### Two id bugs caught before they reached the fold

1. `read_fastas` appended `_s<N>` to a stem that already carried one, producing
   `p_0_s1_s1` and breaking the backbone join in `06_manifest.py`.
2. **Worse:** both pools design the *same* backbones, so `neutral/p_0_s1.fa` and
   `basic/p_0_s1.fa` collided on one dict key — one pool would have silently overwritten the
   other and every `mpnn_pool` label downstream would have been wrong. Ids are now namespaced
   (`neutral_p_0_s1`, `basic_p_31_s2`).

Both were found by reading the code while generation ran, not by a test. Neither would have
raised an error — they would have produced a plausible-looking manifest.

---

## Stage 06 — three-state Boltz-2 and GATE 1

| field | value |
|---|---|
| start | 2026-11-07 02:33 | | end | 2026-11-07 03:12 |
| GPU-h | **0.58** (1587 s main + 508 s decoy) |
| inputs | 32 peptide-targeted + 8 β2m framework decoys |
| outputs | `pilot/boltz/{main,decoy}_three_state.csv`, `results/processed/manifest*.csv/.parquet` |
| status | **complete — 120/120 predictions, all `ok`** |

Measured **16.5 s/prediction** (main, 96 predictions), 21.2 s (decoy, 24 — model load
amortises over fewer). Earlier in-run estimates of 49 s were wrong: Boltz writes outputs in
a batch at the end, so counting files mid-run badly understates throughput.

### GATE 1 VERDICT — 2 PASS, 2 FAIL

| criterion | result |
|---|---|
| 1. runs end to end, no manual intervention | **PASS** (after the `set -e` fix) |
| 2. some designs reach `peptide_contact_fraction` ≥ 0.25 | **FAIL** — 0/32 |
| 3. three-state differential discriminates at all | **FAIL** — see below |
| 4. throughput fits the Stage 3 matrix | **PASS** |

### Criterion 2 — the gate is above its own ceiling

| | median | range |
|---|---|---|
| `bsa_total_target` | 849.2 Å² | 245.3–1554.6 |
| `bsa_peptide` | 68.3 Å² | 0.0–215.2 |
| `peptide_contact_fraction_area` | 0.095 | 0.000–**0.247** |

Max is 0.247 against a 0.25 gate. This is not a design failure — it is the threshold being
set above what the target physically allows (239.4 Å² of peptide; ceiling 0.199 at the
Baker-reference 1200 Å² interface). **The threshold must be recalibrated from negatives.**

Median total interface of 849 Å² is healthy and in the expected minibinder range.
**11 of 32 designs have zero peptide burial** after folding — i.e. a third of backbones that
looked docked drifted off the peptide once a sequence was threaded and the complex refolded.
That gap between backbone geometry and folded pose is itself a useful filter.

### Criterion 3 — ipTM cannot grade peptide dependency, but geometry can

**The two pair scores are not independent.** binder↔peptide vs binder↔mhc ipTM correlate at
**r = 0.9956**, constant ratio ≈ 0.92. `pair_chains_iptm` reports one docking-confidence
number spread across pairs; it does **not** resolve which part of the target is engaged. Any
plan to use the pep/mhc ratio as a peptide-focus metric is invalid.

Clean comparison, controlling for decoys that drifted onto the peptide:

| group | n | `iptm_margin_empty` median |
|---|---|---|
| main, touching peptide | 22 | +0.032 |
| main, **zero** peptide burial | 10 | −0.219 |
| decoy, **true framework** (0 peptide burial) | 5 | **−0.002** |
| decoy, drifted onto peptide | 3 | +0.033 |

Separation between real peptide-binders and true framework decoys: **+0.034 — nothing.**
By `bsa_peptide` quartile the relationship is **non-monotonic** (−0.223, +0.056, −0.091,
+0.125), so it is not a dose–response.

Read correctly: the decoys sit at −0.002, which is *exactly* the ≈0 a framework binder should
give — the control behaves as designed. The failure is that peptide-engaging designs also sit
near 0, so `iptm_margin_empty` only flags "ignores the peptide region entirely" and cannot
grade dependency.

**What does work:** geometry. `bsa_peptide` median 68.3 Å² (peptide-targeted) vs 0.0 Å²
(decoys); 21/32 vs 3/8 contact the peptide at all; 9/32 vs 1/8 hit all four hotspots.

### Criterion 4 — fits, for ONE model

- Stage 3 Boltz leg, 600 designs × 3 states = 1,800 predictions ≈ **8.3 GPU-h**
- Triage pass (~1,500, g12d only) ≈ **6.9 GPU-h**
- **≈15.3 GPU-h** against a ~20 GPU-h block → fits, ~4.7 GPU-h spare.

That spare will not cover the Chai-1 leg (≈8 GPU-h), which the ranking's "differentials
across ≥2 models" requirement depends on — and Chai-1 is not installed.

### Not yet interpretable

Margins sit near zero (`margin_wt` median +0.028, sd 0.157; `margin_empty` median +0.015,
sd 0.322) and **there is no run-to-run variance estimate**, so none of that spread can be
called signal. The seed replicate is the first thing the overnight run measures. Until it
lands, no design should be called a hit and no margin threshold should be set.

---

## Stage 07 — seed-replicate noise floor — **THE MOST IMPORTANT PILOT RESULT**

| field | value |
|---|---|
| start | 2026-11-07 03:10 | | end | 2026-11-07 03:22 |
| GPU-h | 0.14 (8 designs × 3 states, seed 1234 vs the pilot's unseeded run) |
| outputs | `pilot/boltz_rep/rep2_three_state.csv` |
| status | complete |

### `iptm_margin_wt` is indistinguishable from noise

| design | margin_wt seed A | seed B | \|Δ\| |
|---|---|---|---|
| basic_smoke_1_s1 | −0.353 | **+0.185** | **0.538** |
| neutral_p_21_s4 | +0.056 | −0.240 | 0.295 |
| neutral_p_24_s2 | +0.263 | +0.021 | 0.242 |
| neutral_p_6_s1 | +0.063 | −0.080 | 0.143 |
| neutral_p_22_s2 | +0.268 | +0.145 | 0.123 |
| neutral_p_4_s2 | −0.038 | −0.117 | 0.078 |
| basic_p_29_s4 | −0.312 | −0.256 | 0.056 |
| basic_p_14_s1 | −0.029 | −0.016 | 0.013 |

| metric | within-design seed noise (mean \|Δ\|) | between-design spread (sd, n=32) | **S/N** |
|---|---|---|---|
| `iptm_margin_wt` | **0.186** | 0.157 | **0.84×** |
| `iptm_margin_empty` | 0.173 | 0.322 | 1.86× |

**The seed-to-seed variance of a single design is LARGER than the spread across all 32
designs.** Ranking on `iptm_margin_wt` ranks the random number generator. `basic_smoke_1_s1`
moved from the worst design in the set to above the median on a seed change alone;
`neutral_p_24_s2` went from joint-best (+0.263) to nothing (+0.021).

`iptm_margin_empty` at 1.86× is better but still too weak to gate on without replication.

### Why this outranks every other pilot finding

**Mutation specificity is objective 1 and ranks FIRST** in the official lexicographic order
(`docs/CHALLENGE_SPEC.md §4`). `docs/PIPELINE.md` Stage 4 names `iptm_margin_wt` as its proxy,
required "in both models". The pilot says that proxy carries no usable signal at one sample
per design. **The ranking plan does not currently work**, independent of the disk, the
schedule, or Chai-1.

This also retrospectively invalidates the Gate 1 criterion-3 analysis: the +0.034 "clean
separation" and the non-monotonic quartiles were never going to be interpretable against a
noise floor of 0.17–0.19. The correct reading is not "the differential fails" but
**"the differential was never measured"** — one sample per design is below the resolution
needed to measure it at all.

### Two ways out, and what is being measured tonight

| option | cost | status |
|---|---|---|
| (a) average *n* diffusion samples in ONE run | model loads once; ~√n noise cut if samples are independent | **`scripts/21_sampletest.sh`, chained after the overnight run** |
| (b) *n* independent seeded runs | n× full cost — at n=8 this is ~66 GPU-h for Stage 3 alone, infeasible | not attempted |
| (c) drop `iptm_margin_wt` as the objective-1 proxy | free, but leaves objective 1 with no predictor | **decision for Harjas** |

`scripts/21_sampletest.sh` runs 12 designs spanning the margin range at
`diffusion_samples=5 --average-samples`, **twice with different seeds**, so the noise of the
*averaged* estimator is measured exactly as the single-sample noise was. Verdict threshold
set at S/N > 2.

If (a) works, Stage 3 stays affordable and the ranking survives. If it does not, objective 1
has no model-based proxy and the honest options are to rank on objective 2 (peptide
dependency, where geometry does discriminate: 68.3 vs 0.0 Å² buried) and declare objective 1
unpredicted — which is defensible, since `docs/PIPELINE.md §3` already says "no method
reliably predicts the WT-vs-G12D window for a de novo binder".

---

## Stage 08 — definitive noise floor, n=32 (REVISES Stage 07)

| field | value |
|---|---|
| start | 2026-11-07 05:27 | | end | 2026-11-07 05:54 |
| GPU-h | 0.44 (32 designs × 3 states, seed 1234) |
| outputs | `pilot/boltz_rep/rep32_three_state.csv` |
| status | complete |

The n=8 estimate in Stage 07 **over-stated the noise for `iptm_margin_wt`**. Corrected figures:

| metric | noise mean \|Δ\| | between-design sd | **S/N** | **Spearman ρ between seeds** |
|---|---|---|---|---|
| `iptm_margin_wt` | 0.136 (was 0.186) | 0.157 | **1.16×** | **+0.239** |
| `iptm_margin_empty` | 0.136 (was 0.173) | 0.322 | **2.37×** | **+0.853** |
| `g12d_iptm_binder_peptide` | 0.116 | 0.251 | 2.17× | +0.703 |

### The rank correlation is the number that matters, not S/N

We rank designs; we do not use the absolute values. So the question is whether the *ordering*
survives a seed change.

- **`iptm_margin_wt`: ρ = +0.239.** Rank the 32 designs by it under one seed and re-rank under
  another and you get a nearly unrelated ordering. **Confirmed unusable for ranking**, even
  though S/N rose to 1.16×. The Stage 07 conclusion stands; only the magnitude changes.
- **`iptm_margin_empty`: ρ = +0.853.** Strong rank agreement. This is **more usable than
  Stage 07 claimed** and it clears the S/N > 2 bar. Demoting it to "weak tiebreak" was too
  harsh.

### Reproducible is not the same as informative — the distinction that matters here

`iptm_margin_empty` is now shown to be **reproducible** (ρ 0.85). It was separately shown in
Stage 06 **not to discriminate** framework binders from peptide binders (main-touching-peptide
+0.032 vs true-framework decoys −0.002, separation +0.034).

Both are true. It reliably measures *something* about a design; that something is not peptide
dependency. A reproducible metric that does not track the property you care about is a
precise instrument pointed in the wrong direction, and it is more dangerous than a noisy one
because the stability invites trust.

**Consequence for ranking:** the geometry-first order stays. `iptm_margin_empty` is promoted
from "tiebreak" to a legitimate ranking input **for interface confidence**, but it must not be
presented as evidence of peptide dependency — `bsa_peptide` and hotspot coverage are the only
measured discriminators of that (68.3 vs 0.0 Å²).

`iptm_margin_wt` remains recorded, not ranked on. `scripts/21_sampletest.sh` will show whether
diffusion-sample averaging lifts its ρ; S/N alone is not the test — **ρ > 0.7 is.**
