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

---

## Stage 09 — shuffled-sequence controls: the first clean positive

| field | value |
|---|---|
| start | 2026-11-07 05:54 | | status | g12d complete; wt/empty in progress |

32 real designs vs 32 **composition-matched** scrambles (every control has the exact same
amino-acid content as its parent, asserted in `scripts/07_controls.py`), `g12d` state:

| | n | median | range |
|---|---|---|---|
| real designs | 32 | **0.462** | 0.125–0.962 |
| shuffled controls | 32 | **0.203** | 0.087–0.565 |

Separation **+0.260**. Only **2 of 32** scrambles reach the real median.

**The absolute binder↔peptide ipTM is sequence-specific.** It is not merely reporting "a
75-residue chain of this composition fits in this groove" — if it were, composition-matched
scrambles would score the same, and they do not. The ProteinMPNN sequences are contributing
real signal, and `g12d_iptm_binder_peptide` now has a measured negative floor rather than a
guessed one.

### How this reconciles with the margin failures

The three results are consistent once absolute scores and differences of scores are kept
apart:

| quantity | reproducible? | discriminates? | usable for ranking |
|---|---|---|---|
| `g12d_iptm_binder_peptide` (absolute) | ρ +0.703 | **yes**, +0.260 vs scrambles | **yes** |
| `iptm_margin_empty` (difference) | ρ +0.853 | no, +0.034 vs framework decoys | confidence only |
| `iptm_margin_wt` (difference) | ρ +0.239 | not measurable at this noise | **no** |
| `bsa_peptide` (geometry) | deterministic | **yes**, 68.3 vs 0.0 Å² | **yes** |

Differencing two ipTM values subtracts most of the signal and adds both their noises, which is
why the absolute score survives and the margins do not. Nothing here contradicts Stage 08 —
the margins remain unfit for ranking while the score they are built from is fine.

### Resulting ranking scheme

1. **Gate on geometry** — `bsa_peptide` / hotspot coverage: is the binder on the peptide at all?
   (The only measured discriminator of peptide focus.)
2. **Gate on calibrated absolute confidence** — `g12d_iptm_binder_peptide` above the shuffled-control
   percentile: is the interface real rather than a composition artefact?
3. **Rank within survivors** by interface size and contact density.
4. **Record but do not rank on** `iptm_margin_wt` (noise) and `iptm_margin_empty` (reproducible
   but not discriminating). Neither may be cited as evidence of mutation specificity or
   peptide dependency in the methods note.

---

## Stage 10 — the two-gate filter: 0/32 scrambles survive, 10/32 real designs do

Cuts taken from the **negative** distributions, not chosen in advance:
`iptm_binder_peptide_g12d > 0.499` (shuffled-control p95) **AND**
`bsa_peptide > 97.2 Å²` (β2m-decoy p95).

| set | n | ipTM gate only | geometry gate only | **both gates** |
|---|---|---|---|---|
| real designs | 32 | 14 | 14 | **10 (31%)** |
| shuffled controls | 32 | 3 | 5 | **0** |
| β2m framework decoys | 8 | 2 | 1 | 1 |

### Why the conjunction is necessary, and why neither gate alone would do

Each gate has a different failure mode, and they do not overlap:

- **ipTM alone** admits framework binders — the β2m decoys reach p95 0.506, essentially the
  same as the shuffled p95 of 0.499. A confident interface is not a peptide interface.
- **Geometry alone** admits scrambles — shuffled `bsa_peptide` has median 2.95 Å² but
  **p95 212.5 Å², above the real designs' p95 of 204**. Boltz docks a scrambled chain
  somewhere, and once it lands on the peptide the buried area is set by shape and size, not by
  sequence quality. A p95-calibrated cut on geometry alone would have reproduced the original
  0.25-contact-fraction mistake: a threshold above where the real designs live.

Together they are complementary: geometry asks *is it on the peptide*, confidence asks *is the
interface real rather than a composition artefact*. **All 32 composition-matched scrambles are
eliminated** while 31% of real designs survive.

### What this does and does not do to the Gate 1 verdict

- Criterion 2 is **satisfied in substance**: 10/32 designs clear a *calibrated* peptide-contact
  gate. The original 0.25 figure was unusable because it sat above the physical ceiling, not
  because the designs failed.
- Criterion 3 **still fails as literally written** — the g12d-minus-empty ipTM spread does not
  separate peptide-targeted from framework designs, and no amount of calibration changes that.
  What is now true is that its *purpose* — telling peptide binders from framework binders — is
  achieved by other means, with a measured 0% false-positive rate on sequence controls.

That distinction is worth keeping precise rather than declaring Gate 1 retroactively passed.

### Projected yield for Stage 2

900 backbones × 4 sequences × 2 pools, ~80% composition pass, 1 per backbone per pool
→ ~1,800 candidates → at the pilot's 31% two-gate rate, **~550 survivors** for a 20-design
submission. Ample margin, and the constraint shifts from yield to diversity.

---

## Stage 11 — formal calibration (CORRECTS Stage 09 and Stage 10)

| field | value |
|---|---|
| end | 2026-11-07 07:15 | GPU-h | ~1.8 total for the overnight block |
| outputs | `configs/thresholds_calibrated.yaml`, `results/processed/calibration_report.json`, `logs/calibration.log` |
| inputs | 32 real · 64 sequence controls · 21 framework decoys · 32-design seed replicate |

### Correction: Stages 09–10 calibrated on the weaker control family

Stage 09 reported "+0.260 separation, the absolute ipTM is sequence-specific" and Stage 10
reported "0/32 scrambles survive, 10/32 real pass". **Both used shuffled controls only.**
With the reversed family included:

| family | n | median | p95 | max |
|---|---|---|---|---|
| real designs | 32 | 0.462 | 0.894 | 0.962 |
| shuffled | 32 | 0.221 | 0.499 | 0.598 |
| **reversed** | 32 | **0.379** | **0.777** | **0.849** |
| β2m decoys | 21 | 0.234 | 0.565 | 0.570 |

**Reversal is a much stronger negative than shuffling.** It preserves composition *and* local
run-length and hydrophobic patterning, and Boltz rewards that structure: reversed controls
reach ipTM 0.849. Median separation from real designs is **+0.083**, not the +0.260 claimed
against shuffled alone. Including it was the right call — a control family that is too easy
flatters the metric and would have set the threshold far too low.

Revised two-gate performance at the calibrated cuts (ipTM > 0.7746, `bsa_peptide` > 97.2 Å²):

| set | pass both |
|---|---|
| real designs | **5/32 (16%)** |
| shuffled | 0/32 |
| reversed | **2/32** |
| β2m framework decoys | 0/21 |

So **2 of 64 sequence controls survive — a 3% false-positive rate, not 0%** — and real yield
is 16%, not 31%. The filter still works, and still eliminates every framework decoy, but it is
weaker than Stage 10 claimed.

Real designs do separate in the **upper tail**, which is where selection happens: real p95
0.894 vs reversed p95 0.777. Selecting the top of the distribution remains meaningful.

### Calibrated cuts and yields

| metric | cut | source | real designs passing |
|---|---|---|---|
| `bsa_peptide` | 97.2 Å² | decoy p95 | 14/32 |
| `n_peptide_contacts` | 4 | decoy p95 | 12/32 |
| `iptm_binder_peptide_g12d` | 0.7746 | control p95 | 8/32 |
| `iptm_g12d` | 0.7009 | control p95 | 8/32 |
| `peptide_contact_fraction_area` | 0.1751 | decoy p95 | 5/32 |
| `iptm_margin_empty` | 0.3352 | control p95 | 4/32 |
| **`iptm_margin_wt`** | **0.4041** | **seed-noise p95** | **0/32** |

`iptm_margin_wt` is now dead by construction: the calibrator detected that the sequence-control
p95 (0.399) sat *below* the seed-replicate noise p95 (0.404), raised the cut to the noise
floor as designed, and **no real design clears it**. That is the correct outcome — a threshold
beneath its own noise floor would have been measuring the RNG.

Note also `peptide_contact_fraction_area` calibrates to 0.1751, against the original
provisional 0.25 that sat above the physical ceiling. 5/32 pass. The provisional value was not
merely wrong, it was unreachable.

### Gate 1 re-run: still FAIL

Criterion 3 fails as written, unchanged and uncalibratable — the g12d-minus-empty ipTM spread
does not separate peptide-targeted from framework designs. Criterion 2 now passes in substance
(5/32 clear the calibrated contact gate; 14/32 clear `bsa_peptide`).

### Revised Stage 2 projection

~1,800 candidates × 16% two-gate survival ≈ **290 survivors** for a 20-design submission
(was ~550 on the over-optimistic rate). Still ample; diversity remains the binding constraint.

---

## Stage 11b — operational failure: self-match deadlock in the chained waiters

**~30 minutes of idle GPU**, 17:15–17:45, between the overnight block finishing and the sample
test being restarted by hand.

The chain scripts gated on:

```bash
while pgrep -f "20_overnight.sh" >/dev/null 2>&1; do sleep 60; done
```

`pgrep -f` matches **full command lines**. The harness shell wrappers that created those
scripts carry the heredoc body — including the literal string `20_overnight.sh` — in their own
command lines, and they outlive the command. So each waiter was blocked on a process whose
command line contained the name of the thing it was waiting for, and would have waited
indefinitely. `chain_stage2` was additionally gated on `chain_sampletest`, so the whole
remaining queue was stuck behind it.

The same trap bit the cleanup: `pkill -f "chain_sampletest|chain_stage2"` matched the shell
executing it and killed itself (exit 144).

**Rule going forward:** do not gate background work on a `pgrep -f` pattern that can appear in
the command line of the process doing the matching. Long jobs are launched directly through
the harness background runner, which tracks them and notifies on exit — the 32-backbone
generation and all Boltz stages ran that way without trouble. The `nohup` + `pgrep` wrapper
pattern bought nothing and cost half an hour.

No data was lost; the overnight block had already completed and committed before the stall.

---

## Stage 12 — diffusion-sample averaging RESCUES objective 1

| field | value |
|---|---|
| start | 2026-11-07 07:45 | | end | 2026-11-07 08:32 |
| GPU-h | 0.80 (12 designs × 3 states × 5 samples, run twice at different seeds) |
| outputs | `pilot/boltz_samp/samp5{,b}_three_state.csv`, `logs/sampletest.log` |
| status | **PASS — both margins now usable** |

| metric | noise mean \|Δ\| | Spearman ρ | S/N |
|---|---|---|---|
| `margin_wt` 1 sample | 0.136 | +0.239 | 1.15× |
| **`margin_wt` 5 samples** | **0.053** | **+0.846** | **2.94×** |
| `margin_empty` 1 sample | 0.136 | +0.853 | 2.37× |
| **`margin_empty` 5 samples** | **0.073** | **+0.881** | **4.40×** |

Noise fell **2.57×** on 5 samples against √5 = 2.24 expected for independent draws — slightly
better than independent, so the samples are not correlated in a way that defeats averaging.

**Objective 1 is no longer unpredicted.** `iptm_margin_wt` is restored to the ranking, with a
hard condition: **valid only at ≥5 averaged samples.** Any single-sample margin in this repo
is noise. The manifest carries `n_samples_averaged` per row so this cannot be lost track of.

### The cost works out because samples share the trunk

Measured **39.7 s/prediction at 5 samples vs 16.5 s at 1 — a 2.41× cost for 5× the samples**,
because Boltz recomputes only the diffusion module per sample while the trunk is shared. Naive
5× costing would have made this look impossible.

| stage | design × states | @1 sample | @5 samples |
|---|---|---|---|
| triage (g12d only) | 1500 × 1 | 6.9 GPU-h | 16.6 GPU-h |
| full matrix | 600 × 3 | 8.2 GPU-h | 19.9 GPU-h |
| **total** | | 15.1 | **36.4 GPU-h** |

36.4 h does not fit the ~20 h block, so samples are spent only where they buy something:

- **Triage at 1 sample.** It ranks on absolute `iptm_binder_peptide` (ρ +0.703 unaveraged) and
  geometry; no margin is consulted, so averaging would buy nothing. → 8.3 GPU-h at ~1,800
  candidates.
- **Full matrix at 5 samples**, cut from 600 to **450** designs. → 14.9 GPU-h.
- **Total ≈ 23 GPU-h**, against ~60 h remaining to the freeze with Stage 2 (16 h) also to run.

### Outstanding, and it matters

`configs/thresholds_calibrated.yaml` sets `iptm_margin_wt > 0.4041`, derived from
**single-sample** noise. At 5 samples the noise floor is ~2.6× lower, so that cut is far too
strict and currently passes 0/32 designs. **It must be recalibrated against 5-sample negative
controls before it gates anything.** Scheduled with the Stage 3 triage, where 5-sample
controls become available as a by-product.

---

## Stage 13 — 1 TB disk mounted; bulk outputs relocated

Harjas ran the `mkfs.ext4` + mount by hand (the sandbox blocks device formatting from here,
correctly). `/dev/sdb` — a 1000 GB GCP PersistentDisk provisioned 2026-11-06, verified zeroed
and signature-free before formatting — is now **984 G mounted at `/mnt/comp5data`**, writable.

Idle Boltz output trees moved there and symlinked back, so every relative path in the scripts
keeps working untouched:

| moved | size |
|---|---|
| `pilot/boltz` | 659 M |
| `pilot/boltz_rep` | 225 M |
| `pilot/boltz_samp` | 637 M |

Verified after the move: CSVs read through the symlinks, `scripts/11_gate1.py` still resolves
its inputs, and Stage 2 generation was unaffected. `pilot/af2_all` was deliberately left in
place because it was mid-run.

**fstab entry added by Harjas 2026-11-07** — `UUID=812876ce.../mnt/comp5data ext4 defaults,noatime,nofail 0 2`, UUID verified against the device and `nofail` set so a detached disk cannot block boot. The mount is now permanent.
persistence. The mount is live for this session, which covers the run to Monday, but it will
not survive a reboot. To make it permanent:

```
! UUID=$(sudo blkid -s UUID -o value /dev/sdb); echo "UUID=$UUID /mnt/comp5data ext4 defaults,noatime,nofail 0 2" | sudo tee -a /etc/fstab
```

Disk pressure is no longer a constraint: 25 G free on `/` plus 982 G on `/mnt/comp5data`,
against a Stage 3 requirement of ~11 G.

---

## Stage 14 — AF2 screen on all 272 pilot sequences, and a failed validation attempt

### Pass rate is 0.7%, not the 3.1% estimated from n=32

| | n=32 (pre-selected) | **n=272 (all MPNN sequences)** |
|---|---|---|
| pass all criteria | 1 (3.1%) | **2 (0.7%)** |
| median `pae_interaction` | 27.03 | 26.78 |
| median `target_aligned_rmsd` | 40.7 | **43.6** |
| failure mode | 30/32 `pae>=10` | **269/272 `pae>=10`** |

The n=32 figure was optimistic because those 32 were selected by MPNN score, and the whole
point of Stage 14 is that MPNN score does not predict AF2 success. 0.7% is the honest number.

Both passers — `neutral_p_3_s2` (pae 4.92) and `basic_p_3_s2` (pae 5.21) — come from the
**same backbone, `p_3`**. So 1 backbone in 34 produced anything viable, and the two hits are
not independent. Clustering by backbone is itself evidence the signal is structured rather
than random noise.

Composition correlates with AF2 interface quality, monotonically across quartiles:

| quartile by pae | median pae | hydrophobic | charged |
|---|---|---|---|
| best | 25.22 | 0.324 | 0.566 |
| Q2 | 26.53 | 0.304 | 0.587 |
| Q3 | 26.99 | 0.286 | 0.611 |
| worst | 27.30 | 0.256 | 0.671 |

Hydrophobic ≥0.38 → median pae 24.15; <0.30 → 26.89. The direction is consistent and argues
for dropping SolubleMPNN, whose 0.52–0.67 charged fraction is roughly double the typical de
novo minibinder. But the effect is small and **every quartile still fails**, so composition is
a contributing factor, not the cause.

### The positive control I built was invalid — my error

Before committing ~20 GPU-h to the full screen I tried to validate the setup by running
AF2-ig on the real crystallographic TCR–pMHC interface from 7OW6 (TCR chains D+E merged as
the binder, platform+peptide as the target). Result:

```
pae_interaction 28.11   plddt_binder 29.2   binder_aligned_rmsd 28.78   plddt_target 94.23
```

A validated binder scoring as badly as our worst designs would mean the setup is broken. It
does not mean that, because **the control is not a fair test and that is my mistake**:
`plddt_binder` of 29.2 says AF2 could not fold the 426-residue two-domain TCR at all from a
single sequence with no MSA. AF2-initial-guess works precisely because the binder is a small
de novo miniprotein foldable single-sequence while the target is templated. A TCR is outside
that regime, so the number says nothing about our designs.

What the control *did* confirm: `plddt_target` 94.23, so target templating via the
superposed side-chain crop is working correctly.

### The correct control, now running

Composition-matched **shuffled** sequences threaded onto the same backbones (n=32). The
question is whether AF2 separates real designs from scrambles on this target:

- if scrambles score clearly worse than 26.8, `pae_interaction` has signal here and our
  designs are simply mediocre — the two passers are real and the screen is worth running;
- if scrambles score the same ~27, `pae_interaction` is saturated for this target and must
  **not** be used as the primary gate, in which case the pipeline falls back to the calibrated
  Boltz geometry + confidence scheme.

**The 20 GPU-h full screen is held until this returns.** Committing to it on an unvalidated
metric is exactly the mistake the earlier `iptm_margin_wt` work was supposed to teach.

---

## Stage 15 — the correct negative control: AF2-ig setup VALIDATED, screen released

32 composition-matched shuffled sequences threaded onto the same backbones, vs the 272 real
designs:

| metric | real (n=272) | scramble (n=32) | separation |
|---|---|---|---|
| `pae_interaction` median | 26.78 | 27.14 | **+0.36** |
| `pae_interaction` best | **4.92** | **19.56** | — |
| `pae_interaction < 10` | **3/272** | **0/32** | — |
| `plddt_binder` | **87.9** | 68.8 | **−19.1** |
| `binder_aligned_rmsd` | **2.03** | 23.30 | **+21.3** |
| `target_aligned_rmsd` | 43.6 | 45.9 | +2.27 |

### The setup works; `pae_interaction` is stringent, not broken

`plddt_binder` (−19.1) and `binder_aligned_rmsd` (+21.3) separate real sequences from
composition-matched scrambles decisively. Real sequences fold into the backbone they were
designed for; scrambles do not. **AF2-ig is functioning correctly on this target**, which is
what the invalid TCR control in Stage 14 failed to establish.

`pae_interaction` behaves differently from the others and the distinction matters:

- **Saturated in the bulk** — real 26.78 vs scramble 27.14 is no separation at all. Judged on
  medians alone the metric looks useless.
- **Discriminating in the tail** — 0 of 32 scrambles reach <10 and the best scramble manages
  only 19.56, while 3 of 272 real designs clear it.

Both statements are true, and only the second is relevant to how the metric is used. We are
not averaging `pae_interaction`, we are selecting the extreme tail of it. A median comparison
would have wrongly condemned the gate. This is the same lesson as Stage 08 in a different
guise: match the statistic to how the number is actually consumed.

**Decision: the full AF2 screen is released.** `pae_interaction < 10` stays the hard gate,
with the honest caveat that it is stringent and the expected yield is ~0.7–1.1%.

### A second, cheaper filter emerges

`binder_aligned_rmsd` and `plddt_binder` discriminate strongly and are available for free in
the same scorefile. They give a graded ranking among designs that fail `pae<10`, which matters
because at 0.7% yield we may not fill 20 slots on the hard gate alone. Fallback order, to be
used only if the hard gate under-fills:

1. `pae_interaction < 10` — primary, published, tail-validated here
2. `binder_aligned_rmsd < 1.5` **and** `plddt_binder > 80` — fold self-consistency, separation
   +21.3 / −19.1 against scrambles
3. Boltz three-state differential at 5 samples on whatever survives

### Expected yield for Stage 2

900 backbones × 8 sequences = 7,200 candidates. At the measured 0.7–1.1%: **~50–80 survivors**,
clustered on perhaps 25–30 distinct backbones (both pilot passers shared backbone `p_3`).
Enough for 20 submissions, thin for diversity — the topology monoculture noted in Stage 16
remains the larger worry.

---

## Stage 16 — Stage 2 generation complete; timeline becomes the binding constraint

| field | value |
|---|---|
| end | 2026-11-07 21:33 Melbourne |
| GPU-h | ~16.5 (900 backbones at 64.9 s, minus the restart) |
| outputs | `gen/rfdiff/g_*.{pdb,trb}` — **900 backbones, 129 MB total** |

The `write_trajectory=False` fix held: zero `traj` directories, 129 MB for 900 backbones
against the 6.0 GB that 532 backbones had consumed before it.

### Remaining budget, measured not estimated

~44.5 h to the Monday 18:00 freeze.

| stage | cost | basis |
|---|---|---|
| beta diversity arm, 150 designs | ~4.5 h | measured 107 s/backbone |
| MPNN, 4 pools × 2 seqs on 1,050 backbones | ~1–2 h | running |
| **AF2 screen, ~8,400 candidates** | **~23 h** | measured 10 s/design post-compile |
| Boltz 3-state ×5 on survivors | ~7 h | measured 39.7 s/prediction |
| ranking + submission audit | ~3 h | |
| **total** | **~39 h** | against 44.5 h |

**The AF2 screen is now the critical path and it has no slack.** Risks and the chosen
mitigations:

- *Partial completion is length-biased.* `predict_pdb.py` sorts by sequence length to
  minimise JAX compiles (21 distinct lengths × ~2 min each), which is the right trade while
  the run finishes, but means a run stopped early has seen only short designs. Accepted, with
  the consequence recorded here: if the screen is cut short, the surviving set is biased
  towards short binders and that must be stated rather than ranked over.
- *It is resumable.* Every design is checkpointed, so an interruption costs one design.
- *Fallback if it overruns:* drop to 1 sequence per pool (4,200 candidates, ~12 h). Not taken
  now because sequences-per-backbone is the main lever on hit rate at a 0.7% pass rate, and
  upstream `dl_binder_design` uses 8–48.

### Beta arm deferred behind MPNN

The 150-design beta arm is queued after MPNN rather than run concurrently: RFdiffusion and
AF2 coexisted happily earlier, but MPNN is short and serialising it removes any chance of the
two writing interleaved output into the same manifest.

---

## Stage 17 — MPNN at scale: vanilla weights are clearly better, correcting Stage 16's n=4 read

7,200 sequences, 900 backbones × 4 pools × 2 sequences, 7,140 s (2.0 h).

| pool | n | hydrophobic | charged | net |
|---|---|---|---|---|
| sol_neutral | 1800 | 0.311 | 0.590 | −8 |
| sol_basic | 1800 | 0.301 | 0.612 | −7 |
| **van_neutral** | 1800 | **0.413** | **0.457** | **−4** |
| **van_basic** | 1800 | **0.388** | **0.493** | **−3** |

**This reverses the earlier reading.** On n=4 I measured vanilla at charged 0.464 vs soluble
0.490 and concluded "most of the charge comes from the backbones and the Ala de-bias, not the
weight set". At n=1800 the gap is 0.457 vs 0.590 charged and 0.413 vs 0.311 hydrophobic —
substantial, not marginal. The n=4 sample was worthless and I should not have drawn a
directional conclusion from it.

Vanilla lands hydrophobic at 0.413, at or above the 0.30–0.40 target band, against soluble's
0.311. Given that composition tracked AF2 interface quality monotonically (Stage 14), the
vanilla pools are the better prospect — but both are screened, so this is now a measurement
rather than a prediction.

Net charge also halves (−4/−3 vs −8/−7), which matters for the electrostatics argument: the
discriminating feature is a buried Asp making the HLA surface more negative, and a less
negative binder is better placed to read it.

### Threading bug caught at launch

The first screen launch aborted immediately: *"7200 skipped: backbone not found"*. The
candidate-id parser stripped the pool prefix with `cid.split("_", 1)[0]`, which yields `sol`
for `sol_neutral_g_0_s1` — matching no pool, so every backbone lookup failed. Pool names
introduced in Stage 17 contain underscores; the pilot's names (`neutral`, `basic`) did not,
which is why it only surfaced now. Fixed to strip by full pool name, longest first. All 7,200
threaded on the retry.

The abort was loud and immediate rather than silent, which is the one good thing about it: a
parser that had half-worked would have screened a biased subset without saying so.

### Now running concurrently

- **AF2 screen**, 7,200 designs, ~20 h, resumable and checkpointed per design
- **Beta diversity arm**, 150 designs, ~4.5 h

RFdiffusion and AF2 coexisted without trouble earlier in the pilot, and the two write to
separate trees.

---

## Stage 18 — beta diversity arm delivers real topology spread (n=150)

150 backbones, ~4 h, 22 MB.

| | default `Complex_base` (n=900) | **beta `Complex_beta_ckpt`** (n=150) |
|---|---|---|
| helical fraction, median | 0.901 | **0.772** |
| extended fraction, median | 0.014 | **0.065** |
| extended >10% | 4 (**0%**) | 46 (**31%**) |
| extended >20% | 0 (0%) | 30 (**20%**) |
| mostly-helical (>80%) | 843 (**94%**) | 61 (**41%**) |

**This corrects the Stage 17 assessment, which called the shift "modest" on n=2.** At n=150 it
is substantial: a third of beta designs carry real extended content against none in the
default arm, and the mostly-helical share falls from 94% to 41%. The monoculture flagged
earlier is genuinely broken by this arm.

That is the third time today a sample under n=30 produced a directional conclusion that larger
n reversed — the MPNN weight-set comparison (n=4 said marginal, n=1800 said substantial), the
AF2 pass rate (n=32 said 3.1%, n=272 said 0.7%), and now this. The pattern is consistent
enough to be a rule: **do not state a direction from n<30 in this project.** Each of those
small samples was cheap to take and each cost a retraction.

What this does *not* establish: that beta-topology designs make better interfaces. Topology
diversity is a submission-diversity objective (`docs/PIPELINE.md §5`), not an affinity
argument. Whether these designs survive the AF2 screen at a better or worse rate than the
helical arm is a separate question, and the screen will answer it — the beta designs are
labelled and screened as their own cohort so the comparison is measurable.

### Screening order

The 7,200-design helical screen is already mid-flight and `predict_pdb.py` enumerated its
input directory at startup, so the 1,200 beta designs are screened as a **second, sequential
run** rather than being dropped into the live directory. Adding files underneath a running
screen would either be ignored or force a restart that discards the current compile state.

Projected: helical screen ends ~Sun 19:40 Melbourne, beta screen +3.3 h to ~23:00, Boltz on
survivors ~7 h to Mon ~06:00, ranking and audit ~3 h to ~09:00, against an 18:00 freeze.

---

## Stage 19 — interim screen (2,000/7,200): the composition hypothesis does not survive

Screened in length-sorted order, so these are the shortest designs (65–70 aa).

| pool | n | pae median | pae best | pae<10 | plddt median | binder_rmsd median |
|---|---|---|---|---|---|---|
| sol_neutral | 500 | 26.74 | **5.05** | **2** | 90.1 | 1.54 |
| sol_basic | 500 | 26.80 | 12.05 | 0 | 90.2 | 1.58 |
| van_neutral | 500 | 26.83 | 8.02 | **2** | 88.2 | 1.66 |
| van_basic | 500 | 26.83 | 15.82 | 0 | 89.1 | 1.74 |

**Vanilla's better composition has not translated into better interfaces.** The pae medians are
identical to two decimal places across all four pools (26.74–26.83), and both full-criteria
passes so far are `sol_neutral`. This is a clean correlation-versus-causation failure and it
is mine: composition *correlated* with pae in the pilot (best quartile 0.324 hydrophobic vs
worst 0.256), so I used it to justify adding vanilla pools — but *intervening* on composition
moved composition (0.413 vs 0.311 hydrophobic, net −4 vs −8) and did not move pae at all.
A correlate of quality is not a lever on quality.

The decision to run both pools rather than switch to vanilla outright is what makes this
visible instead of invisible. Had I switched, the pipeline would have shipped vanilla-only
designs on a rationale the data does not support.

**The `basic` pools have zero passes in 1,000 designs** against 4 in 1,000 for the neutral
pools. Numbers far too small to conclude from, but the positional R/K bias is a HYPOTHESIS
(so labelled in the config since the start) and nothing yet supports it. Watch, do not act.

Binder length does not shift the median: 26.75–26.87 flat across 65–70 aa, with hits appearing
at 67/69/70 as rare tail events rather than a trend. So the length-sorted ordering is less of a
bias on the *rate* than feared.

### Yield risk is now the live problem

Interim full-criteria rate is **2/2000 = 0.10%**, against 0.7% in the pilot. The pilot figure
came from 272 designs on 34 backbones and was itself revised down from 3.1%; the honest range
for the final 7,200 is wide, roughly **7–50 survivors**, most plausibly 15–25.

**That may not fill 20 submission slots on the hard gate alone.** The fallback is already
validated and does not need inventing:

1. `pae_interaction < 10` — primary, published, tail-validated against scrambles (0/32)
2. then rank by `binder_aligned_rmsd` and `plddt_binder`, which separated real designs from
   composition-matched scrambles by **+21.3 Å** and **−19.1 points** respectively (Stage 15)
3. then the Boltz three-state differential at 5 samples

Designs promoted via step 2 rather than step 1 must be labelled as such in the submission
manifest and the methods note. They clear a validated fold-self-consistency bar but **not** the
published interface criterion, and saying otherwise would misrepresent them.
