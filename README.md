# kras-g12d-pmhc-comp5

De novo binder design against **HLA-A\*11:01 · β2-microglobulin · VVVGADGVGK** (KRAS G12D,
peptide position 6). Anthropic × Adaptyv 2026, Challenge 05, **Track 3 — submitted, 20 designs**.

The challenge asks for three things, in this order: bind the G12D complex, **do not** bind the
wild-type complex (`VVVGAGGVGK`), and **do not** bind the same HLA with an empty groove.
Ranking is lexicographic, so a modest clean discriminator beats a tight pan-HLA binder.

---

## Result

20 single-chain de novo minibinders, 69–95 aa, all inside the published 40–100 aa molecule
category, over **12 distinct backbones** from four independent generation arms.

| objective | metric | result |
|---|---|---|
| **1. mutation specificity** | `iptm_margin_wt` (G12D − WT) | all 20 positive, median +0.010, **but only 2/20 clear the 0.053 noise floor** |
| **2. peptide dependency** | `iptm_margin_empty`, buried peptide area | all 20 positive, 10/20 above noise; **71.5–199.1 Å² buried peptide** vs **0 Å²** for framework decoys |
| **3. affinity proxy** | interface size, AF2 confidence | **660–1,322 Å²** total interface (median 989); AF2 `pae_interaction` **3.84–9.70**, median 5.97 |

**We are strongest on objective 2, credible on 3, and objective 1 is a long shot for most of
the set.** That is mechanistic, not a pipeline failure: the mutant Asp at p6 is buried in the
F-pocket, salt-bridged to HLA Arg114 (2.95 / 3.04 Å) and H-bonded to Gln70, so no binder reads
it directly. Discrimination can only be indirect.

## Pipeline

```
RFdiffusion  ->  ProteinMPNN / SolubleMPNN  ->  AF2 initial-guess  ->  Boltz-2 three-state
 1,250 backbones      11,152 sequences          158 pass (1.4%)        g12d / wt / empty
```

RFdiffusion was the **only** generative model used — standard complex checkpoint plus
`Complex_beta_ckpt` for topology diversity. FreeBindCraft, RFantibody and RFD4-Proteina were
evaluated and not used.

**AF2 initial-guess is the selection filter**, not Boltz. It applies the published criteria
(`pae_interaction < 10`, `plddt_binder > 80`, `target_aligned_rmsd < 5 Å`) and was calibrated
against 32 composition-matched scrambled sequences, **none of which reached pae < 10**.
Boltz-2 supplies the three-state differentials that AF2 structurally cannot, since it never
sees the wild-type peptide or the empty groove.

## Three things measurement changed

1. **Deep resampling of proven backbones hit 10.4% against 0.39% for fresh geometry** — ~27×,
   at a fifth of the GPU cost. Once a backbone supports a solution, more sequence draws on it
   beat more geometry. This was the single most effective intervention of the campaign.
2. **Single-sample ipTM margins are noise.** Between-seed Spearman ρ for `iptm_margin_wt` is
   **+0.239** at one diffusion sample and **+0.846** at five. Averaging costs 2.41×, not 5×,
   because Boltz shares the trunk. Every margin here is a 5-sample average.
3. **Composition is a marker, not a lever.** Hits have better composition than non-hits
   (0.372 vs 0.341 hydrophobic at n=7,200), but deliberately engineering composition — vanilla
   vs soluble ProteinMPNN weights — moved composition substantially and the hit rate not at all
   (z = +0.76). Selecting on a correlate works; intervening on it did not.

## Layout

| path | contents |
|---|---|
| `configs/challenge.yaml` | target sequences, hotspots, thresholds — single source of truth |
| `configs/thresholds_calibrated.yaml` | thresholds derived from negative controls, not guessed |
| `scripts/` | numbered stages, each independently runnable and resumable |
| `structures/` | the three prepared assay states built from 7OW6 |
| `results/RUNLOG.md` | **what actually ran, in order, including every failure and retraction** |
| `results/processed/` | manifests, AF2 screens, calibration report |
| `results/rejected/` | every rejected candidate, with a reason |
| `submission/` | final CSV, pre-submission audit, methods note |

`results/RUNLOG.md` is the honest record and the most useful file here — it includes the
hypotheses that failed, a self-inflicted 30-minute GPU stall, a 6 GB disk leak from
RFdiffusion trajectory dumps, and three occasions where a sample under n=30 produced a
conclusion that larger n reversed.

## Reproducing

```bash
source env.sh                                   # or activate the boltz/mpnn/af2ig conda envs
python scripts/verify_targets.py                # sequence + hotspot assertions, must pass
python scripts/01_prep_targets.py               # build g12d / wt / empty from 7OW6
python scripts/02_crop_target.py                # alpha1/alpha2 platform crop for generation
bash   scripts/30_generate.sh 900 60            # RFdiffusion backbones
python scripts/04_mpnn.py --rfdiff gen/rfdiff   # four labelled sequence pools
bash   scripts/41_af2_screen.sh                 # AF2 initial-guess, the primary filter
bash   scripts/80_finalise.sh                   # Boltz three-state, ranking, submission audit
```

Nothing in this repo submits. `scripts/70_submission.py` writes the CSV and audits it;
submission is manual.

## Caveats

Everything here is **predicted**, never measured. No K_D is claimed, and Boltz-2's affinity
head was not used. PyRosetta, PRODIGY, ThermoMPNN and APBS were unavailable, so
`rosetta_ddg`, `prodigy`, `thermompnn_score` and `electrostatic_delta` are recorded as
`not_run` rather than estimated. There is no Rosetta FastRelax between sequence design and
AF2, which upstream `dl_binder_design` performs — a real deviation.

All designs are de novo. No known binder sequence entered the pipeline at any stage; novelty
was verified after the fact against the JDI TCR lineage in 7OW6 (maximum 20-residue window
identity **0.40**, threshold 0.70).
