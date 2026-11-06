# CLAUDE.md — comp5 / KRAS G12D pMHC

Adaptyv × Anthropic 2026, **Challenge 05**. Design de novo binders that bind
HLA-A\*11:01 + β2m + **VVVGADGVGK** (KRAS G12D), and do **not** bind the wild-type
complex (VVVGAGGVGK) or the peptide-free MHC.

Read `docs/CHALLENGE_SPEC.md` and `docs/PIPELINE.md` before doing anything.

---

## 1. Hard deadline and time budget

- Official close: **Sunday 8 Nov 2026, 23:59 AoE** = **Monday 9 Nov 2026, 22:59 Australia/Melbourne**.
- **Internal freeze: Monday 9 Nov, 18:00 Melbourne.** Submission CSV done, audited, nothing running.
- Total budget from kickoff: **~71 h wall clock, one A100 80 GB.**

Time is the binding constraint, not ideas. Rules:

- Never start a stage that cannot finish inside its block in `docs/PIPELINE.md`. Cut scope instead.
- Every long run is checkpointed and resumable. Assume the session can drop at any point.
- Every stage writes a timestamped row to `results/RUNLOG.md`: stage, start, end, GPU-h, inputs, outputs, exit status.
- If a stage overruns its block by >25%, stop it, record what you have, and move to the next stage. Report the cut.
- No new model, framework or dependency gets introduced after **H+20** unless it is the only way to unblock a failed stage.

## 2. Non-negotiables

1. **De novo only.** Never seed, thread, graft or mutate from a known binder — no JDI/JDIa41b1/JDIa96b35 TCRs, no published TCR-mimic antibodies, no scDb/CAR sequences. These are *reference and control material only*. Any design whose sequence is substantially similar to a known binder is disqualified.
2. **Never fabricate.** No invented K_D, ipTM, ddG, PDB ID, benchmark number, or novelty claim. If a metric was not computed, write `not_run`. If a run failed, write `failed` and the reason. Never fill a manifest column with a plausible-looking number.
3. **Predictions are not measurements.** Nothing in this repo is a validated binder. Language in reports: "predicted", "modelled", "ranked". Never "binds", "affinity of X".
4. **Never overwrite raw outputs.** `results/raw/` is append-only. Processed metrics go to `results/processed/`. Rejected candidates go to `results/rejected/` with a reason column — never deleted.
5. **No external actions without explicit approval from Harjas:** no `git push`, no commits to remote, no uploads, no submission to the ProteinBase form, no posting anywhere.
6. **Verify the spec, don't remember it.** If you need a rule, read `docs/CHALLENGE_SPEC.md`. If the answer isn't there, fetch the live challenge page. Do not infer a rule.

## 3. Environment — check before you assume

First action of the session, before any design work:

```
nvidia-smi                      # confirm A100 80GB, free VRAM, no other tenant
git status && git log --oneline -10
conda env list / pip list       # what is actually installed
ls -la weights/ checkpoints/    # which model weights are already on disk
df -h                           # scratch space for structures + co-folds
curl -s 10.148.0.19:8080/health # self-hosted MSA server for Boltz-2
```

Do not assume a model, checkpoint or GPU is available. Report what is actually present, then
propose the smallest set of installs needed. Reuse what is already working from
`tnf-a-comp2`, `glp1r-comp3`, `tfr1-comp4` rather than reimplementing.

## 4. Repo layout

```
configs/challenge.yaml       # target sequences, hotspots, thresholds — single source of truth
configs/                     # hydra configs per stage
scripts/                     # 01_prep … 09_submit, numbered, each independently runnable
src/comp5/                   # importable modules, tested
structures/                  # prepared target states (G12D / WT / empty, bulge-down / bulge-up)
results/raw/                 # append-only model outputs
results/processed/           # metrics tables, manifests
results/rejected/            # rejected candidates + reason
results/RUNLOG.md            # what actually ran
docs/                        # spec, pipeline, decisions, methods
submission/                  # final CSV + audit
```

## 5. Manifest contract

Every candidate carries, in `results/processed/manifest.parquet` (+ CSV mirror):

`candidate_id, sequence, length, molecule_class, track, design_method, model_version,
params_hash, seed, target_state, input_structure, peptide_contact_fraction,
n_peptide_contacts, iptm_g12d, iptm_wt, iptm_empty, iptm_margin_wt, iptm_margin_empty,
pae_interface_g12d, rosetta_ddg_g12d, rosetta_ddg_wt, prodigy_g12d, mpnn_peptide_recovery_gap,
electrostatic_delta, thermompnn_score, anarci_class, qc_status, reject_reason,
diversity_cluster, rank, artifact_path`

Unknown = `not_run`. Never null-filled with a guess.

## 6. Scientific guardrails specific to this target

- **The mutant Asp (peptide p6) is buried in the F-pocket in the TCR-bound state** (7OW6; Poole & Karuppiah, *Nat Commun* 2022). You cannot design a binder that directly reads the Asp side chain. Specificity is **indirect** — altered surface electrostatics and bulge conformation. Any design rationale that claims "contacts Asp12 directly" is wrong unless the model actually shows it; check, don't assert.
- **Framework binders fail 2 of 3 objectives.** A binder on the α1/α2 helices or β2m binds all three targets and ranks last. `peptide_contact_fraction` is the primary hard gate, not a nice-to-have.
- **Co-folders are unreliable here.** Boltz-2 drops to ~0.70 DockQ on unseen pMHC; AF3-class models over-call non-binders as binders. Therefore: rank on **differentials across ≥2 models** (G12D minus WT, G12D minus empty), never on an absolute score from one model.
- **Boltz-2's affinity head is weak for protein–protein.** Use it for interface geometry and confidence only. Do not report its affinity output as a K_D.

## 7. Working style

- Report what actually ran, the numbers, what failed, remaining uncertainty, next action. In that order.
- Flag every assumption you make that is not in `docs/CHALLENGE_SPEC.md`.
- Prefer a small validated pilot over a large speculative run. Always smoke-test end-to-end before bulk.
- Concise. No restating the task back. No generic background.
- If a decision materially changes the plan (dropping a track, changing a target state, missing a gate), stop and ask Harjas rather than deciding alone.