# comp5 pipeline — 71-hour plan

Kickoff: **Fri 6 Nov 2026, ~24:00 Melbourne**. Freeze: **Mon 9 Nov 18:00**. Deadline Mon 9 Nov 22:59.
One A100 80 GB. H+N below is hours from kickoff.

The whole design is shaped by one fact: **objectives 1 and 2 are differentials, and they rank
above affinity.** A framework binder is worthless here no matter how tight. So the pipeline
spends its compute on discrimination, not on affinity maturation.

---

## 0. Stack decision

### Keep (proven in comp2/3/4, zero setup risk)
| Component | Role |
|---|---|
| **RFdiffusion** (Docker, from `tnf-a-comp2/rfdiffusion.sh`) | Track A1 backbone generation, hotspot-conditioned |
| **ProteinMPNN / SolubleMPNN** | sequence design, two pools |
| **Boltz-2** + self-hosted MSA (10.148.0.19:8080) | primary co-folder, three-state |
| **PyRosetta FastRelax + InterfaceAnalyzer, PRODIGY** | interface energetics, contact composition |
| **ANARCI, IgFold/ABodyBuilder2** | antibody-class annotation and modelling |
| **comp2 `04b/04c` negative-control + cutoff calibration** | threshold setting |
| **comp3 submission-verification harness** | final CSV audit |

### Add (high value, low setup risk)
| Component | Role | Why |
|---|---|---|
| **FreeBindCraft** | Track A2 generator | PyRosetta-free BindCraft; end-to-end AF2 hallucination with hotspots, self-filtering. Different failure mode from RFdiffusion — genuine ensemble diversity, not just more seeds. |
| **Chai-1** (`chai-lab`, Apache-2.0) | second co-folder | Differentials must be multi-model. Chai-1 is the cheapest orthogonal AF3-class check. |
| **OpenDDE** (`opendde_abag.pt`) | antibody-track co-folder | Already in the Actimo stack. Co-folders systematically underperform on Ab–Ag interfaces; the abag checkpoint is the mitigation. |

### Gated (only if Track A clears its gate)
| Component | Role |
|---|---|
| **RFantibody** (MIT) | Track B — nanobody CDR design on fixed VHH framework. ANARCI-clean by construction, hits a separate molecule category. |

### Rejected for this run
| Component | Why not |
|---|---|
| **RFD4-Proteina** | Released ~2026-10-19, under 3 weeks old. Design quality not independently validated, and there are reports of pathological charged/repetitive output on unconditional generation. No pMHC track record. In a 71-hour window there is no budget to debug a new generator — the cost of it quietly producing junk is a wasted day we cannot recover. **Optional diversity batch only, if and only if Track A finishes early and the GPU would otherwise idle.** |
| **Chai-2** | Closed-source, partner-only. Unavailable. |
| **Mol-1, OpenFold3, Genie3, BoltzGen** | Each is defensible but none displaces an incumbent. Two co-folders plus OpenDDE is already enough consensus for the time available; a third generator adds integration risk without adding discrimination power. |
| **Engineered TCR format** | Not de novo, not in the allowed format list. |
| **scFv / Fab** | Highest literature ceiling for pMHC specificity, but multi-chain, disulfide-dependent, and the expression system is undisclosed. Too much risk at 71 h. |

**Answer on RFD4: not as the spine.** It is the most interesting model on the list and the
wrong one for this week.

---

## 1. Stage plan

### Stage 0 — environment + target prep (H+0 → H+4)
1. Environment audit per `CLAUDE.md §3`. Report actual GPU, weights, MSA server health.
2. Fetch **7OW6**, strip TCR chains, keep A (heavy) / B (β2m) / C (peptide). This is the **primary target state** — the only structure with a resolved peptide centre.
3. Build three systems at that state:
   - **S1 `g12d`** — 7OW6 A/B/C as-is.
   - **S2 `wt`** — in-silico D6→G on chain C, local relax only (do not re-pack the whole groove).
   - **S3 `empty`** — chain C deleted, light relax of groove.
4. Build a **secondary "bulge-up" state** for S1/S2 by modelling, since 7OW4's centre is unresolved. **Mark this explicitly as a modelled hypothesis, lower confidence.** If it costs more than 1 h, drop it and run single-state.
5. Sanity: NetMHCpan-4.1 on both peptides against A\*11:01. Confirm both plausibly load — the WT counter-target must actually exist as a complex. Record the %Rank for each; if WT scores as a non-binder, note the conflict rather than resolving it silently.
6. Define hotspots on **chain C**: `C4, C5, C7, C8` (exposed bulge), with `C6` tracked in contact accounting but **not** a primary hotspot (it is buried in this state). Explicitly exclude α1/α2 helix residues.
   - **Gotcha from comp3 round 2:** hotspots must name the *target* chain's renumbered residues after cropping, not the binder's. Assert this in code, don't eyeball it.
7. Crop a 15 Å shell around hotspots for RFdiffusion (reuse `02_crop_target`).

### Stage 1 — pilot smoke test (H+4 → H+7) — **GATE 1**
End-to-end on **32 designs**: generate → MPNN → three-state Boltz-2 → metrics → manifest row.

**Gate 1 pass criteria:**
- Pipeline runs end to end with no manual intervention.
- At least some designs show `peptide_contact_fraction` ≥ 0.25.
- The three-state differential *discriminates at all* — i.e. the ipTM spread between `g12d` and `empty` is wider for high-peptide-contact designs than for framework-dominated ones.
- Measured throughput lets you project the Stage 3 matrix inside its block.

**If Gate 1 fails:** stop, report, do not proceed to bulk. Diagnose hotspot placement first — it is the most likely cause.

### Stage 2 — bulk generation + sequence design (H+7 → H+26)
- **Track A1 — RFdiffusion**, hotspot-conditioned over the bulge. Over-generate: target **3,000–4,000 backbones**. Expect only a small fraction to make extensive peptide contact — that is the known behaviour of this protocol, not a bug.
- **Track A2 — FreeBindCraft**, same hotspots, run in parallel where VRAM allows. Target **600–1,000 designs**.
- Length mix: majority in **40–100 aa** (the one published molecule category), plus a smaller block outside it for category hedging.
- **Sequence design:** ProteinMPNN + **SolubleMPNN** (expression-safe), two pools:
  - *neutral pool* — default.
  - *basic-biased pool* — mild enrichment for Arg/Lys in positions projecting over p6. **This is a hypothesis** derived from the R50α/K70α electrostatic mechanism in the 2022 structure, not a validated design rule. Keep the pools separately labelled so we can tell afterwards whether it did anything.
- **Cheap pre-filter before co-folding:** ProteinMPNN peptide-recovery test — thread `VVVGADGVGK` vs `VVVGAGGVGK` through each designed interface and keep designs scoring the mutant more favourably. This is the single highest-leverage filter per GPU-second; it runs before the expensive stage and is the Baker-protocol specificity proxy.

**GATE 2 at H+14:** if Track A1 + A2 are both producing designs that clear Gate 1 criteria,
**start Track B (RFantibody nanobodies)** in whatever GPU headroom exists. If not, Track B is
cancelled — do not start it late.

### Stage 3 — three-state co-folding matrix (H+26 → H+48) — the expensive stage
This is the budget bottleneck: *states × designs × models*. Manage it by triage, not by hope.

1. **Triage pass:** Boltz-2, `g12d` state only, on all surviving designs. Rank by interface confidence + peptide contact. Cut to **≤600 designs**.
2. **Full matrix** on survivors: **3 states × 2 models** (Boltz-2 + Chai-1). Track B designs additionally through **OpenDDE `opendde_abag.pt`**.
3. Cache the shared MHC MSA once — do not re-derive it per design.
4. Checkpoint after every batch. This stage will be interrupted; assume it.

**If projected runtime exceeds the block:** cut the design count, not the number of states or
models. Three states across two models is the product; one model or two states makes the
differential meaningless.

### Stage 4 — metrics and the three objective proxies (H+48 → H+56)

| Objective | Proxy | Notes |
|---|---|---|
| **1. Mutation specificity** (ranks first) | `iptm_margin_wt` = score(g12d) − score(wt), required in **both** models; plus Rosetta ddG differential; plus `mpnn_peptide_recovery_gap` | The hardest and least reliable proxy. No current method predicts the WT/mutant K_D window for a de novo binder. Treat as enrichment, not prediction. |
| **2. Peptide dependency** | `iptm_margin_empty` = score(g12d) − score(empty), both models; **and** `peptide_contact_fraction` ≥ threshold over p4–p8 | Hard gate. Baker's crystallised mart1-3 buried ~1,200 Å² with ~31% from peptide — use that as the sanity scale, not as a cutoff. |
| **3. Affinity** | Rosetta ddG + PRODIGY on `g12d`, interface confidence | **Rank signal only. Never report as K_D.** Boltz-2's affinity head is weak for protein–protein. |

Tiebreak / hypothesis metric: **electrostatic complementarity** — APBS surface potential of the
paratope over p6, and the g12d-minus-wt potential difference. Flagged as hypothesis, not a
validated metric for de novo binders.

**Threshold calibration (do not invent cutoffs):** reuse the comp2 `04b/04c` pattern —
shuffled-sequence negatives plus a deliberately-constructed framework-binder decoy placed on
the α1/α2 helix. Set the "no detectable binding" thresholds for `wt` and `empty` from where
those negatives land, mirroring the organisers' 300%-over-negative-control rule. Every
threshold in `configs/challenge.yaml` marked `provisional` must be replaced by a calibrated
value before ranking, or the ranking is meaningless.

### Stage 5 — ranking, diversity, selection (H+56 → H+64)
Hard gates first, in this order (anything failing is moved to `results/rejected/` with a reason):
1. Length 10–250; molecule class valid; ANARCI-clean for Track B.
2. `peptide_contact_fraction` above calibrated threshold.
3. `iptm_margin_empty` positive in both models.
4. `iptm_margin_wt` above calibrated threshold in both models.
5. Developability: SolubleMPNN pass, no aggregation-prone or liability motifs, ThermoMPNN stability.
6. **De novo check:** BLAST/identity screen against known TCRs and TCR-mimics. Hard fail on substantial similarity.

Then rank in the **competition's own lexicographic order** — specificity margin, then peptide-dependency margin, then affinity proxy. Do not build a weighted composite score; it would contradict how designs are actually ranked.

Diversity: TM-score clustering, cap per-scaffold representation, keep both tracks and multiple
docking geometries represented. Avoid the common failure where all submissions converge on one
tool and one binding mode.

### Stage 6 — submission (H+64 → H+66, freeze H+66 = Mon 18:00)
1. **Confirm the track with Harjas** (prior rounds Track 3 → cap 20). Do not guess the cap.
2. **Re-check the live submit form for the exact CSV columns** before writing the file. Our prior rounds used `name,sequence,molecule_class`; verify, do not assume.
3. Run the comp3 verification harness: length bounds, duplicate check, invalid residues, de-novo run-length/identity check, column names, row count ≤ cap, rank order.
4. Write the methods note: **state that we target the peptide bulge, and link the predicted complexes.** The page explicitly rewards declaring the surface — *"Tell us which surface you target and link your evidence. If you don't, we take it from a predicted complex."*
5. Produce the pre-submission audit for Harjas's approval. **Do not submit.**

Buffer H+66 → H+71 is for failure, not for extra designs.

---

## 2. Compute budget

| Stage | Block | GPU-h (est.) |
|---|---|---|
| 0 Target prep | H+0–4 | ~0.5 |
| 1 Pilot | H+4–7 | ~2 |
| 2 Generation + MPNN | H+7–26 | ~16 |
| 3 Co-folding matrix | H+26–48 | ~20 (the bottleneck) |
| 4 Metrics | H+48–56 | ~4 |
| 5 Ranking | H+56–64 | ~1 |
| 6 Submission | H+64–66 | 0 |

Estimates are order-of-magnitude and must be replaced by measured throughput from the pilot.
If the pilot says the Stage 3 matrix does not fit, cut the design count at Stage 3 step 1
immediately — not at H+45.

---

## 3. Known failure modes

- **Mutation specificity is the ceiling.** Discrimination is electrostatic and dynamic because the Asp is buried. No method reliably predicts the WT-vs-G12D window for a de novo binder. The differential filters raise the odds; they do not guarantee anything. This is the most likely place to fail in the wet lab.
- **Absolute K_D is not predictable.** We rank; we do not promise affinities.
- **"No binding" calls are model-dependent.** AF3-class models over-call non-binders as binders — hence two-model consensus plus explicit negative controls.
- **Peptide-free HLA-A\*11:01 is physically unstable.** Our `empty` model is an approximation of a target that may itself behave oddly in the assay.
- **Induced fit (bulge up ↔ down) is poorly captured** by static co-folding. The bulge-up state is modelled, not observed.
- **Expression risk** is unquantified — the expression system is undisclosed. This is why single-chain minibinders and nanobodies are preferred over scFv/Fab.