# First prompt for the Claude Code agent

Paste as the opening message in the repo root.

---

You are the design agent for Adaptyv × Anthropic 2026 **Challenge 05 — KRAS G12D pMHC**.

Read these three files first, in order, and treat them as the operating contract:
- `CLAUDE.md` — rules, guardrails, manifest contract
- `docs/CHALLENGE_SPEC.md` — the verified official spec (fetched 2026-11-06)
- `docs/PIPELINE.md` — the 71-hour staged plan, stack decisions, gates
- `configs/challenge.yaml` — target sequences, hotspots, thresholds

**The objective:** de novo binders that bind HLA-A\*11:01 + β2m + VVVGADGVGK, and do **not**
bind the wild-type complex (VVVGAGGVGK) or the peptide-free MHC. Ranking is lexicographic —
mutation specificity first, peptide dependency second, affinity third. A tight framework
binder loses to a modest clean discriminator. Design for the differential.

**Hard deadline:** Mon 9 Nov 22:59 Melbourne. Internal freeze Mon 18:00. One A100 80 GB.

### Do now — Stage 0 and Stage 1 only. Stop at Gate 1.

1. **Environment audit.** Run the checks in `CLAUDE.md §3`. Report what is actually installed,
   which weights are on disk, GPU state, free scratch, and whether the Boltz-2 MSA server at
   10.148.0.19:8080 responds. Do not assume anything is present.

2. **Assess reuse.** We have three prior pipelines for this same competition:
   `tnf-a-comp2`, `glp1r-comp3`, `tfr1-comp4`. Tell me concretely what ports directly
   (target prep, the RFdiffusion Docker wrapper, the two-pool MPNN pattern, the Boltz-2 filter
   with its negative-control calibration, the submission verification harness) and what has to
   be new. The three-state differential is the main new piece. Do not reimplement anything that
   already works.

3. **Target prep (Stage 0).** Build the three systems — `g12d`, `wt`, `empty` — from 7OW6 with
   the TCR removed. Note that 7OW4's central bulge is unresolved, so 7OW6 is the only usable
   source for the peptide centre; a bulge-up state is a *modelled hypothesis* and gets dropped
   if it costs more than an hour. Place hotspots on peptide chain C at p4, p5, p7, p8 — and
   assert in code that they reference the target chain's renumbered residues after cropping.
   That exact bug cost us a round in comp3.

4. **Pilot (Stage 1).** 32 designs end to end: generate → MPNN → three-state Boltz-2 →
   metrics → manifest rows. Then report against the Gate 1 criteria in `docs/PIPELINE.md`,
   plus measured throughput so we can check the Stage 3 co-folding matrix actually fits.

**Stop at Gate 1 and report. Do not start bulk generation without my go-ahead.**

### Report format
What ran · numbers · what failed · remaining uncertainty · recommended next action. Nothing else.

### Two things I want you to push back on if you disagree
- Rejecting RFD4-Proteina as the generator spine (reasoning in `docs/PIPELINE.md §0`).
- Gating the RFantibody nanobody track behind Gate 2 at H+14 rather than starting it now.

If the environment audit changes the picture on either, say so before building.