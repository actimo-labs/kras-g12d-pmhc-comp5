# Challenge 05 — KRAS G12D pMHC: verified specification

Source: https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/pmhc
Fetched and verified: **2026-11-06**. Status at fetch: **Submissions open**.

Everything in §1–§5 is official text from the challenge page. §6 is marked unknown.
§7 is scientific background, not rules. Do not blur these.

---

## 1. Target construct (OFFICIAL)

| Field | Value |
|---|---|
| Peptide | `VVVGADGVGK` — KRAS P01116 residues 7–16, 10 aa |
| Mutation | KRAS Gly12 → Asp, **peptide position 6** |
| Allele | HLA-A\*11:01, four-field `A*11:01:01:01`, IPD-IMGT/HLA `HLA00043` |
| Heavy chain | UniProt Q5S3G3, residues 25–305, 281 aa |
| β2-microglobulin | UniProt P61769, residues 21–119, 99 aa |
| Assay construct | Complex of the three chains · **tags not shown** |
| Structure reference | **PDB 7OW6, chains A–C, TCR removed** |
| Recommended epitope | "The peptide and the groove around it" |

**HLA-A\*11:01 heavy chain**

```
GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPT
```

**β2-microglobulin**

```
IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM
```

Organisers' note: *"The rest of the complex is the same for both peptides, and the same again
with no peptide, so a binder that sits there cannot tell them apart."*

## 2. Counter-targets (OFFICIAL)

1. **Wild-type complex** — same HLA-A\*11:01 + β2m carrying `VVVGAGGVGK` (Gly at p6). Measured in the same conditions.
2. **Peptide-free MHC** — HLA-A\*11:01 + β2m, empty groove. Measured in the same conditions.

No other alleles or irrelevant-peptide controls are scored.

## 3. Assay and the binding call (OFFICIAL)

All three targets measured in the same conditions, **design immobilised, target in solution**.

> "A design counts as binding a target when its sensorgram gives a fittable curve and a K_D,
> or, where the curve cannot be fitted, when the association shift clears 300% over the
> negative control. Otherwise we call it no detectable binding."

Instrument is not named on the page (BLI vs SPR → unknown). The immobilised-design /
solution-target geometry is explicit.

## 4. Ranking (OFFICIAL — lexicographic, not weighted)

1. **Mutation specificity** — measurable K_D on G12D **and** no detectable binding on WT. Designs that bind both rank after these, ordered by K_D(WT)/K_D(G12D), largest first.
2. **Peptide-dependent binding** — measurable K_D on G12D **and** no detectable binding on peptide-free MHC. Designs that bind it rank after these, by K_D(empty)/K_D(G12D), largest first. *"A binder that reads the HLA or the β2-microglobulin instead of the peptide binds all three targets and ranks there."*
3. **Affinity** — lowest K_D on G12D. *"A design that binds the peptide and the groove around it has an advantage. Tell us which surface you target and link your evidence. If you don't, we take it from a predicted complex."*

Winner must clear all three steps. Each step also awards an honourable mention.

**Implication:** affinity is the *third* tiebreak. A 100 nM clean discriminator beats a 1 nM
pan-HLA binder. Optimise the differential first.

## 5. Submission requirements (OFFICIAL)

| Field | Value |
|---|---|
| Length | **10–250 aa** |
| Format | Single chain protein, nanobody, scFv, or Fab |
| Designs per participant | **Up to 40 (Track 1), up to 20 (Tracks 2 and 3)** |
| Design | **De novo designs only** |
| Deadline | **Sunday 8 Nov 2026, 23:59 AoE** (= Mon 9 Nov 22:59 Australia/Melbourne) |
| Window | Oct 26 – Nov 8, 2026 |

FAQ detail:
- De novo proteins (minibinders, larger proteins, microbinders), nanobodies and antibodies allowed.
- Single-chain proteins: min 10 aa, max 250 aa.
- *"we consider a design to be a nanobody/antibody if it is annotated as such by [ANARCI]"*. Advice: use common framework regions as a starting point.
- scFvs may be reformatted to the organisers' preferred linker.
- Performance is **stratified into five molecule categories**. Only category 1 is published: **Protein minibinders, 40–100 aa inclusive.** Categories 2–5 are not listed on the page.

Post-deadline: Anthropic + Adaptyv select designs for synthesis using a Claude-based workflow
(published after the competition), then Adaptyv wet-lab validation.

## 6. UNKNOWN — do not invent

- Instrument (BLI vs SPR) and buffer composition.
- Tag / biotinylation strategy ("tags not shown").
- Whether the pMHC is refolded monomer or single-chain trimer.
- Expression system for synthesised designs (E. coli vs mammalian) — **affects disulfide-dependent scFv/Fab viability**.
- Molecule categories 2–5.
- Submission CSV column names for this round. Our prior rounds used `name,sequence,molecule_class`; **re-check the live submit form before freezing the CSV.**
- Which track Harjas is entering (prior rounds: Track 3 → cap 20). **Confirm before selecting the final set.**

## 7. Scientific background (NOT rules — references and controls only)

Key structural fact, from Poole & Karuppiah et al., *Nat Commun* 13:5333 (2022),
DOI 10.1038/s41467-022-32811-1; PDB **7OW4** (pMHC alone, 1.81 Å) and **7OW6** (TCR–pMHC):

- In the TCR-bound state the mutant **Asp at peptide p6 is buried in the HLA groove**, salt-bridging R114 in the F-pocket (stabilised by D116) plus an H-bond to Q70. The shorter WT Gly cannot make these.
- **Induced fit:** in free pMHC (7OW4) the central peptide region (residues 4–6) faces away from the HLA and is **disordered — no observable electron density**. On TCR binding it flips toward the F-pocket.
- Specificity is therefore **indirect**: the buried Asp shifts HLA surface electrostatics negative, which the TCR reads via R50α and K70α.
- Measured affinities for calibration: parent JDI TCR K_D 63 µM on G12D, no measurable WT binding; affinity-enhanced JDIa41b1 743 pM (G12D) vs 3.0 µM (WT) ≈ 4,000×; JDIa96b35 64 pM vs 425 nM ≈ 6,000×.
- Alanine scan: peptide D6, G4, A5, G9 most critical for TCR binding.

**Consequence for 7OW4:** its central bulge is unresolved, so it cannot be used directly as a
design target for the peptide centre. **7OW6 (TCR removed) is the only structure with a
resolved peptide centre and is what the organisers point at.** A "bulge-up" alternative state
must be *modelled*, and is a lower-confidence hypothesis.

Prior art — **reference and control only, never a design seed:**
- Rosenberg/NCI HLA-C\*08:02-restricted G12D TCRs (different allele); murine A\*11:01 G12D TCR (NCT03745326).
- Immunocore affinity-enhanced A\*11:01 G12D TCR → IMC-KRAS-G12D (source of 7OW4/7OW6).
- TCR-mimic antibodies vs KRAS-G12D/HLA-A\*11 (phage-derived, TCRm-ADC).
- Baker lab de novo pMHC-I binders, *Science* 389:386–391 (2025), DOI 10.1126/science.adv0185 — the only published de novo (non-antibody) pMHC-I campaign; 11 targets, peptide-focused RFdiffusion, nM affinities.
- Chai-2 reported 2/50 hits on HLA-A\*03:01-KRAS G12V (one at 1.5 nM apparent, with single-residue peptide discrimination). **Chai-2 is closed-source — not available to us.**