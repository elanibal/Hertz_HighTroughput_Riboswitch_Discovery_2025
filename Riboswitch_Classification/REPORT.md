# REPORT.md — TPP (RF00059) classification pilot

**Family:** Rfam **RF00059** (TPP / THI element) — **SEED alignment, N = 115 sequences**
**Output:** [`tpp_RF00059_classification.csv`](tpp_RF00059_classification.csv) (34-column schema, `DESIGN.md`)
**Status:** **computational prediction only — no wet-lab data.** These are hypotheses to be
validated, not measured results. Read the *Honest limits* section before quoting any number.

TPP was chosen as the pilot because it already has a covariation model and appears in the
paper's Figure 2, so it is the cleanest first test of the generalized classifier.

---

## 1. What was run

1. Downloaded the RF00059 **seed** Stockholm from Rfam (115 sequences, 104 unique
   accessions, +68 / −47 genomic strands, aptamer length 77–148 nt).
2. For each aptamer: extended downstream via NCBI to capture the expression platform (EP),
   located the first downstream CDS (product + start-codon offset), and classified along the
   two computational axes (mechanism, direction) plus leader micro-ORFs.
3. Rate-limited NCBI pull was **authorised before running** (per task constraint).

### ARNold substitution (documented)

ARNold is a Selenium-driven **webserver** (`CDS_ARNold.ipynb`) and cannot run headlessly in
this environment. It was replaced by a **local intrinsic-terminator detector** =
**ViennaRNA hairpin ΔG** (`RNA.fold`, v2.7.2) + the paper's own **poly-U regex** (from
`Read_Classification_RegExp.ipynb`). A terminator is called when a poly-U tract has an
upstream hairpin closing within 6 nt of it and stem ΔG ≤ −8 kcal/mol. This is a documented
substitution, **not** an exact reproduction of ARNold's trained model; the ΔG threshold is a
tunable heuristic that should be calibrated against the fluoride `Measured_Terminators.fasta`
ground truth.

---

## 2. How the TPP family splits (predictions)

### Mechanism — transcriptional vs translational

| call | n | % | basis |
|------|---|---|-------|
| **transcriptional** | 34 | 30 % | intrinsic terminator, no sequestered SD |
| **ambiguous** | 29 | 25 % | terminator **and** MFE-sequestered SD both present |
| **translational** | 15 | 13 % | sequestered SD, no terminator |
| **other** | 37 | 32 % | neither signal (see §3 — mostly data limits) |

- **Intrinsic terminator found in 63 / 115 (55 %)** sequences — the reliable, specific
  positive signal. This is the transcriptional read-out and the one relevant to
  eROSALIND / EcRNAP porting.
- **MFE-sequestered SD in 44 / 115 (38 %)** — the translational signal, but **low
  specificity** (see limits): treat translational + the SD half of ambiguous as an
  **upper bound** on translational membership.

### Terminator strength (ΔG distribution, n = 63)

Median **−13.0 kcal/mol**, range −27.1 … −8.2. Bins: ≤ −20: 4 · ≤ −15: ~17 · ≤ −12: ~36 ·
≤ −10: ~46 · ≤ −8: 63. Strong, well-populated hairpins — consistent with genuine terminators
rather than incidental stems. (ΔG will later double as a **sensitivity** proxy once measured
Δ%AT exists to calibrate it.)

### Direction — ON vs OFF

| call | n | % | interpretation |
|------|---|---|----------------|
| **OFF** | 52 | 45 % | biosynthesis / salvage downstream (thiC/thiD/thiE/thiM/thiO, thiamin-phosphate pyrophosphorylase, TenA/thiaminase) → feedback repression |
| **ON** | 18 | 16 % | transport downstream (thiamine-binding periplasmic / TbpA / YkoF / ABC & MFS transporters, TonB receptor) → import when ligand present |
| **unknown** | 45 | 39 % | no informative annotation (see §3) — **not guessed** |

Among annotated sequences TPP is **predominantly OFF (biosynthesis feedback)** with a
transport-driven **ON minority** — the textbook expectation for a thiamine riboswitch, which
is a sanity check that the annotation rubric is behaving.

### Mechanism × direction

```
                 OFF  ON  unknown
transcriptional   19   5    10
ambiguous          7   6    16
translational     10   3     2
other             16   4    17
```

### Confidence

`high` 37 · `low` 78 (no `medium`: a clean mechanism call with a specific direction is
graded high, everything ambiguous/other/unknown-direction is low — the tiering is honest but
coarse and could be refined).

### New leader/CDS fields (requested by A. Arce)

- **`cds_extends_beyond_window` = 82 / 115.** The downstream CDS *start* was found but the
  ORF does not *close* within the 600 nt search window. This is **expected for most bacterial
  genes** (typical ORFs exceed the ~530 nt of CDS that fits after the aptamer), so the flag
  mainly earmarks sequences for **full-CDS retrieval** if the complete protein is wanted; it
  is not an anomaly. Only ~10 annotated genes closed inside the window (short ORFs).
- **micro-ORFs (uORFs) = 96 / 115 leaders carry ≥ 1** candidate (≥ 8 aa, ATG/GTG/TTG →
  in-frame stop; total 270, mean 2.3, max 11). These are **candidate** leader peptides /
  SD-overlapping ORFs for follow-up — the count is sensitive to the length/alt-start settings
  (ATG-only ≥ 8 aa → 53/115), so treat as a screen, not a call.

---

## 3. Sequences that could not be classified (with reasons — no silent loss)

Every unresolved sequence is **kept** in the CSV with `mechanism_call = other` /
`predicted_direction = unknown` and a reason in `notes`.

- **7 have no EP at all.** 6 are **RNAcentral `URS…` identifiers** (not NCBI nucleotide
  accessions — cannot be extended by efetch); 1 (`AF159589.1`) has the aptamer at the **3′ end
  of a short clone** with no downstream sequence deposited.
- **23 have an EP but no annotated downstream CDS within 600 nt.** Inspection shows these are
  largely **eukaryotic and environmental** accessions (e.g. *Arabidopsis*/rice genomic & cDNA
  `AK…`/`AC…`, marine metagenome `AACY0…`/`AARF0…`). **TPP riboswitches in plants and fungi sit
  in UTRs/introns**, where the bacterial "intrinsic terminator vs SD-sequestration" and
  "biosynthesis vs transport" framework **does not directly apply.** This is a real limitation
  of classifying an all-kingdoms seed with a bacterial model, not a pipeline failure.
- The remainder of `other` are short EPs or sequences with a poly-U tract but no hairpin
  stable enough to pass the ΔG threshold.

---

## 4. Honest limits — do NOT overclaim

1. **Transcription-only measurability.** The IVT→NGS assay (and eROSALIND / EcRNAP, which is
   transcriptional) reads **transcription only**. The 15 `translational` (and the SD-driven
   part of the 29 `ambiguous`) are **flagged, not characterized** — they would look like the
   flat controls in the assay. For downstream porting they are **excluded, not lost**.
2. **Translational calls are an upper bound.** SD sequestration is scored from a **single MFE
   fold** (`RNA.fold`). In GC-rich aptamers an AGGAGG-like motif is very often base-paired by
   chance; requiring an **upstream** pairing partner and ≥ 60 % paired nt trims the worst
   false positives but cannot separate a *regulatory* sequestering helix from an incidental
   one. A rigorous call needs the **apo vs holo** structures (the R-scape/CaCoFold apo/holo
   models already exist for TPP) or wet-lab. **The terminator (transcriptional) call is the
   trustworthy one; translational is a hypothesis.**
3. **Direction is annotation-driven.** ON/OFF comes from the downstream gene's NCBI product
   via a keyword rubric, not from EP architecture. Mis-/under-annotation → `unknown`. Which
   aptamer helix the terminator overlaps (the architectural ON/OFF signal) is **not** yet
   computed per-sequence; that needs the covariation model's helix annotation mapped onto each
   sequence (a follow-up).
4. **Kingdom mismatch.** ~¼ of the seed is eukaryotic/environmental; the bacterial mechanism
   model does not apply to those and they land in `other`/`unknown` honestly.
5. **No measured data.** `pct_AT_*`, `delta_pct_AT`, `greB_delta`, `sensitivity_bin` are empty
   placeholders — the schema is ready to absorb IVT→NGS results without a re-run.

**One-line takeaway:** *Of 115 TPP seed sequences, ~55 % carry a predicted intrinsic
terminator (transcriptional-compatible) and TPP reads predominantly OFF/biosynthesis; a
translational subset is flagged but not measurable here, and ~¼ of the seed is
eukaryotic/environmental and outside the bacterial model. All calls are computational
predictions awaiting validation.*

---

## 5. Reproduce

```bash
cd Riboswitch_Classification
python tests/test_classifiers.py                      # 19/19 offline unit tests

# full pilot (rate-limited NCBI pull — authorise first):
python run_classification.py --rfam RF00059 --fetch-seed \
    --email you@example.com --out tpp_RF00059_classification.csv

# iterate on classifiers without re-pulling (offline, from stored sequences):
python rebuild_from_csv.py tpp_RF00059_classification.csv out.csv \
    --email you@example.com --seed data/RF00059_seed.sto
```

Environment: Python 3.11.7 (paper declares 3.10.9; pure-Python logic is version-agnostic),
ViennaRNA 2.7.2, Biopython 1.85, pandas 2.3.3. No new dependencies required.
