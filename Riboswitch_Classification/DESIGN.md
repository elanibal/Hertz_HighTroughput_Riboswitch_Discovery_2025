# DESIGN.md — Per-sequence riboswitch classification schema

**Module:** `Riboswitch_Classification/`
**Scope:** Extend the Hertz et al. (2026, *NAR* gkag542) covariation/classification
pipeline from the fluoride riboswitch (RF01734) to *any* metabolite-sensing
riboswitch class, producing one machine-readable row per riboswitch sequence
along three mechanistic axes.

This document defines a **single table specification that is identical across all
riboswitch classes**, so that adding cobalamin, FMN, THF, SAM, etc. later requires
**no schema change and no re-run** of already-classified families.

---

## 0. Design principles

1. **One row = one riboswitch sequence** (aptamer + expression platform), keyed by
   `(rfam_family, aptamer_id)`.
2. **Positive, mechanistic calls** — not negative proxies. The paper labels a
   variant "translational" only when `%AT > 95%` (i.e. it *never* terminates in
   IVT, an *inert*/negative proxy). Here the mechanism call is made from sequence
   architecture (terminator present → transcriptional; SD sequestrable by an
   aptamer helix → translational), so it works with **no wet-lab data**.
3. **Evidence, not just labels** — every call carries an evidence string and a
   confidence level so a human (or a later covariation/wet-lab result) can audit it.
4. **No silent truncation** — sequences that cannot be classified are *kept* with
   `mechanism_call = other` / `predicted_direction = unknown` and a `notes` reason.
5. **Forward-compatible with measured data** — placeholder columns hold IVT→NGS
   `%AT` / `Δ%AT` / GreB fields so the same table absorbs experimental results when
   they exist, without a re-run.
6. **Honest limits** — the IVT→NGS assay (and eROSALIND/EcRNAP porting) reads
   **transcription only**. Translational members are *flagged, not characterized*.
   Any report built from this table must state this.

---

## 1. Column specification

Column order is canonical. Types: `str` unless noted. Missing/not-applicable values
are the empty string `""` (never a silent `0`/`NaN`), except numeric ΔG which is
`""` when no terminator was found.

### 1.1 Identity & provenance

| # | Column | Type | Justification (one line) |
|---|--------|------|--------------------------|
| 1 | `rfam_family` | str | Rfam accession of the aptamer class (e.g. `RF00059`); makes the table self-describing when families are concatenated. |
| 2 | `aptamer_id` | str | Stable per-sequence key (Rfam seq name, e.g. `AE017225.1/…`); primary key together with `rfam_family`. |
| 3 | `accession` | str | NCBI nucleotide accession the sequence was pulled from; needed to re-fetch / audit coordinates. |
| 4 | `genome_coords` | str | `start-stop(strand)` of the aptamer in `accession`; anchors every downstream coordinate and lets anyone reproduce the EP extension. |
| 5 | `strand` | str | `+` / `-`; the EP extension and CDS search run in opposite genomic directions, so strand must be explicit (mirrors `Prepping_sequences.ipynb`). |

### 1.2 Sequence content

| # | Column | Type | Justification |
|---|--------|------|---------------|
| 6 | `aptamer_seq` | str | The Rfam aptamer sequence; the ligand-binding module whose helices the EP competes with. |
| 7 | `ep_seq` | str | Expression-platform sequence downstream of the aptamer (aptamer 3′ end → downstream start codon / terminator), captured exactly as in `Prepping_sequences.ipynb`. |
| 8 | `ep_len` | int | Length of `ep_seq`; a very short EP means the terminator/SD may be truncated → lowers confidence. |

### 1.3 Downstream gene (direction evidence source)

| # | Column | Type | Justification |
|---|--------|------|---------------|
| 9 | `downstream_gene_product` | str | NCBI CDS `product` qualifier of the first downstream ORF (from `NCBI_CDS.ipynb`); the primary evidence for ON vs OFF. |
| 10 | `downstream_gene` | str | Gene symbol/locus tag when available; disambiguates generic products (e.g. `thiC` vs "hypothetical protein"). |
| 11 | `nts_to_start_codon` | int | Distance (nt) from aptamer 3′ end to the downstream start codon; defines the window in which the SD and terminator must fall. |

### 1.4 Mechanism axis — transcriptional vs translational

| # | Column | Type | Justification |
|---|--------|------|---------------|
| 12 | `has_intrinsic_terminator` | bool | True if an intrinsic terminator (stable hairpin + poly-U tract) is predicted in the EP; the *positive* signal for **transcriptional** control. |
| 13 | `terminator_dG` | float\|"" | Hairpin stem ΔG (kcal/mol, ViennaRNA) — stability of the terminator; also a **sensitivity** proxy and feeds `confidence`. |
| 14 | `terminator_coords` | str | `start-stop` of the predicted terminator within `ep_seq`; lets us later ask *which aptamer helix it overlaps* (ON/OFF architecture). |
| 15 | `polyU_seq` | str | The matched poly-U tract (repo regex); records exactly what triggered the terminator call, for audit. |
| 16 | `has_SD` | bool | True if a Shine–Dalgarno motif (AGGAGG-like) is found in the window upstream of the start codon; prerequisite for a **translational** call. |
| 17 | `sd_seq` | str | The matched SD motif + its offset from the start codon; audit trail for `has_SD`. |
| 18 | `anti_SD_overlap` | bool | True if the SD is predicted to be **base-paired (sequestered)** by an aptamer/EP helix in the MFE fold — the mechanistic hallmark of a translational riboswitch. |
| 19 | `mechanism_call` | enum | `transcriptional` \| `translational` \| `ambiguous` \| `other`; the axis-1 result (rules in §2). |
| 20 | `mechanism_evidence` | str | Human-readable justification (e.g. `terminator ΔG=-12.3, no sequestered SD`); never just the label. |

### 1.5 Direction axis — ON vs OFF

| # | Column | Type | Justification |
|---|--------|------|---------------|
| 21 | `predicted_direction` | enum | `ON` \| `OFF` \| `unknown`; ligand-induced expression change (sign of Δ%AT in measured data). |
| 22 | `direction_evidence` | str | Matched rubric keyword + gene product (e.g. `biosynthesis:"thiamine biosynthesis protein ThiC"→OFF`); the *reason*, per task requirement. |
| 23 | `direction_basis` | enum | `annotation` \| `architecture` \| `both` \| `none`; records *how* the direction was called so annotation-only calls can be down-weighted. |

### 1.6 Confidence & notes

| # | Column | Type | Justification |
|---|--------|------|---------------|
| 24 | `confidence` | enum | `high` \| `medium` \| `low`; combines signal strength across axes so downstream users can threshold. |
| 25 | `notes` | str | Free text: dropped-reason, ambiguity, truncation, multiple terminators, etc. Guarantees **no silent loss**. |

### 1.7 Placeholders for future measured data (empty until IVT→NGS exists)

| # | Column | Type | Justification |
|---|--------|------|---------------|
| 26 | `pct_AT_ligand` | float\|"" | Measured % anti-termination **with** ligand (e.g. +NaF); filled from NGS later. |
| 27 | `pct_AT_apo` | float\|"" | Measured % anti-termination **without** ligand; the apo baseline. |
| 28 | `delta_pct_AT` | float\|"" | `%AT_ligand − %AT_apo`; **sign = measured direction**, **magnitude = sensitivity**. Validates axes 2 & 3. |
| 29 | `greB_delta` | float\|"" | Shift in Δ%AT caused by GreB; a >20 % shift with no sequence motif implicated cotranscriptional **kinetics** (not aptamer K_D) in 26 % of fluoride variants — a sensitivity determinant. |
| 30 | `sensitivity_bin` | enum\|"" | `low` (20–40 %) \| `mid` (40–50 %) \| `high` (50–100 %); dynamic-range bin from `abs(delta_pct_AT)` (paper's bins). Empty until measured. |
| 31 | `data_source` | str | Where measured values came from (dataset/run id) when filled; `""` = computational-only row. |

---

## 2. Decision rules (implemented in `classify/`)

### 2.1 Mechanism (`mechanism_call`)

Let `T` = `has_intrinsic_terminator`, `S` = `has_SD and anti_SD_overlap`
(a *sequestered* SD, not merely a present one).

| T | S | `mechanism_call` | rationale |
|---|---|------------------|-----------|
| ✔ | ✗ | `transcriptional` | intrinsic terminator = Rho-independent transcription attenuation. |
| ✗ | ✔ | `translational` | SD sequestered by a ligand-responsive helix = translation initiation control. |
| ✔ | ✔ | `ambiguous` | both signals present (dual / tandem control, or mis-prediction) — flag, do not force. |
| ✗ | ✗ | `other` | neither signal — truncated EP, novel architecture, or fetch gap; kept with a `notes` reason. |

> A bare SD **without** predicted sequestration (`has_SD` ✔, `anti_SD_overlap` ✗) is
> *not* sufficient for `translational` — nearly every ORF has an SD. Sequestration
> by an aptamer-coupled helix is what makes it *regulatory*.

### 2.2 Direction (`predicted_direction`)

Annotation-first rubric on `downstream_gene_product` / `downstream_gene`
(keyword lists live in `classify/direction.py`, easily extended per class):

- **OFF** (feedback repression) — biosynthesis/salvage: `synthase`, `synthetase`,
  `biosynthesis`, `biosynthetic`, and class-specific enzyme families (e.g. TPP:
  `thiC/thiD/thiE/thiM`, "thiamine biosynthesis"). Ligand present ⇒ pathway product
  present ⇒ shut biosynthesis off.
- **ON** — transport/efflux/detox/alarmone: `transporter`, `permease`, `efflux`,
  `importer`, `ABC transporter`, `channel`, `MFS`. Ligand present ⇒ import/expel to
  restore balance ⇒ keep expression on. (Fluoride = ON: `crcB`/`eriC` efflux.)
- **unknown** — no product, `hypothetical protein`, or no rubric hit. **Not guessed.**

`direction_basis` = `annotation` when only the gene rubric fired; `architecture`
when only terminator/SD-vs-aptamer-helix overlap informed it; `both` when they
agree; `none` when direction is `unknown`. (Architecture-based direction — which
helix the terminator overlaps — is recorded as a `notes`/evidence hint in this
version; the annotation rubric is the primary caller. See §4.)

### 2.3 Confidence

- `high` — mechanism signal is unambiguous (clean terminator ΔG ≤ −8 kcal/mol *or*
  clearly sequestered SD) **and** direction rubric hit a specific keyword.
- `medium` — one axis strong, the other weak/annotation-only, or `ep_len` short.
- `low` — `mechanism_call ∈ {ambiguous, other}` *or* `predicted_direction = unknown`.

---

## 3. Cross-class invariants

- Columns 1–31, names, order, and enums are **frozen**. New classes populate the
  same table; concatenation across classes is a plain `pd.concat`.
- Class-specific knowledge lives **only** in the keyword lists of
  `classify/direction.py` and (optionally) the SD window width — never in the schema.
- Rfam accessions to be added later (confirm before pulling): cobalamin `RF00174`,
  FMN `RF00050`, THF/tetrahydrofolate `RF01831`, SAM `RF00162`, purine `RF00167`,
  ZTP `RF01750`, lysine `RF00168`, TPP `RF00059`, glmS `RF00083`.
  **These are to be verified against Rfam before any pull; only RF00059 is used this session.**

---

## 4. Things I am unsure how to compute (flagged for review)

1. **`anti_SD_overlap` (SD sequestration).** Computed here as: fold
   `aptamer_3′tail + ep_seq` with ViennaRNA and test whether the SD nucleotides are
   base-paired in the MFE structure to an upstream (aptamer/EP) segment. This is a
   **single-structure MFE proxy**, not a two-state (apo/holo) model. A rigorous call
   needs the ligand-bound vs free structures (as the R-scape apo/holo models give for
   fluoride). Treat `translational` calls as *hypotheses*. **Flagged.**

2. **Which aptamer helix the terminator overlaps (ON/OFF architecture).** Truly
   deciding ON vs OFF from architecture requires the covariation model's helix
   annotation (P1/anti-terminator) aligned to each sequence. Without per-sequence
   structure mapping, this version uses architecture only as a *secondary hint* and
   relies on the annotation rubric as the primary direction caller. **Flagged** —
   full architecture-based direction should consume the existing R-scape/CaCoFold
   models (a follow-up, not this session).

3. **SD consensus & spacing per organism.** SD strength and optimal spacing vary by
   species/16S tail. A single AGGAGG-like regex + fixed window will miss weak/kingdom-
   variant SDs. Window width and motif are parameters in `sd_detector.py`. **Flagged.**

4. **Terminator ΔG threshold.** The ΔG cutoff separating "real" terminators from
   incidental hairpins is a tunable heuristic (default ≤ −8 kcal/mol for the hairpin
   stem, plus a poly-U tract). Substitutes for ARNold's trained model; see REPORT.md.
   **Flagged** — should be calibrated against the fluoride measured terminators
   (`Measured_Terminators.fasta`) as ground truth.

5. **ARNold substitution.** ARNold is a webserver (Selenium in `CDS_ARNold.ipynb`)
   and is **not scriptable headlessly here** (selenium/chromedriver absent). It is
   replaced by a local ViennaRNA hairpin-ΔG + repo poly-U-regex detector. Documented
   in REPORT.md. **Flagged** as a substitution, not an exact reproduction.
