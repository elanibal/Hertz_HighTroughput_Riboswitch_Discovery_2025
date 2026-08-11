# REPORT — SAM-I (RF00162) classification pilot

**Family:** Rfam **RF00162** (SAM riboswitch / S box leader) — **SEED, first N = 150** of 457 sequences
**Output:** [`sam_RF00162_classification.csv`](sam_RF00162_classification.csv) (same 34-column schema as TPP)
**Status:** **computational prediction only — no wet-lab data.** Hypotheses to validate.

SAM-I was chosen as a **bacterial-only** contrast to TPP (which is eukaryote-capable), to test
whether the pipeline cleanly recovers a family known to be predominantly **transcriptional**.

---

## 1. Result — SAM-I is strongly transcriptional

| axis | result |
|------|--------|
| **Intrinsic terminator** | **136 / 150 (91 %)** — the headline; SAM-I reads as a transcriptional (attenuation) family, as expected for the S box |
| terminator ΔG | median **−14.8** kcal/mol (range −29.1 … −8.0) — strong hairpins |
| mechanism | transcriptional 90 · ambiguous 46 · other 11 · translational 3 |
| direction | OFF 60 · ON 38 · unknown 52 |

- Only **3 translational** and **11 other** — the family resolves cleanly, unlike TPP.
- **OFF calls are solid SAM biology:** MetK / S-adenosylmethionine synthetase, methionine
  synthase (MetE/MetH), O-acetylhomoserine sulfhydrylase, cystathionine synthase, MetF,
  methylthio-ribose salvage (Mtn) — all methionine/SAM/sulfur biosynthesis → feedback OFF.

## 2. Honest caveat — ON is over-called

Of the **38 ON** calls, only **5 are SAM-specific** (metN / metQ / D-methionine ABC
transporters — genuine SAM-ON methionine-import targets). The other **33 are generic
transporters** (ribose, oligopeptide, metal-ion, Na⁺/H⁺ antiporter, multidrug) that the
keyword rubric flags ON because they are transport annotations, but which **may not be the
SAM-regulated gene** (the first downstream CDS is not always the regulated one). Treat ON as
an **upper bound**; the high-confidence ON set is the 5 methionine-specific transporters
(`direction_evidence` contains `RF00162` for those). This is a limitation of annotation-only
direction calling, not a folding error.

`unknown` = 52 (22 no annotated CDS + 30 uninformative products, e.g. "hypothetical protein",
chemotaxis, luciferase-family — genuinely not methionine-related or unannotated).

## 3. TPP vs SAM — the pipeline separates them as it should

| | **TPP (RF00059)** | **SAM-I (RF00162)** |
|---|---|---|
| N (this pilot) | 115 (full seed) | 150 (of 457) |
| **intrinsic terminator** | **55 %** | **91 %** |
| `other` (unresolved) | 32 % | **7 %** |
| translational (flagged) | 13 % | 2 % |
| kingdom | bacteria + **10 % eukaryotic** | bacterial (Firmicutes-heavy) |
| character | **mixed**, eukaryote-capable | **cleanly transcriptional** |

This is the intended contrast: SAM-I is a tight, bacterial, terminator-driven family, and the
much lower `other` fraction supports the earlier finding that **TPP's high `other`/`unknown`
was partly its eukaryotic members and sparser annotation, not a pervasive pipeline failure.**

## 4. Same limits as the TPP report apply

Transcription-only measurability (translational flagged not characterized), MFE-SD
sequestration is an upper bound, direction is annotation-driven (here that specifically
inflates ON), no measured data. Leader micro-ORFs (149/150 leaders carry ≥1 candidate) and
`cds_extends_beyond_window` (126/150) behave as in TPP — a screen and an
expected-common flag, respectively. See `REPORT.md` §4 for the full limits discussion.

## 5. Reproduce
```bash
python run_classification.py --rfam RF00162 --fetch-seed \
    --email you@example.com --limit 150 --out sam_RF00162_classification.csv
```
