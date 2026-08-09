# Riboswitch_Classification

Extends the Hertz et al. (2026, *NAR* gkag542) fluoride covariation/classification
pipeline to **any** metabolite-sensing riboswitch class. Produces one machine-readable
row per riboswitch sequence (aptamer + expression platform) along three axes:

1. **Mechanism** — transcriptional vs translational (positive, sequence-based call).
2. **Direction** — ON vs OFF (annotation-driven, evidence recorded).
3. **Sensitivity** — placeholder columns for IVT→NGS Δ%AT / GreB data (not measured here).

See [`DESIGN.md`](DESIGN.md) for the frozen per-sequence schema and decision rules,
and `REPORT.md` (after the TPP pilot) for results and honest limits.

## Layout
```
classify/
  schema.py       canonical frozen column schema (shared across all classes)
  terminator.py   intrinsic-terminator detector  (ViennaRNA hairpin ΔG + repo poly-U regex; ARNold substitute)
  sd_detector.py  Shine-Dalgarno + sequestration test  (ViennaRNA MFE)
  mechanism.py    combines terminator + SD -> transcriptional/translational/ambiguous/other
  direction.py    ON/OFF keyword rubric on downstream gene product (class keywords extensible)
  classifier.py   classify_sequence(): pure, network-free row builder
  ncbi.py         Rfam parse/download + NCBI EP extension & gene annotation (the only network I/O)
run_classification.py   CLI: Rfam family -> classification CSV
tests/          offline unit tests (no network)  -> `python tests/test_classifiers.py`
```

## Reuse of the paper's code
- EP extension & Rfam header parsing: from `Bioinformatics_Covariation/Prepping_sequences.ipynb`.
- Downstream-gene / start-codon lookup: from `Bioinformatics_Covariation/NCBI_CDS.ipynb`.
- Poly-U terminator regex: from `RNA-seq/Read_Classification_RegExp.ipynb`.
- ARNold (webserver, `CDS_ARNold.ipynb`) is replaced by a local ViennaRNA detector
  (documented in `DESIGN.md` §4 and `REPORT.md`).

## Quick start
```bash
python tests/test_classifiers.py          # offline sanity check (no network)
# pilot (rate-limited NCBI pull — authorise first, keep --limit small):
python run_classification.py --rfam RF00059 --fetch-seed --email you@example.com \
    --limit 5 --out pilot.csv
```
