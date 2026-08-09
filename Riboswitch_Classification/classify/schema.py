"""Canonical per-sequence riboswitch classification schema.

The column list, order, and enum vocabularies defined here are FROZEN and shared
across every riboswitch class (see ../DESIGN.md). Adding a new metabolic class must
never require editing this file. `empty_row()` returns a dict with every column
present and defaulted, so no downstream code can silently drop a column.
"""

from __future__ import annotations

# ---- Frozen column order (see DESIGN.md §1) -------------------------------------
COLUMNS = [
    # identity & provenance
    "rfam_family", "aptamer_id", "accession", "genome_coords", "strand",
    # sequence content
    "aptamer_seq", "ep_seq", "ep_len",
    # downstream gene
    "downstream_gene_product", "downstream_gene", "nts_to_start_codon",
    # mechanism axis
    "has_intrinsic_terminator", "terminator_dG", "terminator_coords", "polyU_seq",
    "has_SD", "sd_seq", "anti_SD_overlap", "mechanism_call", "mechanism_evidence",
    # direction axis
    "predicted_direction", "direction_evidence", "direction_basis",
    # confidence & notes
    "confidence", "notes",
    # placeholders for future measured data
    "pct_AT_ligand", "pct_AT_apo", "delta_pct_AT", "greB_delta",
    "sensitivity_bin", "data_source",
]

# ---- Controlled vocabularies ---------------------------------------------------
MECHANISM_VALUES = {"transcriptional", "translational", "ambiguous", "other"}
DIRECTION_VALUES = {"ON", "OFF", "unknown"}
DIRECTION_BASIS_VALUES = {"annotation", "architecture", "both", "none"}
CONFIDENCE_VALUES = {"high", "medium", "low"}
SENSITIVITY_BINS = {"low", "mid", "high", ""}  # "" until measured


def empty_row() -> dict:
    """A schema-complete row with safe defaults (missing = '' , never NaN/0)."""
    row = {c: "" for c in COLUMNS}
    row.update(
        ep_len=0,
        nts_to_start_codon="",
        has_intrinsic_terminator=False,
        has_SD=False,
        anti_SD_overlap=False,
        mechanism_call="other",
        predicted_direction="unknown",
        direction_basis="none",
        confidence="low",
    )
    return row


def sensitivity_bin_from_delta(delta_pct_at) -> str:
    """Paper's dynamic-range bins from |Δ%AT| (empty if no measurement)."""
    if delta_pct_at == "" or delta_pct_at is None:
        return ""
    d = abs(float(delta_pct_at))
    if d < 40:
        return "low"      # 20-40 %
    if d < 50:
        return "mid"      # 40-50 %
    return "high"         # 50-100 %


def validate_row(row: dict) -> list[str]:
    """Return a list of schema violations (empty = valid). Used by tests."""
    problems = []
    missing = [c for c in COLUMNS if c not in row]
    if missing:
        problems.append(f"missing columns: {missing}")
    if row.get("mechanism_call") not in MECHANISM_VALUES:
        problems.append(f"bad mechanism_call: {row.get('mechanism_call')!r}")
    if row.get("predicted_direction") not in DIRECTION_VALUES:
        problems.append(f"bad predicted_direction: {row.get('predicted_direction')!r}")
    if row.get("direction_basis") not in DIRECTION_BASIS_VALUES:
        problems.append(f"bad direction_basis: {row.get('direction_basis')!r}")
    if row.get("confidence") not in CONFIDENCE_VALUES:
        problems.append(f"bad confidence: {row.get('confidence')!r}")
    if row.get("sensitivity_bin") not in SENSITIVITY_BINS:
        problems.append(f"bad sensitivity_bin: {row.get('sensitivity_bin')!r}")
    return problems
