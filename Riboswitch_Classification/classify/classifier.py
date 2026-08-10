"""Top-level per-sequence classifier: sequence + annotation -> one schema row.

This is the network-free heart of the pipeline. It takes already-fetched fields
(aptamer/EP sequence, downstream gene product, start-codon offset) and returns a
schema-complete dict. All I/O (Rfam/NCBI) lives in `ncbi.py`; keeping this pure
makes the classification logic unit-testable with no network (see tests/).
"""

from __future__ import annotations

from . import schema
from .terminator import find_terminator
from .sd_detector import detect_sd
from .mechanism import call_mechanism
from .direction import predict_direction
from .microorf import find_microorfs


def _confidence(mech_call: str, direction, term, ep_len: int) -> str:
    strong_mech = (
        (term.found and term.dG != "" and float(term.dG) <= -8.0)
        or mech_call == "translational"
    )
    specific_dir = direction.direction in {"ON", "OFF"}
    if mech_call in {"ambiguous", "other"} or direction.direction == "unknown":
        base = "low"
    elif strong_mech and specific_dir and ep_len >= 40:
        base = "high"
    else:
        base = "medium"
    return base


def classify_sequence(
    *,
    rfam_family: str,
    aptamer_id: str,
    accession: str,
    genome_coords: str,
    strand: str,
    aptamer_seq: str,
    ep_seq: str,
    downstream_gene_product: str = "",
    downstream_gene: str = "",
    nts_to_start_codon="",
    cds_extends_beyond_window: bool = False,
    dg_threshold: float = -8.0,
) -> dict:
    """Classify one riboswitch sequence into a schema-complete row dict."""
    row = schema.empty_row()
    row.update(
        rfam_family=rfam_family,
        aptamer_id=aptamer_id,
        accession=accession,
        genome_coords=genome_coords,
        strand=strand,
        aptamer_seq=aptamer_seq,
        ep_seq=ep_seq,
        ep_len=len(ep_seq or ""),
        downstream_gene_product=downstream_gene_product,
        downstream_gene=downstream_gene,
        nts_to_start_codon=nts_to_start_codon,
        cds_extends_beyond_window=bool(cds_extends_beyond_window),
    )

    notes = []
    if not ep_seq:
        notes.append("no EP sequence (extension failed)")
    if row["ep_len"] and row["ep_len"] < 30:
        notes.append(f"short EP ({row['ep_len']} nt): terminator/SD may be truncated")
    if cds_extends_beyond_window:
        notes.append("downstream CDS not closed within search window (candidate for extension)")

    # --- Leader micro-ORFs (uORFs): scan aptamer + EP up to the main start codon ----
    try:
        start_idx = int(nts_to_start_codon)
    except (TypeError, ValueError):
        start_idx = -1
    leader = (aptamer_seq or "") + (ep_seq[:start_idx] if start_idx >= 0 else (ep_seq or ""))
    # min_aa=8 suppresses the many chance 2-4 aa ORFs; a candidate-generation setting
    # for leader peptides / SD-overlapping uORFs (validate downstream). See REPORT.md.
    uorfs = find_microorfs(leader, min_aa=8, allow_alt_start=True)
    row["microORF_count"] = uorfs.count
    row["microORF"] = uorfs.summary()

    # --- Axis 1: mechanism -------------------------------------------------------
    term = find_terminator(ep_seq, dg_threshold=dg_threshold)
    sd = detect_sd(ep_seq, nts_to_start_codon, aptamer_seq=aptamer_seq)
    mech = call_mechanism(term, sd)

    row.update(
        has_intrinsic_terminator=term.found,
        terminator_dG=term.dG,
        terminator_coords=term.coords(),
        polyU_seq=term.polyU,
        has_SD=sd.found,
        sd_seq=(f"{sd.motif}@-{sd.offset_to_start}" if sd.found and sd.offset_to_start != "" else sd.motif),
        anti_SD_overlap=bool(sd.found and sd.sequestered),
        mechanism_call=mech.call,
        mechanism_evidence=mech.evidence,
    )

    # --- Axis 2: direction -------------------------------------------------------
    direction = predict_direction(downstream_gene_product, downstream_gene, rfam_family)
    row.update(
        predicted_direction=direction.direction,
        direction_evidence=direction.evidence,
        direction_basis=direction.basis,
    )

    # --- Confidence & notes ------------------------------------------------------
    row["confidence"] = _confidence(mech.call, direction, term, row["ep_len"])
    if mech.call == "other":
        notes.append("mechanism unresolved: neither terminator nor sequestered SD")
    if mech.call == "ambiguous":
        notes.append("both terminator and sequestered SD detected")
    row["notes"] = "; ".join(notes)

    return row
