"""Intrinsic (Rho-independent) terminator detection — local ARNold substitute.

ARNold (used in the paper's `CDS_ARNold.ipynb`) is a Selenium-driven webserver and
cannot be run headlessly here. This module reproduces the *intent* — find a stable
stem-loop immediately followed by a poly-U tract — with two local components:

  1. the paper's own poly-U regex (from `RNA-seq/Read_Classification_RegExp.ipynb`),
  2. a ViennaRNA hairpin ΔG for the stem just upstream of that poly-U tract.

A terminator is called when a poly-U tract is found AND the sequence immediately
upstream folds into a hairpin whose closing stem sits within `max_gap` nt of the
tract, with stem ΔG at or below `dg_threshold`.

This is a documented substitution, not an exact ARNold reproduction (see DESIGN.md §4,
REPORT.md). ΔG threshold and window are tunable and should be calibrated against the
fluoride `Measured_Terminators.fasta` ground truth.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import RNA  # ViennaRNA

# The exact poly-U pattern used in the paper (DNA/T alphabet).
POLYU_PATTERN = re.compile(r"TTTTT|T[AGC]TTTT|TT[AGC]TTT|TTT[AGC]TT")


@dataclass
class TerminatorHit:
    found: bool = False
    dG: float | str = ""              # kcal/mol of the folded hairpin window ("" if none)
    start: int | str = ""             # 0-based start of hairpin in ep_seq
    end: int | str = ""               # 0-based end (poly-U end) in ep_seq
    polyU: str = ""                   # matched poly-U tract
    structure: str = ""               # dot-bracket of the stem window (audit)
    reason: str = ""                  # why not found, when found is False

    def coords(self) -> str:
        if not self.found:
            return ""
        return f"{self.start}-{self.end}"


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U").replace(" ", "").replace("\n", "")


def _count_pairs(structure: str) -> int:
    return structure.count("(")


def _last_close_index(structure: str) -> int:
    idx = structure.rfind(")")
    return idx


def find_terminator(
    ep_seq: str,
    dg_threshold: float = -8.0,
    upstream: int = 45,
    max_gap: int = 6,
    min_pairs: int = 4,
) -> TerminatorHit:
    """Detect the first credible intrinsic terminator in `ep_seq`.

    Parameters
    ----------
    ep_seq        expression-platform (extended downstream) sequence, DNA or RNA.
    dg_threshold  max stem ΔG (kcal/mol) to accept as a real hairpin.
    upstream      how many nt upstream of the poly-U to fold as the stem window.
    max_gap       max nt between the hairpin's closing pair and the poly-U start.
    min_pairs     minimum base pairs required in the stem.
    """
    if not ep_seq:
        return TerminatorHit(reason="empty ep_seq")

    dna = ep_seq.upper().replace("U", "T").replace(" ", "").replace("\n", "")

    best: TerminatorHit | None = None
    for m in POLYU_PATTERN.finditer(dna):
        pu_start, pu_end = m.start(), m.end()
        win_start = max(0, pu_start - upstream)
        window = dna[win_start:pu_start]
        if len(window) < 2 * min_pairs + 3:  # too short to hold a stem+loop
            continue

        structure, mfe = RNA.fold(_to_rna(window))
        n_pairs = _count_pairs(structure)
        last_close = _last_close_index(structure)
        gap = (len(window) - 1 - last_close) if last_close >= 0 else 999

        is_term = (
            n_pairs >= min_pairs
            and last_close >= 0
            and gap <= max_gap
            and mfe <= dg_threshold
        )
        if is_term:
            hit = TerminatorHit(
                found=True,
                dG=round(float(mfe), 2),
                start=win_start,
                end=pu_end,
                polyU=m.group(0),
                structure=structure,
            )
            # Keep the most stable (most negative ΔG) terminator.
            if best is None or hit.dG < best.dG:
                best = hit

    if best is not None:
        return best

    # No terminator: say why, but still report the best poly-U tract seen (if any).
    first_pu = POLYU_PATTERN.search(dna)
    if first_pu is None:
        return TerminatorHit(reason="no poly-U tract")
    return TerminatorHit(
        polyU=first_pu.group(0),
        reason="poly-U present but no stable upstream hairpin",
    )
