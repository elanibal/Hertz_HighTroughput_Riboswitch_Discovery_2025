"""Mechanism axis: combine terminator + SD signals into a positive call.

Decision table (DESIGN.md §2.1), with T = intrinsic terminator present,
S = SD present AND sequestered:

    T & !S -> transcriptional
   !T &  S -> translational
    T &  S -> ambiguous
   !T & !S -> other
"""

from __future__ import annotations

from dataclasses import dataclass

from .terminator import TerminatorHit
from .sd_detector import SDHit


@dataclass
class MechanismCall:
    call: str
    evidence: str


def call_mechanism(term: TerminatorHit, sd: SDHit) -> MechanismCall:
    T = bool(term.found)
    S = bool(sd.found and sd.sequestered)

    term_ev = (
        f"terminator ΔG={term.dG} at {term.coords()}" if T
        else f"no terminator ({term.reason})"
    )
    if sd.found:
        sd_ev = (
            f"SD '{sd.motif}' sequestered (paired={sd.paired_fraction})" if sd.sequestered
            else f"SD '{sd.motif}' present but accessible (paired={sd.paired_fraction})"
        )
    else:
        sd_ev = f"no SD ({sd.reason})"

    if T and not S:
        call = "transcriptional"
    elif S and not T:
        call = "translational"
    elif T and S:
        call = "ambiguous"
    else:
        call = "other"

    return MechanismCall(call, f"{term_ev}; {sd_ev}")
