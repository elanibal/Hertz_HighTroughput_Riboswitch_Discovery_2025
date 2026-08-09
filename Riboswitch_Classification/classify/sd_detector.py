"""Shine-Dalgarno (SD) detection and sequestration test for translational calls.

Two questions, answered separately:

  * `has_SD`         — is there an AGGAGG-like ribosome-binding motif in the window
                       upstream of the downstream start codon?
  * `anti_SD_overlap`— is that SD *base-paired (sequestered)* by an aptamer/EP helix
                       in the MFE fold? This is the mechanistic hallmark of a
                       translational riboswitch. A bare, accessible SD is NOT enough
                       (every ORF has one) — see DESIGN.md §2.1.

Sequestration is an MFE single-structure proxy (DESIGN.md §4, flag 1): a rigorous
call needs apo vs holo structures. Translational calls are therefore hypotheses.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import RNA  # ViennaRNA

# SD variants, longest/strongest first so we prefer the most complete motif.
# Complementary to the 16S 3' anti-SD tail (3'-...UCCUCC-5').
_SD_VARIANTS = ["AGGAGG", "AAGGAG", "GGAGGU", "AGGAGA", "AGGAG", "GGAGG", "GAGGA", "GGAG"]
_SD_REGEX = re.compile("|".join(_SD_VARIANTS))


@dataclass
class SDHit:
    found: bool = False
    motif: str = ""
    offset_to_start: int | str = ""   # nt between SD end and start codon (spacer)
    sequestered: bool = False
    paired_fraction: float | str = ""  # fraction of SD nt paired in MFE fold
    structure: str = ""
    reason: str = ""


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U").replace(" ", "").replace("\n", "")


def detect_sd(
    ep_seq: str,
    nts_to_start_codon,
    aptamer_seq: str = "",
    window_lo: int = 20,
    window_hi: int = 3,
    fold_context: int = 60,
    seq_frac_for_sequestered: float = 0.5,
) -> SDHit:
    """Find an SD upstream of the start codon and test whether it is sequestered.

    Parameters
    ----------
    ep_seq              expression platform (starts at aptamer 3' end).
    nts_to_start_codon  distance (nt) from ep_seq start to the downstream start codon.
    aptamer_seq         aptamer, used only to provide 5' folding context.
    window_lo/window_hi SD is searched in [start-window_lo, start-window_hi].
    fold_context        nt of aptamer 3' tail prepended when folding for sequestration.
    """
    if not ep_seq:
        return SDHit(reason="empty ep_seq")

    ep = ep_seq.upper().replace("U", "T")

    # Locate the start codon within ep_seq.
    try:
        start_idx = int(nts_to_start_codon)
    except (TypeError, ValueError):
        start_idx = -1

    if start_idx < 0 or start_idx > len(ep):
        # Unknown / out-of-range start: search the whole EP but flag it.
        search_region = (0, len(ep))
        note = "start codon position unknown; searched whole EP"
    else:
        lo = max(0, start_idx - window_lo)
        hi = max(0, start_idx - window_hi)
        search_region = (lo, hi)
        note = ""

    region = ep[search_region[0]:search_region[1]]
    # Prefer the SD closest to the start codon (last match in the region).
    matches = list(_SD_REGEX.finditer(region))
    if not matches:
        return SDHit(reason="no SD motif in window" + (f"; {note}" if note else ""))
    m = matches[-1]
    sd_abs_start = search_region[0] + m.start()
    sd_abs_end = search_region[0] + m.end()
    offset = (start_idx - sd_abs_end) if start_idx >= 0 else ""

    # --- Sequestration test: fold aptamer tail + EP up to the start codon ---------
    tail = (aptamer_seq or "").upper().replace("U", "T")[-fold_context:]
    fold_stop = start_idx if start_idx >= 0 else len(ep)
    fold_dna = tail + ep[:fold_stop]
    sd_in_fold_start = len(tail) + sd_abs_start
    sd_in_fold_end = len(tail) + sd_abs_end

    sequestered = False
    paired_fraction: float | str = ""
    structure = ""
    if len(fold_dna) >= 4 and sd_in_fold_end <= len(fold_dna):
        structure, _mfe = RNA.fold(_to_rna(fold_dna))
        sd_struct = structure[sd_in_fold_start:sd_in_fold_end]
        paired = sum(1 for ch in sd_struct if ch in "()")
        paired_fraction = round(paired / max(1, len(sd_struct)), 2)
        sequestered = paired_fraction >= seq_frac_for_sequestered

    return SDHit(
        found=True,
        motif=m.group(0),
        offset_to_start=offset,
        sequestered=sequestered,
        paired_fraction=paired_fraction,
        structure=structure,
        reason=note,
    )
