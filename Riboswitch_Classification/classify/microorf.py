"""Micro-ORF (uORF) scanner for the riboswitch leader region.

Small ORFs in the 5' leader (aptamer + expression platform, upstream of the main
downstream gene) are functionally relevant to riboswitch regulation — e.g. leader
peptides / uORFs whose translation couples to the aptamer's structural switch, and
short ORFs that overlap a sequestered Shine-Dalgarno. Requested by A. Arce for the
TPP pilot; recorded as a first-class schema field so it is captured for every class.

Definition: an ORF is a start codon (ATG, and the bacterial alternatives GTG/TTG)
followed by an in-frame stop (TAA/TAG/TGA), with peptide length in [min_aa, max_aa].
Scanned on the sense strand in all three frames, within the leader only.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

_START = {"ATG", "GTG", "TTG"}
_STOP = {"TAA", "TAG", "TGA"}


@dataclass
class MicroORF:
    start: int          # 0-based start (of start codon) in the scanned leader
    end: int            # 0-based end (exclusive, after the stop codon)
    aa_len: int         # peptide length in amino acids (excl. stop)
    start_codon: str

    def __str__(self) -> str:
        return f"{self.start}-{self.end}({self.aa_len}aa,{self.start_codon})"


@dataclass
class MicroORFResult:
    orfs: list = field(default_factory=list)

    @property
    def count(self) -> int:
        return len(self.orfs)

    def summary(self) -> str:
        if not self.orfs:
            return ""
        return f"{self.count} uORF: " + ";".join(str(o) for o in self.orfs)


def find_microorfs(
    leader: str,
    min_aa: int = 2,
    max_aa: int = 50,
    allow_alt_start: bool = True,
) -> MicroORFResult:
    """Find micro-ORFs in `leader` (sense strand, 3 frames).

    Overlapping ORFs in different frames are all reported; within a frame the first
    in-frame stop closes each ORF. Returns them sorted by start position.
    """
    if not leader:
        return MicroORFResult([])
    s = leader.upper().replace("U", "T")
    starts = _START if allow_alt_start else {"ATG"}
    orfs: list[MicroORF] = []

    for frame in range(3):
        i = frame
        while i + 3 <= len(s):
            codon = s[i:i + 3]
            if codon in starts:
                # walk to the next in-frame stop
                j = i + 3
                while j + 3 <= len(s):
                    c2 = s[j:j + 3]
                    if c2 in _STOP:
                        aa = (j - i) // 3            # codons before the stop
                        if min_aa <= aa <= max_aa:
                            orfs.append(MicroORF(i, j + 3, aa, codon))
                        i = j                        # continue scan after this ORF's stop
                        break
                    j += 3
                else:
                    break  # no stop before leader end in this frame
            i += 3

    orfs.sort(key=lambda o: (o.start, o.end))
    return MicroORFResult(orfs)
