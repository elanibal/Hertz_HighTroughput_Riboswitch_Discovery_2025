"""Data-acquisition layer: Rfam sequences + NCBI EP extension and gene annotation.

Refactored (not rewritten) from the paper's `Prepping_sequences.ipynb` and
`NCBI_CDS.ipynb`. All network I/O is confined to this module so the classifiers in
the sibling modules stay pure and unit-testable offline.

NOTHING here runs at import time. The functions that hit NCBI (`extend_and_annotate`)
are rate-limited and gated by an explicit `limit` in the CLI, because per-sequence
efetch over a whole Rfam family is a rate-limited pull that must be authorised first
(see task constraints / REPORT.md).
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field

from Bio import Entrez, SeqIO


@dataclass
class RfamRecord:
    aptamer_id: str          # full Rfam seq name, e.g. "AE017225.1/1234-1300"
    accession: str           # NCBI accession, e.g. "AE017225.1"
    start: int               # genome start (as in the Rfam name)
    stop: int                # genome stop
    aptamer_seq: str         # aptamer nucleotides (Rfam)
    strand: str = "+"        # "+" if start<stop else "-"

    @property
    def genome_coords(self) -> str:
        return f"{self.start}-{self.stop}({self.strand})"


# ---------------------------------------------------------------------------
# Rfam parsing (no network) — mirrors Prepping_sequences.ipynb header parsing
# ---------------------------------------------------------------------------
_HDR = re.compile(r">?(?P<acc>\S+?)/(?P<start>\d+)-(?P<stop>\d+)")


def parse_rfam_fasta(path: str) -> list[RfamRecord]:
    """Parse an Rfam FASTA (headers like '>ACC.V/START-STOP ...'). No network."""
    records: list[RfamRecord] = []
    acc = start = stop = None
    seq_lines: list[str] = []

    def _flush():
        if acc is None:
            return
        seq = "".join(seq_lines).strip().replace(" ", "")
        strand = "+" if start < stop else "-"
        records.append(RfamRecord(
            aptamer_id=f"{acc}/{start}-{stop}",
            accession=acc, start=start, stop=stop,
            aptamer_seq=seq, strand=strand,
        ))

    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            _flush()
            m = _HDR.match(line)
            if not m:
                acc = None
                continue
            acc = m.group("acc")
            start = int(m.group("start"))
            stop = int(m.group("stop"))
            seq_lines = []
        else:
            seq_lines.append(line)
    _flush()
    return records


def stockholm_to_records(path: str) -> list[RfamRecord]:
    """Parse an Rfam Stockholm alignment into RfamRecords (ungapped seqs). No network."""
    seqs: dict[str, list[str]] = {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        name, aln = parts
        seqs.setdefault(name, []).append(aln)
    out: list[RfamRecord] = []
    for name, chunks in seqs.items():
        m = _HDR.match(name)
        if not m:
            continue
        seq = "".join(chunks).replace("-", "").replace(".", "").upper().replace("U", "T")
        start, stop = int(m.group("start")), int(m.group("stop"))
        out.append(RfamRecord(
            aptamer_id=name, accession=m.group("acc"),
            start=start, stop=stop, aptamer_seq=seq,
            strand="+" if start < stop else "-",
        ))
    return out


# ---------------------------------------------------------------------------
# Rfam download (one HTTP call) — seed alignment via the Rfam API
# ---------------------------------------------------------------------------
def fetch_rfam_seed_stockholm(rfam_acc: str, dest: str) -> str:
    """Download the SEED alignment (bounded, curated) for a family. One request."""
    import requests
    # Verified working Rfam SEED endpoints (checked Aug 2026). The plain /alignment
    # route returns the seed Stockholm directly; the gzip route is the documented one.
    urls = [
        f"https://rfam.org/family/{rfam_acc}/alignment",
        f"https://rfam.org/family/{rfam_acc}/alignment/stockholm?gzip=1&download=1",
    ]
    last = None
    for u in urls:
        try:
            r = requests.get(u, timeout=60)
            if r.ok and r.content:
                data = r.content
                if data[:2] == b"\x1f\x8b":  # gzip
                    import gzip
                    data = gzip.decompress(data)
                if not data.lstrip().startswith(b"# STOCKHOLM"):
                    last = RuntimeError(f"unexpected content from {u}")
                    continue
                with open(dest, "wb") as f:
                    f.write(data)
                return dest
        except Exception as e:  # noqa: BLE001
            last = e
    raise RuntimeError(f"Rfam seed download failed for {rfam_acc}: {last}")


# ---------------------------------------------------------------------------
# NCBI extension + downstream-gene annotation (RATE-LIMITED PULL)
# ---------------------------------------------------------------------------
@dataclass
class Annotation:
    ep_seq: str = ""
    downstream_gene_product: str = ""
    downstream_gene: str = ""
    nts_to_start_codon: int | str = ""
    cds_extends_beyond_window: bool = False   # start found in-window, end past it
    error: str = ""


class EndOfSequence(Exception):
    """Requested window is past the 3' end of the deposited accession."""


def _efetch_gb(accession: str, seq_start: int, seq_stop: int, retries: int = 3):
    """efetch a genbank slice, retrying transient errors (NCBI 400/429/500 under load).

    A persistent HTTP 400 after retries on a forward window usually means seq_start is
    beyond the record end (short/partial clone) -> raised as EndOfSequence so the caller
    can record it distinctly from a real network failure.
    """
    from urllib.error import HTTPError

    last = None
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="gb", retmode="text",
                seq_start=seq_start, seq_stop=seq_stop,
            )
            rec = SeqIO.read(handle, "genbank")
            handle.close()
            return rec
        except HTTPError as e:
            last = e
            time.sleep(0.6 * (attempt + 1))  # backoff for transient throttling
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(0.6 * (attempt + 1))
    if isinstance(last, HTTPError) and last.code == 400:
        raise EndOfSequence(f"{accession}:{seq_start}-{seq_stop} likely past end")
    raise last


def extend_and_annotate(
    rec: RfamRecord,
    ep_window: int = 150,
    max_search: int = 600,
    delay: float = 0.34,          # ~3 req/s (NCBI limit without API key)
) -> Annotation:
    """Fetch the EP (downstream) sequence and the first downstream CDS for one record.

    Efficient variant of NCBI_CDS.ipynb: one genbank fetch of a downstream window
    (`ep_window` nt) usually yields both the EP sequence and the first CDS (product,
    gene symbol, start-codon offset). If no CDS is found, the window is extended once
    to `max_search`. Returns an Annotation; on failure, `error` is set (never raises
    into the caller — no silent drops, the pipeline records the reason).
    """
    try:
        for window in (ep_window, max_search):
            if rec.strand == "+":
                s = rec.stop + 1
                e = rec.stop + window
                gb = _efetch_gb(rec.accession, s, e)
                ep = str(gb.seq).upper()
            else:
                # downstream is toward lower coordinates; fetch below the aptamer and RC
                lo = max(1, rec.stop - window)  # rec.stop < rec.start on - strand
                hi = rec.stop - 1
                if hi < lo:
                    return Annotation(error="no downstream room on - strand")
                gb = _efetch_gb(rec.accession, lo, hi)
                ep = str(gb.seq.reverse_complement()).upper()

            # Unlike Prepping_sequences.ipynb (which drops N-containing seqs before
            # oligo synthesis), we KEEP them for classification: an N in the EP must
            # not discard the downstream-gene annotation. It only makes the folding
            # (terminator/SD) less reliable, which we flag.
            notes = []
            if "N" in ep:
                notes.append("EP contains N (mechanism may be unreliable)")

            product = gene = ""
            nts = ""
            truncated = False
            # On - strand the RC flips feature orientation; find CDS on the strand that
            # points away from the aptamer (the downstream gene). We take the first CDS.
            for feat in gb.features:
                if feat.type != "CDS":
                    continue
                if rec.strand == "+":
                    offset = int(feat.location.start)  # 0-based in the downstream window
                    # CDS end reaches the fetched-slice boundary => runs past the window.
                    reaches_edge = int(feat.location.end) >= len(gb.seq)
                else:
                    offset = len(gb.seq) - int(feat.location.end)
                    reaches_edge = int(feat.location.start) <= 0
                if offset < 0:
                    continue
                product = (feat.qualifiers.get("product", [""])[0]) or ""
                gene = (feat.qualifiers.get("gene", [""])[0]) or ""
                nts = offset
                # Partial-feature markers ('<'/'>') also indicate a run-off CDS.
                partial = ("<" in str(feat.location)) or (">" in str(feat.location))
                truncated = bool(reaches_edge or partial)
                break

            time.sleep(delay)
            # Return early only for a CDS that is fully closed within this window.
            # A truncated CDS at the small window falls through to max_search to try
            # to close it; if still open there, we flag it for extension.
            if (product and not truncated) or window == max_search:
                if not product:
                    notes.append("no CDS within search window")
                elif truncated:
                    notes.append(f"downstream CDS '{product}' start found but not closed "
                                 f"within {window} nt (candidate for window extension)")
                return Annotation(
                    ep_seq=ep, downstream_gene_product=product,
                    downstream_gene=gene, nts_to_start_codon=nts,
                    cds_extends_beyond_window=bool(product and truncated),
                    error="; ".join(notes),
                )
        return Annotation(error="unreachable")
    except EndOfSequence:
        return Annotation(error="aptamer at 3' end of accession; no downstream sequence deposited")
    except Exception as e:  # noqa: BLE001
        return Annotation(error=f"{type(e).__name__}: {e}")
