#!/usr/bin/env python3
"""CLI: build a per-sequence classification CSV for one Rfam riboswitch family.

Usage
-----
Offline (classify a local Rfam FASTA/Stockholm you already downloaded):
    python run_classification.py --rfam RF00059 --in data/RF00059_seed.fa \
        --email you@example.com --out tpp_RF00059_classification.csv

Fetch the seed alignment from Rfam, then annotate via NCBI (RATE-LIMITED PULL —
authorise first; use --limit for a small pilot):
    python run_classification.py --rfam RF00059 --fetch-seed \
        --email you@example.com --limit 5 --out pilot.csv

The NCBI extension/annotation step is the only rate-limited part; it is bounded by
--limit and requires --email (NCBI policy). Every dropped/failed sequence is still
written to the CSV with mechanism_call='other' and a reason in `notes` (no silent loss).
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd
from Bio import Entrez

from classify import schema, classify_sequence
from classify import ncbi


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rfam", required=True, help="Rfam accession, e.g. RF00059")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--in", dest="infile", help="local Rfam FASTA or Stockholm file")
    src.add_argument("--fetch-seed", action="store_true", help="download the Rfam SEED alignment first")
    p.add_argument("--email", required=True, help="email for NCBI Entrez (required by NCBI)")
    p.add_argument("--api-key", default=None, help="optional NCBI API key (raises rate limit)")
    p.add_argument("--limit", type=int, default=0, help="max sequences to annotate (0 = all). Use a small value for a pilot.")
    p.add_argument("--ep-window", type=int, default=150, help="downstream nt to fetch as the EP")
    p.add_argument("--dg-threshold", type=float, default=-8.0, help="max terminator hairpin ΔG (kcal/mol)")
    p.add_argument("--out", required=True, help="output CSV path")
    args = p.parse_args(argv)

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    # 1) Load Rfam records (no network unless --fetch-seed).
    if args.fetch_seed:
        dest = f"data/{args.rfam}_seed.sto"
        print(f"[rfam] downloading SEED alignment for {args.rfam} -> {dest}", file=sys.stderr)
        ncbi.fetch_rfam_seed_stockholm(args.rfam, dest)
        records = ncbi.stockholm_to_records(dest)
    else:
        if args.infile.endswith((".sto", ".stk", ".stockholm")):
            records = ncbi.stockholm_to_records(args.infile)
        else:
            records = ncbi.parse_rfam_fasta(args.infile)
    print(f"[rfam] {len(records)} sequences parsed for {args.rfam}", file=sys.stderr)

    if args.limit:
        records = records[: args.limit]
        print(f"[limit] annotating first {len(records)} sequences (pilot)", file=sys.stderr)

    # 2) Annotate (RATE-LIMITED NCBI PULL) + classify each record.
    rows = []
    for i, rec in enumerate(records, 1):
        ann = ncbi.extend_and_annotate(rec, ep_window=args.ep_window)
        row = classify_sequence(
            rfam_family=args.rfam,
            aptamer_id=rec.aptamer_id,
            accession=rec.accession,
            genome_coords=rec.genome_coords,
            strand=rec.strand,
            aptamer_seq=rec.aptamer_seq,
            ep_seq=ann.ep_seq,
            downstream_gene_product=ann.downstream_gene_product,
            downstream_gene=ann.downstream_gene,
            nts_to_start_codon=ann.nts_to_start_codon,
            cds_extends_beyond_window=ann.cds_extends_beyond_window,
            dg_threshold=args.dg_threshold,
        )
        if ann.error:
            row["notes"] = "; ".join(x for x in [row["notes"], f"fetch: {ann.error}"] if x)
        rows.append(row)
        if i % 10 == 0 or i == len(records):
            print(f"[classify] {i}/{len(records)}", file=sys.stderr)

    # 3) Write CSV in canonical schema order.
    df = pd.DataFrame(rows, columns=schema.COLUMNS)
    df.to_csv(args.out, index=False)
    print(f"[done] wrote {len(df)} rows -> {args.out}", file=sys.stderr)

    # Console summary
    print("\n=== summary ===")
    print("mechanism_call:\n", df["mechanism_call"].value_counts().to_string())
    print("predicted_direction:\n", df["predicted_direction"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
