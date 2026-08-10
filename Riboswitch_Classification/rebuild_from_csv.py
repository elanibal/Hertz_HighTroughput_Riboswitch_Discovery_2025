#!/usr/bin/env python3
"""Re-classify an existing classification CSV with the current classifiers.

Runs the pure classifiers OFFLINE on the sequences already stored in the CSV, and
re-fetches ONLY rows whose EP is empty (previous fetch failures) and whose accession
is a real NCBI nucleotide accession (skips RNAcentral URS ids, which NCBI can't serve).
This lets us iterate on the classifiers without re-pulling the whole family.

Usage: python rebuild_from_csv.py in.csv out.csv --email you@example.com
"""
from __future__ import annotations

import argparse
import re
import sys

import pandas as pd
from Bio import Entrez

from classify import schema, classify_sequence
from classify import ncbi


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--email", required=True)
    ap.add_argument("--seed", default=None, help="Stockholm seed to source coords for re-fetch")
    args = ap.parse_args(argv)
    Entrez.email = args.email

    df = pd.read_csv(args.infile).fillna("")
    seed_by_id = {}
    if args.seed:
        seed_by_id = {r.aptamer_id: r for r in ncbi.stockholm_to_records(args.seed)}

    refetched = 0
    rows = []
    for _, r in df.iterrows():
        ep = str(r["ep_seq"])
        product = str(r["downstream_gene_product"])
        gene = str(r["downstream_gene"])
        nts = r["nts_to_start_codon"]
        acc = str(r["accession"])
        note_extra = ""
        # Preserve the truncation flag already in the CSV (a fetch-time observation).
        trunc = str(r.get("cds_extends_beyond_window", "")).lower() == "true"

        # Re-fetch prior failures on real NCBI accessions: rows with no EP OR with an
        # EP but no product (recovers the fixed N-in-EP discards).
        if (not ep or not product) and (not acc.startswith("URS")) and (r["aptamer_id"] in seed_by_id):
            ann = ncbi.extend_and_annotate(seed_by_id[r["aptamer_id"]])
            if ann.ep_seq:
                ep, product, gene, nts = ann.ep_seq, ann.downstream_gene_product, ann.downstream_gene, ann.nts_to_start_codon
                trunc = ann.cds_extends_beyond_window
                refetched += 1
                print(f"[refetch OK] {r['aptamer_id']} -> ep={len(ep)}nt product={product!r}", file=sys.stderr)
            if ann.error:
                note_extra = f"fetch: {ann.error}"
        elif (not ep) and acc.startswith("URS"):
            note_extra = "fetch: RNAcentral URS id (no NCBI nucleotide record; cannot extend)"
        row = classify_sequence(
            rfam_family=str(r["rfam_family"]), aptamer_id=str(r["aptamer_id"]),
            accession=acc, genome_coords=str(r["genome_coords"]), strand=str(r["strand"]),
            aptamer_seq=str(r["aptamer_seq"]), ep_seq=ep,
            downstream_gene_product=product, downstream_gene=gene, nts_to_start_codon=nts,
            cds_extends_beyond_window=trunc,
        )
        if note_extra:
            row["notes"] = "; ".join(x for x in [row["notes"], note_extra] if x)
        rows.append(row)

    out = pd.DataFrame(rows, columns=schema.COLUMNS)
    out.to_csv(args.outfile, index=False)
    print(f"[done] re-fetched {refetched} rows; wrote {len(out)} -> {args.outfile}", file=sys.stderr)
    print("\nmechanism_call:\n", out["mechanism_call"].value_counts().to_string())
    print("\npredicted_direction:\n", out["predicted_direction"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
