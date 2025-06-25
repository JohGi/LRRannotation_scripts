#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare expert and predicted GFF annotations by overlapping genes and summarize CDS stats.
Adds best hit info and identity score from CDScompR CSV output.
"""
import argparse
import pandas as pd
import sys, os
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
)
from CDScompR_lib.gene import Gene
from CDScompR_lib.overlap_group import OverlapGroup
from CDScompR_lib.comparison_utils import add_identity_scores, summarize_overlaps, load_score_file
from CDScompR_lib.gff_utils import build_db

def main():
    parser = argparse.ArgumentParser(description="Compare expert and predicted GFF annotations.")
    parser.add_argument("--ref_gff", help="Reference (expert) GFF file")
    parser.add_argument("--pred_gff", help="Predicted GFF file")
    parser.add_argument("--cdscompr_csv", help="CDScompR csv output file")
    parser.add_argument("--span_type", choices=["gene", "mRNA", "CDS"], default="gene",
                        help="Feature span to use for overlap detection (default: gene)")
    parser.add_argument("-o", "--output", help="Output TSV file")

    args = parser.parse_args()

    print("Building GFF databases...")
    ref_genes_db = build_db(args.ref_gff)
    pred_genes_db = build_db(args.pred_gff)

    print("Parsing genes...")
    ref_genes = [Gene.from_gff(ref_genes_db, g, is_ref=True, span_type=args.span_type) for g in ref_genes_db.features_of_type("gene")]
    pred_genes = [Gene.from_gff(pred_genes_db, g, is_ref=False, span_type=args.span_type) for g in pred_genes_db.features_of_type("gene")]


    if args.cdscompr_csv:
        score_df = load_score_file(args.cdscompr_csv)
        add_identity_scores(ref_genes, score_df, is_ref=True)
        add_identity_scores(pred_genes, score_df, is_ref=False)

    print("Detecting overlapping gene groups...")
    overlap_groups = OverlapGroup.overlap_groups_from_genes(ref_genes, pred_genes)
    print(f"Found {len(overlap_groups)} overlapping groups.")

    summarize_overlaps(overlap_groups, args.span_type, args.output)

if __name__ == "__main__":
    main()
