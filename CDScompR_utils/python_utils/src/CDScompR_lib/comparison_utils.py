import pandas as pd
import polars as pl
from typing import List, Optional
from .gene import Gene
from .overlap_group import OverlapGroup


def load_score_file(csv_path: str) -> pl.DataFrame:
    """
    Load and preprocess a CDScompR CSV score file.
    Keeps only relevant columns and renames them for easier downstream use.
    """
    print("Loading score CSV...")

    needed_cols = ["Reference locus", "Alternative locus", "Identity score (%)"]

    score_df = pl.read_csv(
        csv_path,
        columns=needed_cols,
        null_values=["_", "~"],
        try_parse_dates=False,
    ).rename({
        "Reference locus": "ref_id",
        "Alternative locus": "alt_id",
        "Identity score (%)": "identity_score",
    })

    return score_df


def add_identity_scores(genes: List[Gene], score_df: pl.DataFrame, is_ref: bool) -> None:
    """
    Update all Gene objects in the list with best hit ID and identity score.
    """

    self_col = "ref_id" if is_ref else "alt_id"
    hit_col = "alt_id" if is_ref else "ref_id"

    best_hit_lookup = {
        row[self_col]: (row[hit_col], row["identity_score"])
        for row in score_df.iter_rows(named=True)
    }
    for gene in genes:
        gene.set_identity_scores(best_hit_lookup)


def summarize_overlaps(groups: List[OverlapGroup], span_type: str, output_path: Optional[str] = None) -> None:
    """
    Summarize a list of OverlapGroups into a DataFrame and write it to a TSV file or print it.

    """
    summary_rows = [group.summarize() for group in groups]
    df = pd.DataFrame(summary_rows)
    
    df.insert(0, "span_type", span_type)

    if output_path:
        df.to_csv(output_path, sep="\t", index=False)
        print(f"Results written to {output_path}")
    else:
        print(df)
