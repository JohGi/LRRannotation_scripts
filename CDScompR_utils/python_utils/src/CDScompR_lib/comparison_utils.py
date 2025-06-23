import pandas as pd
from typing import List, Optional
from .gene import Gene
from .overlap_group import OverlapGroup

def load_score_file(csv_path: str) -> pd.DataFrame:
    """
    Load and preprocess a CDScompR CSV score file.
    Keeps only relevant columns and renames them for easier downstream use.
    """
    print("Loading score CSV...")

    cols_needed = [
        "Reference locus",
        "Alternative locus",
        "Identity score (%)"
    ]
    dtypes = {
        "Reference locus": "string",
        "Alternative locus": "string",
        "Identity score (%)": "float"
    }

    score_df = pd.read_csv(
        csv_path,
        usecols=cols_needed,
        na_values=["_", "~"],
        dtype=dtypes,
        low_memory=False
    )

    score_df.rename(columns={
        "Reference locus": "ref_id",
        "Alternative locus": "alt_id",
        "Identity score (%)": "identity_score"
    }, inplace=True)

    return score_df

def add_identity_scores(genes: List[Gene], score_df: pd.DataFrame, is_ref: bool) -> None:
    """
    Update all Gene objects in the list with best hit ID and identity score.
    """
    for gene in genes:
        gene.set_identity_scores(score_df, is_ref)


def summarize_overlaps(groups: List[OverlapGroup], output_path: Optional[str] = None) -> None:
    """
    Summarize a list of OverlapGroups into a DataFrame and write it to a TSV file or print it.

    """
    summary_rows = [group.summarize() for group in groups]
    df = pd.DataFrame(summary_rows)

    if output_path:
        df.to_csv(output_path, sep="\t", index=False)
        print(f"Results written to {output_path}")
    else:
        print(df)
