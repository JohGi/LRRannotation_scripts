import pandas as pd
from typing import List, Optional
from .gene import Gene
from .overlap_group import OverlapGroup

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
