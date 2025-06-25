import polars as pl
import matplotlib.pyplot as plt
import numpy as np


def plot_identity_hist_from_csv(csv, ref_name, alt_name, out_path):
    """
    Plot a histogram of the identity scores computed by CDScompR.

    Args:
        csv: CDScompR csv output
        ref_name (str): Name of the ref gff (for printing)
        alt_name (str): Name of the alt gff (for printing)
        out_path (str): Output file path
    """

    cols = ["Reference locus", "Alternative locus", "Identity score (%)"]
    df = pl.read_csv(csv, null_values=["~"], columns=cols)
    ref_col, alt_col, score_col = cols

    # Count the number of unique non-null genes in ref and alt
    alt_gff_nb_genes = df.filter(pl.col(alt_col).is_not_null()).select(pl.col(alt_col).unique()).height
    ref_gff_nb_genes = df.filter(pl.col(ref_col).is_not_null()).select(pl.col(ref_col).unique()).height

    # Keep only rows where both ref and alt are defined
    filtered = df.filter((pl.col(ref_col).is_not_null()) & (pl.col(alt_col).is_not_null()))

    # Extract identity scores as float values (drop NaN)
    values = filtered[score_col].cast(pl.Float64).drop_nulls()
    values_np = values.to_numpy()

    # Summary statistics
    nb_common_genes = len(values_np)
    nb_genes_ref_only = ref_gff_nb_genes - nb_common_genes
    nb_genes_alt_only = alt_gff_nb_genes - nb_common_genes
    mean_score = values_np.mean()

    # Print summary to terminal
    print(f"\nComparison of {ref_name} with {alt_name}:")
    print(f"--- Common genes: {nb_common_genes}")
    print(f"--- Genes only in {ref_name}: {nb_genes_ref_only}")
    print(f"--- Genes only in {alt_name}: {nb_genes_alt_only}")

    # Plot the histogram
    plt.figure()
    plt.hist(values_np, bins=30, edgecolor="black")

    # Add vertical dashed line at the mean score
    plt.axvline(mean_score, color="red", linestyle="--", linewidth=1.5, label="Mean score")

    # Annotate the mean value next to the line
    plt.text(mean_score + 1, plt.gca().get_ylim()[1] * 0.9, f"{mean_score:.1f}%", color="red")

    # Annotate with gene counts and title
    plt.text(
        0.5,
        1.1,
        f"{ref_name} vs. {alt_name}:\nID scores of their {nb_common_genes} overlapping genes",
        fontsize=12,
        ha="center",
        transform=plt.gca().transAxes,
    )
    plt.text(
        0.5,
        1.02,
        f"[{ref_name}→( {nb_genes_ref_only} ( {nb_common_genes} ) {nb_genes_alt_only} )←{alt_name}]",
        fontsize=10,
        ha="center",
        transform=plt.gca().transAxes,
    )

    plt.xlabel("Identity Score (%)")
    plt.ylabel("Frequency")
    plt.xlim(-5, 105)
    plt.xticks(np.arange(0, 101, 10), rotation=45)
    plt.legend()
    plt.tight_layout()

    plt.savefig(out_path)
    plt.close()

    print(f"[{ref_name} vs. {alt_name}] Histogram saved: {out_path}")
