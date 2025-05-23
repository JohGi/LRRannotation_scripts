import pandas as pd
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

    cols_to_keep = ["Reference locus", "Alternative locus", "Identity score (%)"]
    df = pd.read_csv(csv, usecols=cols_to_keep)
    ref_col, alt_col, score_col = cols_to_keep

    # Nombre total de gènes dans le GFF alternatif (≠ "~" dans alt)
    alt_gff_nb_genes = df[alt_col][df[alt_col] != "~"].nunique()
    ref_gff_nb_genes = df[ref_col][df[ref_col] != "~"].nunique()

    # Filtrage : valeurs sans "~" dans alt et ref
    filtered_df = df[(df[alt_col] != "~") & (df[ref_col] != "~")]

    # Convertir score en numérique
    values = pd.to_numeric(filtered_df[score_col], errors="coerce").dropna()
    nb_common_genes = values.shape[0]
    nb_genes_ref_only = ref_gff_nb_genes - nb_common_genes
    nb_genes_alt_only = alt_gff_nb_genes - nb_common_genes

    # Print info
    print(f"\nComparison of {ref_name} with {alt_name}:")
    print(f"--- Common genes: {nb_common_genes}")
    print(f"--- Genes only in {ref_name}: {nb_genes_ref_only}")
    print(f"--- Genes only in {alt_name}: {nb_genes_alt_only}")

    # Tracer

    plt.figure()
    plt.hist(values, bins=30, edgecolor="black")

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
        f"[{ref_name}\u2192( {nb_genes_ref_only} ( {nb_common_genes} ) {nb_genes_alt_only} )\u2190{alt_name}]",
        fontsize=10,
        ha="center",
        transform=plt.gca().transAxes,
    )

    plt.xlabel("Identity Score (%)")
    plt.ylabel("Fréquence")
    plt.xlim(-5, 105)
    plt.xticks(np.arange(0, 101, 10), rotation=45)
    plt.tight_layout()

    plt.savefig(out_path)
    plt.close()

    print(f"[{ref_name} vs. {alt_name}] Histogramme sauvegardé : {out_path}")
