import os
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
    count = values.shape[0]

    # Print info
    print(f"\nComparison of {ref_name} with {alt_name}:")
    print(f"--- Common genes: {count}")
    print(f"--- Genes only in {ref_name}: {ref_gff_nb_genes - count}")
    print(f"--- Genes only in {alt_name}: {alt_gff_nb_genes - count}")

    # Tracer
    plt.figure()
    plt.hist(values, bins=30, edgecolor="black")
    plt.title(
        f"{ref_name} vs. {alt_name}:\nID scores of their {count} overlapping genes"
    )
    plt.xlabel("Identity Score (%)")
    plt.ylabel("Fréquence")
    plt.xlim(-5, 105)
    plt.xticks(np.arange(0, 101, 10), rotation=45)
    plt.tight_layout()

    plt.savefig(out_path)
    plt.close()

    print(f"[{ref_name} vs. {alt_name}] Histogramme sauvegardé : {out_path}")
