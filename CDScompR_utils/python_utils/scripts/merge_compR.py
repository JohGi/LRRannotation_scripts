#!/usr/bin/env python3
import sys
import os
import argparse
import pandas as pd

sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
)
from CDScompR_lib.plots import plot_identity_hist_from_csv


def extract_gff_info(gff_path, attributes_to_extract):
    """
    Fonction pour extraire dynamiquement les attributs d'un fichier GFF à partir d'une liste.

    Args:
        gff_path (str): Chemin vers le fichier GFF.
        attributes_to_extract (list): Liste des attributs à extraire (en plus de l'ID, obligatoire).

    Returns:
        df: DataFrame contenant l'ID du gène et les attributs demandés.
    """
    gff = pd.read_csv(
        gff_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seqid",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )

    genes_info = {}

    for _, row in gff.iterrows():
        if row["feature"] != "gene":
            continue

        attributes_field = row["attributes"]
        gene_id = None
        extracted = {}

        # Dictionnaire des attributs plats (ID=..., Name=..., etc.)
        flat_attributes = {}
        for attr in attributes_field.split(";"):
            if "=" in attr:
                key, val = attr.strip().split("=", 1)
                flat_attributes[key] = val

        # ID est obligatoire
        gene_id = flat_attributes.get("ID")
        if not gene_id:
            continue

        if attributes_to_extract:
            # Commentaires supplémentaires à parser séparément
            comment_block = flat_attributes.get("comment", "")
            comment_attributes = {}
            if comment_block:
                for part in comment_block.split("/"):
                    part = part.strip()
                    if ":" in part:
                        key, val = part.split(":", 1)
                        if key not in comment_attributes:
                            comment_attributes[key.strip()] = val.strip()

            # Extraire dynamiquement les attributs demandés
            for attr in attributes_to_extract:
                value = flat_attributes.get(attr)
                if value is None:
                    value = comment_attributes.get(attr)
                extracted[attr] = value

            # Stocker l'entrée
            genes_info[gene_id] = extracted
        else:
            genes_info[gene_id] = {"_keep": True}

    df = pd.DataFrame.from_dict(genes_info, orient="index").reset_index()
    df.rename(columns={"index": "gene_id"}, inplace=True)
    if "_keep" in df.columns:
        df.drop(columns=["_keep"], inplace=True)

    return df


# def extract_gff_info(gff_path):
#     """
#     Fonction pour extraire les informations des gènes à partir d'un fichier GFF.

#     Args:
#         gff_path (str): Chemin vers le fichier GFF.

#     Returns:
#         df: dataframe contenant les informations des gènes du gff.
#     """
#     # Charger le fichier GFF en tant que fichier tabulé, et ignorer les commentaires (#)
#     gff = pd.read_csv(gff_path, sep='\t', comment='#', header=None, names=[
#         'seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'
#     ])

#     # Dictionnaires pour stocker les informations des gènes parents et des CDS
#     genes_info = {}
#     cds_info = []

#     # Parcourir les lignes du fichier GFF pour extraire les gènes
#     for index, row in gff.iterrows():
#         source = row['source']
#         attributes = row['attributes']

#         # Si c'est un gène, extraire l'ID, l'origine et la gene_class
#         if row['feature'] == 'gene':
#             gene_id = None
#             origin = None
#             gene_class = None

#             # Parcourir les attributs du gène (colonne 9)
#             for attr in attributes.split(';'):
#                 if attr.startswith('ID='):
#                     gene_id = attr.split('=')[1]
#                 elif attr.startswith('comment='):
#                     # Enlever "comment=" et spliter sur "/"
#                     comment_value = attr.replace('comment=', '')
#                     comment_parts = comment_value.split('/')

#                     # Parcourir les segments du commentaire
#                     for part in comment_parts:
#                         part = part.strip()  # Supprimer les espaces en trop
#                         if part.startswith('Origin:') and origin is None: #je prends la première occurrence car parfois il y en a 2
#                             origin = part.split('Origin:')[1].strip()
#                         elif part.startswith('Gene-Class:'):
#                             gene_class = part.split('Gene-Class:')[1].strip()


#             # Enregistrer les informations du gène dans un dictionnaire
#             if gene_id:
#                 genes_info[gene_id] = {'source': source, 'origin': origin, 'gene_class': gene_class}


#     genes_info_df = pd.DataFrame.from_dict(genes_info, orient='index').reset_index()
#     genes_info_df.rename(columns={'index': 'gene_id'}, inplace=True)  # Renommer la colonne des IDs de gènes

#     return genes_info_df


def load_csv_files_from_list(csv_list_file, summarize, output_dir):
    """
    Fonction pour charger tous les fichiers CSV listés dans un fichier, extraire et renommer les colonnes,
    ajouter un suffixe à chaque fichier CSV pour différencier les colonnes, et plotter les histogrammes
    des scores.

    Args:
        csv_list_file (str): Chemin vers le fichier listant les chemins des CSV et suffixes séparés par des tabulations.
        summarize (bool): Si True, ne garde que les colonnes essentielles de la sortie de CDScompR (locus, identity score)
        output_dir (str): Dossier d'output
    Returns:
        list: Liste de DataFrames pandas des fichiers CSV modifiés avec les colonnes renommées.
    """

    cols_to_keep = ["Reference locus", "Alternative locus", "Identity score (%)"]
    new_col_names = ["gene_id", "Alt_locus", "Identity_Score"]
    if not summarize:
        cols_to_keep += [
            "Comparison matches",
            "Exon_Intron (EI) mismatches",
            "Reading Frame (RF) mismatches",
        ]
        new_col_names += ["Matches", "MismatchesEI", "MismatchesRF"]

    csv_files = []
    with open(csv_list_file, "r") as file:
        csv_paths = [line.strip().split() for line in file if line.strip()]

        for csv_path, suffix in csv_paths:
            print(f"\nLoading {csv_path}...")

            # Plot histogram
            os.makedirs(output_dir, exist_ok=True)
            out_path = os.path.join(output_dir, f"histogram_{suffix}.png")
            # plot_identity_hist_from_csv(df, cols_to_keep, suffix, output_dir)
            plot_identity_hist_from_csv(csv_path, "Ref gff", suffix, out_path)

            # Renommer les colonnes avec les nouveaux noms et ajouter le suffixe (sauf pour 'Reference locus')
            df = pd.read_csv(csv_path, usecols=cols_to_keep)
            renamed_cols = {
                old: (new if old == "Reference locus" else new + f"_{suffix}")
                for old, new in zip(cols_to_keep, new_col_names)
            }
            df.rename(columns=renamed_cols, inplace=True)
            csv_files.append(df)

    return csv_files


def main():
    # Initialiser le parser d'arguments
    parser = argparse.ArgumentParser(description="")

    # Ajouter les arguments pour les fichiers d'entrée
    parser.add_argument(
        "--comp_csv_list",
        required=True,
        help="Chemin vers le fichier listant les fichiers CSV de sortie de CDScompR à merger (1e colonne donne le csv, 2e colonne donne le suffixe correspondant)",
    )
    parser.add_argument(
        "--ref_gff",
        required=True,
        help="Chemin vers le fichier GFF utilisé comme référence lors des lancements de CDScompR",
    )
    parser.add_argument(
        "--output", required=True, help="Chemin pour le fichier de sortie (TSV)"
    )
    parser.add_argument(
        "--summarize",
        action="store_true",
        help="Si utilisé, ne garde que les colonnes essentielles de la sortie de CDScompR (locus, identity score)",
    )
    parser.add_argument(
        "--ref_gff_attributes",
        nargs="*",
        default=[],
        help="Liste des attributs à extraire du champ 9 du GFF en plus de ID (ex: --ref_gff_attributes Origin Gene-Class)",
    )

    # Récupérer les arguments
    args = parser.parse_args()
    output_dir = os.path.dirname(os.path.abspath(args.output))

    # Charger les fichiers CSV à partir de la liste
    comp_csv_files = load_csv_files_from_list(
        args.comp_csv_list, args.summarize, output_dir
    )

    # Extraire les informations du fichier GFF en appelant la fonction dédiée
    genes_info = extract_gff_info(args.ref_gff, args.ref_gff_attributes)
    print(
        f"\nNumber of genes extracted from the reference gff ({args.ref_gff}): {genes_info.shape[0]}"
    )

    # Fusionner chaque fichier CSV avec les informations du GFF
    merged_df = genes_info
    for i, csv_df in enumerate(comp_csv_files):
        merged_df = pd.merge(merged_df, csv_df, on="gene_id", how="left")

    # Trier le DataFrame par ordre alphabétique selon la colonne 'gene_id'
    merged_df.sort_values(by="gene_id", ascending=True, inplace=True)

    # Sauvegarder le résultat final dans le fichier de sortie
    merged_df.to_csv(args.output, sep="\t", index=False, na_rep="NA")

    # Afficher un message de succès
    print(f"\nMerged results saved in: {args.output}")


if __name__ == "__main__":
    main()
