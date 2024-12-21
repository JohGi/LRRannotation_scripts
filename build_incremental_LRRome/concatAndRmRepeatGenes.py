#python concatAndRmRepeatGenes.py --gff_list gff_list.txt --output concat_noDup.gff

import argparse
import os

def check_gene_in_same_gff(gene_ID, gff_gene_list, gff_file):
    if gene_ID in gff_gene_list:
        raise ValueError(f"Erreur : le gène {gene_ID} est présent plusieurs fois dans le fichier {gff_file}")
    gff_gene_list.append(gene_ID)


def append_gff(gff_file, gene_list, out_file):
    """
    Lit un gff ligne par ligne et les ajoute à out_file.
    Elle attend un gff avec les features dans le bon ordre (idéalement sorti de gff_cleaner).
    A chaque nouveau gène elle rajoute le gene_ID à gene_list.
    Si le gène est présent plusieurs fois dans le gff >> erreur.
    Si le gène a déjà été vu dans un gff précédent elle n'écrit pas ses annotations dans out_file.
        >> si on lui passe les gff du plus récent au plus ancien il supprime les doublons et ne garde que les versions les plus récentes dans out_file.
    """
    gff_gene_list=[]
    with open(gff_file, "r") as gff:
       for line in gff:
            seqname, source, feature, start, end, score, strand, frame, attribute = line.split("\t")
            if (feature == "gene"):
                gene_ID=attribute.split(";")[0].replace("ID=", "")
                check_gene_in_same_gff(gene_ID, gff_gene_list, gff_file)
                if gene_ID not in gene_list:
                    gene_list.append(gene_ID)
                    out_file.write(line)
                    skip=False
                else:
                    print(f"INFO: Les annotations du gène {gene_ID} du fichier {gff_file} n'ont pas été conservées car elles étaient présentes dans un gff plus récent.\n")
                    skip=True
            elif (skip==False):
                out_file.write(line)


def main():
    parser = argparse.ArgumentParser(description="Concatène une série de gff. Si le même gène apparaît plusieurs fois dans des gff différents c'est l'annotation issue du gff le plus récent qui conservée. A utiliser après gff_cleaner (les features doivent apparaître dans le bon ordre). Une erreur est levée si un même gff contient plusieurs fois le même gène.")

    parser.add_argument("--gff_list", help="Fichier listant les gff du plus récent au plus ancien. Soit les chemins complets (et laisser préfixe vide), soit juste les noms de fichiers et on peut passer un chemin +et/ou préfixe à '--prefix'. Si on passe juste un chemin il faut inclure le slash final.")
    parser.add_argument("--prefix", help="Préfixe des fichiers ggf. Ex : /home/user/02_build_exp_LRRome/CLEANED_GFF/cleaned_", default="")
    parser.add_argument("-o", "--output", help="Nom du fichier de sortie")

    args = parser.parse_args()
    outfile_name=args.output
    gfflist_name=args.gff_list
    prefix=args.prefix

    gene_list=[]
    with open(gfflist_name, "r") as gff_list:
        with open(outfile_name, "w") as out_file:
            for gff_file in gff_list:
                gff_file = gff_file.strip()
                gff_file = os.path.realpath(prefix + gff_file)
                append_gff(gff_file, gene_list, out_file)


if __name__ == "__main__":
    main()