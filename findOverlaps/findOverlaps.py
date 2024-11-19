import pandas as pd
import argparse
import os


### ------------ Classes ---------------- ###

class Gene:
    instances = []

    def __init__(self, gene_id, start, end, dna_mol, gene_nature=None, gene_source=None):
        self.gene_id = gene_id
        self.start = start
        self.end = end
        self.dna_mol = dna_mol
        self.gene_nature = gene_nature
        self.gene_source = gene_source
        self.overlapping_genes = []  # List of overlapping genes in the other GFF
        self.num_overlapping_genes = 0  # Count of overlapping genes
        Gene.instances.append(self)

    def add_overlap(self, other_gene):
        self.overlapping_genes.append(other_gene)
        self.num_overlapping_genes += 1

    def get_overlapping_gene_ids(self):
        """Returns a human-readable list of the overlapped genes IDs."""
        return [gene.gene_id for gene in self.overlapping_genes]

    def display(self):
        """Shows all attributes of the gene."""
        print(f"Gene ID: {self.gene_id}")
        print(f"Start: {self.start}")
        print(f"End: {self.end}")
        print(f"DNA Mol: {self.dna_mol}")
        print(f"Gene Nature: {self.gene_nature}")
        print(f"Gene Source: {self.gene_source}")
        print(f"Overlapping Genes: {self.get_overlapping_gene_ids()}")
        print(f"Number of Overlapping Genes: {self.num_overlapping_genes}")
        print()

    @classmethod
    def find_GeneOverlapGroups(cls):
        """Searches for all instances with more than one overlapping gene and instantiates a GeneOverlapGroup for each, first sorting the overlapping_genes by start coordinate."""
        for gene in cls.instances:
            if gene.num_overlapping_genes > 1:
                sorted_overlapping_genes = sorted(gene.overlapping_genes, key=lambda x: x.start)
                overlap_group = GeneOverlapGroup(gene, sorted_overlapping_genes)


class GeneOverlapGroup:
    instances = []

    def __init__(self, main_gene, overlapping_genes):
        self.main_gene = main_gene
        self.overlapping_genes = overlapping_genes
        GeneOverlapGroup.instances.append(self)

    def display(self):
        """Shows all attributes of the group."""
        print(f"Main Gene ID: {self.main_gene.gene_id}")
        print(f"Overlapping Gene IDs: {[gene.gene_id for gene in self.overlapping_genes]}")
        print()

    @classmethod
    def display_all(cls):
        """Display all GeneOverlapGroup instances"""
        for group in cls.instances:
            group.display()

    def encompasses_all_overlaps(self, overreach_threshold):
        """
        Determines whether no overlapping gene exceeds the boundaries of the main gene by a length greater than the provided overreach_threshold.
        Returns True if, for each overlapping gene, both at the start and the end of the main gene:
        - either there's no overreach
        - or the overreach_length/overlapping_gene_length <= overreach_threshold
        Returns False otherwise.
        """
        # For the start we only need to test the first gene because the overlapping_genes are ordered by their start values
        first_gene = self.overlapping_genes[0]
        start_diff = first_gene.start - self.main_gene.start
        if start_diff < 0:
            length = first_gene.end - first_gene.start + 1
            start_ratio = (-1 * start_diff) / length
            if start_ratio > overreach_threshold:
                return False

        for overlapping_gene in self.overlapping_genes:
            end_diff = self.main_gene.end - overlapping_gene.end
            if end_diff < 0:
                length = overlapping_gene.end - overlapping_gene.start + 1
                end_ratio = (-1 * end_diff) / length
                if end_ratio > overreach_threshold:
                    return False

        return True

    def overlaps_dont_overlap(self, overlap_threshold):
        """
        Determines whether no pair of consecutive genes overlapping the main gene overlap with each other by a length greater than the provided overlap_threshold.
        Returns True if, for each pair:
        - either there's no overlap
        - or the overlap_length/main_gene_length <= overlap_threshold
        Returns False otherwise.
        """
        main_gene_length = self.main_gene.end - self.main_gene.start + 1

        for i in range(len(self.overlapping_genes) - 1):
            gene1 = self.overlapping_genes[i]
            gene2 = self.overlapping_genes[i + 1]
            diff = gene2.start - gene1.end
            if diff < 0:
                overlap_ratio = (-1 * diff) / main_gene_length
                if overlap_ratio > overlap_threshold:
                    return False

        return True


    def all_overlaps_have_single_overlap(self):
        """Returns True if all overlapping genes have only one overlap, False otherwise."""
        return all(gene.num_overlapping_genes == 1 for gene in self.overlapping_genes)


    @classmethod
    def find_valid_groups(cls, overreach_threshold, overlap_threshold):
        """
        Goes through all instances of GeneOverlapGroup and returns a dataframe of all the groups meeting the three conditions:
        - all_overlaps_have_single_overlap(),
        - encompasses_all_overlaps(overreach_threshold)
        - and overlaps_dont_overlap(overlap_threshold).
        """

        results = []

        for group in cls.instances:
            if (group.all_overlaps_have_single_overlap() and group.encompasses_all_overlaps(overreach_threshold) and group.overlaps_dont_overlap(overlap_threshold)):
                main_gene_id = group.main_gene.gene_id
                main_gene_source = group.main_gene.gene_source
                num_overlapping_genes = group.main_gene.num_overlapping_genes
                overlapping_genes_ids = ', '.join(gene.gene_id for gene in group.overlapping_genes)
                overlapping_genes_sources = ', '.join(gene.gene_source for gene in group.overlapping_genes)

                results.append({
                    "mainGene_id": main_gene_id,
                    "mainGene_source": main_gene_source,
                    "number_overlapping_genes": num_overlapping_genes,
                    "overlapping_genes_id": f"[{overlapping_genes_ids}]",
                    "overlapping_genes_source": f"[{overlapping_genes_sources}]"
                })

        return pd.DataFrame(results)



### ------------------------------------- ###



def read_gff(gff_path, source):
    """Reads a GFF file and returns a dictionary of genes by DNA molecule (= chromosome strand)."""
    gff_data = pd.read_csv(gff_path, sep='\t', comment='#', header=None)
    gff_dict = {}

    for index, row in gff_data.iterrows():
        seqid, _, feature, start, end, _, strand, _, attributes = row
        if feature == 'gene':
            gene_id = attributes.split('ID=')[1].split(';')[0]

            dna_mol = f"{seqid}_{strand}"

            # Extract extra info from the attributes field (ie the nature and confidence if they can be found)
            gene_nature = None
            gene_source = source
            for attr in attributes.split(';'):
                if attr.startswith('comment='):
                    comment_content = attr.replace('comment=', '')
                    for item in comment_content.split(' / '):
                        if item.startswith('Gene-Class:'):
                            gene_nature = item.split(':')[1].strip()
                if attr.startswith('confidence='):
                    confidence = attr.replace('confidence=', '')
                    gene_source = f"{gene_source}_{confidence}"

            # Instantiate the gene
            gene = Gene(gene_id, int(start), int(end), dna_mol, gene_nature, gene_source)

            # Update gff_dict
            if dna_mol not in gff_dict:
                gff_dict[dna_mol] = []
            gff_dict[dna_mol].append(gene)

    return gff_dict


def print_gff_dict(gff_dict):
    for dna_mol, genes in gff_dict.items():
        print(f"--- DNA Molecule: {dna_mol} ---")
        for gene in genes:
            gene.display()


def overlap(gene1, gene2):
    """Returns True if the two genes overlap."""
    return (gene2.start <= gene1.start <= gene2.end) or (gene1.start <= gene2.start <= gene1.end)


def find_overlaps(ref_gff_dict, alt_gff_dict, alt_prefix, verbose=False):
    """Detects all overlaps bewteen two sets of Gene instances (REF and ALT), each stored in a dictionary. Returns a dataframe of all the REF genes overlaps."""
    if verbose: print(f"\nDetecting overlaps...")

    results = []

    for dna_mol in ref_gff_dict.keys() | alt_gff_dict.keys():
        for ref_gene in ref_gff_dict.get(dna_mol, []):
            for alt_gene in alt_gff_dict.get(dna_mol, []):
                if overlap(ref_gene, alt_gene):
                    if verbose: print(f" --> {ref_gene.gene_id} and {alt_gene.gene_id} overlap")
                    # Update both genes to keep track of overlaps with the other gff
                    ref_gene.add_overlap(alt_gene)
                    alt_gene.add_overlap(ref_gene)

        # For each REF gene, list its overlapping genes, and how many genes each of them overlaps
        for ref_gene in ref_gff_dict.get(dna_mol, []):
            overlapping_genes=ref_gene.overlapping_genes
            num_overlaps_alt_genes=[gene.num_overlapping_genes for gene in overlapping_genes]
            overlapping_gene_ids = ', '.join(ref_gene.get_overlapping_gene_ids())
            results.append((
                ref_gene.gene_id,
                ref_gene.gene_nature,
                ref_gene.num_overlapping_genes,
                f"[{overlapping_gene_ids}]",
                num_overlaps_alt_genes
            ))


    results_df = pd.DataFrame(results, columns=['REF_gene_id', 'REF_gene_nature', f'REF_num_overlaps_{alt_prefix}', f'{alt_prefix}_overlaps_id', f'{alt_prefix}_num_overlaps'])

    # delete potential duplicated rows (it can happen when an input gff contains duplicated genes by mistake)
    results_df.drop_duplicates(subset=['REF_gene_id', 'REF_gene_nature'], inplace=True)

    return results_df



def main():

    ## RETRIEVE ARGUMENTS
    parser = argparse.ArgumentParser(description='Compare two GFF files to find overlaps between genes.')
    parser.add_argument('--ref_gff', required=True, help='Path to the reference GFF file.')
    parser.add_argument('--alt_gff', required=True, help='Path to the alternative GFF file.')
    parser.add_argument('--overlaps_output', required=True, help='Name of the output file for the overlaps results.')
    parser.add_argument('--groups_output', help='Name of the output file for the overlapping groups results (optional).')
    parser.add_argument('--overreach_thr', type=float, help='Overreach max threshold to keep a group (if you provided a --groups_output). Set to 0.05 for 5% of the main gene length.')
    parser.add_argument('--overlap_thr', type=float, help='Overlap max threshold to keep a group (if you provided a --groups_output). Set to 0.05 for 5% of the main gene length.')
    parser.add_argument('--prefix', default='', help='Prefix for the ALT gene id column name.')
    parser.add_argument('--verbose', action='store_true', help="Display more information.")
    parser.add_argument('--show_all_genes', action='store_true', help="Display information for all genes in both GFF files.")
    parser.add_argument('--show_all_groups', action='store_true', help="Display information for all overlapping groups detected between the two gff files.")
    args = parser.parse_args()
    alt_prefix = args.prefix if args.prefix else "ALT"


    # Check threshold arguments in case a groups_output was provided
    if args.groups_output:
        if args.overreach_thr is None or args.overlap_thr is None:
            parser.error("--overreach_thr and --overlap_thr are required when a --groups_output is provided.")


    ## READ INPUT GFF FILES
    ref_gff_dict = read_gff(args.ref_gff, "REF")
    alt_gff_dict = read_gff(args.alt_gff, "ALT")


    ## DETECT OVERLAPS
    overlaps_results = find_overlaps(ref_gff_dict, alt_gff_dict, alt_prefix, args.verbose)

    if args.show_all_genes:
        print(f"\n\n######## {os.path.basename(args.ref_gff)} CONTAINS THE FOLLOWING GENES: ########\n")
        print_gff_dict(ref_gff_dict)
        print(f"\n\n######## {os.path.basename(args.alt_gff)} CONTAINS THE FOLLOWING GENES: ########\n")
        print_gff_dict(alt_gff_dict)

    # write output file
    overlaps_results.to_csv(args.overlaps_output, sep='\t', index=False)
    print(f"\nResults saved to {args.overlaps_output}\n")


    ## DETECT OVERLAPPING GROUPS OF INTEREST
    if args.groups_output:
        # detect gene groups (= one main gene overlapping several genes)
        Gene.find_GeneOverlapGroups()

        # identify gene groups of interest
        valid_groups_results = GeneOverlapGroup.find_valid_groups(overreach_threshold=args.overreach_thr, overlap_threshold=args.overlap_thr)

        if args.show_all_groups:
            print(f"\n\n######## THE FOLLOWING OVERLAPPING GROUPS WERE DETECTED: ########\n")
            GeneOverlapGroup.display_all()

        # write output file
        valid_groups_results.to_csv(args.groups_output, sep='\t', index=False)
        print(f"\nResults saved to {args.groups_output}\n")


if __name__ == '__main__':
    main()
