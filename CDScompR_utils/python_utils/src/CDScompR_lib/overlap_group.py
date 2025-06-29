from typing import List, Dict, Tuple
from attrs import define, field
from intervaltree import Interval, IntervalTree
from collections import defaultdict
from .gene import Gene


@define
class OverlapGroup:
    ref_genes: List["Gene"] = field(factory=list)
    pred_genes: List["Gene"] = field(factory=list)

    def get_type(self) -> str:
        """
        Infer the type of annotation change based on gene counts in the group.
        """
        n_ref = len(self.ref_genes)
        n_pred = len(self.pred_genes)

        if n_ref == 0 and n_pred > 0:
            return "appearance"
        elif n_ref > 0 and n_pred == 0:
            return "disappearance"
        elif n_ref == 1 and n_pred == 1:
            return "match"
        elif n_ref == 1 and n_pred > 1:
            return "split"
        elif n_ref > 1 and n_pred == 1:
            return "fusion"
        else:
            return "complex"

    def summarize(self) -> Dict:
        """
        Build a dictionary summary of the group, including gene IDs, coordinates,
        CDS lengths and best hits.
        """
        return {
            "ref_gene_ids": [g.id for g in self.ref_genes],
            "pred_gene_ids": [g.id for g in self.pred_genes],
            "ref_span_coords": [f"{g.span_start}-{g.span_end}" for g in self.ref_genes],
            "pred_span_coords": [f"{g.span_start}-{g.span_end}" for g in self.pred_genes],
            "ref_cumul_cds_lengths": [g.protein.cds_length() for g in self.ref_genes],
            "pred_cumul_cds_lengths": [g.protein.cds_length() for g in self.pred_genes],
            "ref_cds_counts": [g.protein.cds_count() for g in self.ref_genes],
            "pred_cds_counts": [g.protein.cds_count() for g in self.pred_genes],
            "ref_best_hits": [
                f"{g.best_hit_id or '<NA>'}: {g.identity_score}" for g in self.ref_genes
            ],
            "pred_best_hits": [
                f"{g.best_hit_id or '<NA>'}: {g.identity_score}" for g in self.pred_genes
            ],
            "type": self.get_type(),
        }

    @staticmethod
    def _group_by_chrom_and_strand(genes: List["Gene"],) -> Dict[Tuple[str, str], List["Gene"]]:
        """
        Group genes by (chromosome, strand).
        """
        chrom_strand_dict = defaultdict(list)
        for gene in genes:
            chrom = gene.protein.feature.seqid
            strand = gene.protein.feature.strand
            chrom_strand_dict[(chrom, strand)].append(gene)
        return chrom_strand_dict

    @staticmethod
    def _find_root(id: str, parents: Dict[str, str]) -> str:
        """
        Find the root of a node in a union-find structure, with path compression.

        Parameters:
            id: Gene unique identifier.
            parents: Disjoint-set forest mapping.

        Returns:
            The root ID for the given node.
        """
        if id not in parents:
            parents[id] = id
        while parents[id] != id:
            parents[id] = parents[parents[id]]
            id = parents[id]
        return id

    @staticmethod
    def _unite_roots(id1: str, id2: str, parents: Dict[str, str]) -> None:
        """
        Merge two disjoint sets in a union-find structure.

        Parameters:
            id1, id2: Gene identifiers to merge.
            parents: Disjoint-set forest mapping.
        """
        root1 = OverlapGroup._find_root(id1, parents)
        root2 = OverlapGroup._find_root(id2, parents)
        if root1 != root2:
            parents[root2] = root1

    @staticmethod
    def _build_tree(ref_genes: List["Gene"], pred_genes: List["Gene"]) -> IntervalTree:
        """
        Create an interval tree from reference and predicted genes.
        """
        tree = IntervalTree()
        for gene in ref_genes + pred_genes:
            start, end = sorted((gene.span_start, gene.span_end))
            tree.add(Interval(start, end, gene))
        return tree

    @staticmethod
    def _build_parents(tree: IntervalTree) -> Dict[str, str]:
        """
        Construct union-find structure to track connected overlapping genes.

        Returns:
            A mapping of each gene ID to its root in the disjoint-set forest.
        """
        parents: Dict[str, str] = {}

        for interval in tree:
            for overlap in tree.overlap(interval.begin, interval.end):
                if interval.data.uid != overlap.data.uid:
                    OverlapGroup._unite_roots(
                        interval.data.uid, overlap.data.uid, parents
                    )

        return parents

    @staticmethod
    def _build_groups_from_parents(
        tree: IntervalTree, parents: Dict[str, str]
    ) -> List["OverlapGroup"]:
        """
        Build a list of OverlapGroup objects from connected genes in the parents disjoint-set structure.
        """
        groups: Dict[str, OverlapGroup] = {}
        for interval in tree:
            group_id = OverlapGroup._find_root(interval.data.uid, parents)
            group = groups.setdefault(group_id, OverlapGroup())
            if interval.data.is_ref:
                group.ref_genes.append(interval.data)
            else:
                group.pred_genes.append(interval.data)
        return list(groups.values())

    @staticmethod
    def _build_groups(tree: IntervalTree) -> List["OverlapGroup"]:
        """
        Build a list of OverlapGroup objects from overlapping genes in the interval tree.
        """
        parents = OverlapGroup._build_parents(tree)
        return OverlapGroup._build_groups_from_parents(tree, parents)

    @staticmethod
    def overlap_groups_from_genes(ref_genes: List["Gene"], pred_genes: List["Gene"]) -> List["OverlapGroup"]:
        """
        Identify groups of overlapping genes between two annotations.

        Parameters:
            ref_genes: List of reference Gene objects.
            pred_genes: List of predicted Gene objects.

        Returns:
            A list of OverlapGroup instances.
        """
        grouped_ref_genes = OverlapGroup._group_by_chrom_and_strand(ref_genes)
        grouped_pred_genes = OverlapGroup._group_by_chrom_and_strand(pred_genes)
        all_groups = []

        for chrom_strand in set(grouped_ref_genes.keys()) | set(grouped_pred_genes.keys()):
            interval_tree = OverlapGroup._build_tree(
                grouped_ref_genes.get(chrom_strand, []), grouped_pred_genes.get(chrom_strand, [])
            )
            chrom_groups = OverlapGroup._build_groups(interval_tree)
            all_groups.extend(chrom_groups)

        return all_groups
