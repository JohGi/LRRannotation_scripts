import pytest
import os, sys

script_dir = os.path.dirname(__file__)
script_dir = "/".join(script_dir.split("/")[:-1]) + "/"
sys.path.append(script_dir)

from findOverlaps import Gene, GeneOverlapGroup


# Fixture pour créer des groupes d'overlap dynamiquement
@pytest.fixture
def overlap_group_factory():
    def _create(main_gene, *overlapping_genes):
        for gene in overlapping_genes:
            main_gene.add_overlap(gene)
        return GeneOverlapGroup(main_gene, list(overlapping_genes))
    return _create

# Paramétrisation avec les coordonnées des gènes
@pytest.mark.parametrize(
    "main_gene_data, overlap_genes_data, threshold, expected_result", [
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt2", 550, 580, "chr1")], 0, True),
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt3", 550, 630, "chr1")], 0, False),
        (("ref1", 500, 600, "chr1"), [("alt4", 450, 530, "chr1"), ("alt2", 550, 580, "chr1")], 0, False),
        (("ref1", 500, 600, "chr1"), [("alt4", 450, 530, "chr1"), ("alt2", 550, 580, "chr1")], 0.61, False),
        (("ref1", 500, 600, "chr1"), [("alt4", 450, 530, "chr1"), ("alt2", 550, 580, "chr1")], 0.62, True),
        (("ref1", 500, 600, "chr1"), [("alt4", 450, 530, "chr1"), ("alt3", 550, 630, "chr1")], 0, False)
    ]
)
def test_encompasses_all_overlaps(main_gene_data, overlap_genes_data, threshold, expected_result, overlap_group_factory):
    main_gene = Gene(*main_gene_data)
    overlap_genes = [Gene(*gene_data) for gene_data in overlap_genes_data]

    overlap_group = overlap_group_factory(main_gene, *overlap_genes)

    result = overlap_group.encompasses_all_overlaps(threshold)
    assert result == expected_result


# Paramétrisation avec les coordonnées des gènes pour le test overlaps_dont_overlap
@pytest.mark.parametrize(
    "main_gene_data, overlap_genes_data, threshold, expected_result", [
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt2", 550, 580, "chr1")], 0, True),
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt5", 520, 550, "chr1")], 0, False),
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt5", 520, 550, "chr1")], 0.09, False),
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt5", 520, 550, "chr1")], 0.11, True),
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt6", 510, 530, "chr1")], 0, False),
        (("ref1", 500, 600, "chr1"), [("alt1", 510, 530, "chr1"), ("alt7", 515, 525, "chr1")], 0, False)
    ]
)
def test_overlaps_dont_overlap(main_gene_data, overlap_genes_data, threshold, expected_result, overlap_group_factory):
    main_gene = Gene(*main_gene_data)
    overlap_genes = [Gene(*gene_data) for gene_data in overlap_genes_data]

    overlap_group = overlap_group_factory(main_gene, *overlap_genes)

    result = overlap_group.overlaps_dont_overlap(threshold)
    assert result == expected_result
