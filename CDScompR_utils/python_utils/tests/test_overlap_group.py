import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import pytest
from CDScompR_lib.gene import Gene
from CDScompR_lib.overlap_group import OverlapGroup


class DummyFeature:
    def __init__(self, seqid: str):
        self.seqid = seqid


class DummyProtein:
    def __init__(self, seqid: str):
        self.feature = DummyFeature(seqid)

    def cds_length(self) -> int:
        return 100

    def cds_count(self) -> int:
        return 2


def make_gene(id: str, start: int, end: int, is_ref: bool, seqid: str = "chr") -> Gene:
    return Gene(
        id=id,
        start=start,
        end=end,
        protein=DummyProtein(seqid),
        is_ref=is_ref,
        uid=f"{'ref' if is_ref else 'pred'}:{id}"
    )

@pytest.mark.parametrize("ref_genes, pred_genes, expected_nb_groups, expected_types", [
    # Ref and pred overlap → match
    ([make_gene("ref1", 100, 200, True)],
     [make_gene("pred1", 150, 250, False)],
     1, ["match"]),

    # Ref without overlap → disappearance
    ([make_gene("ref2", 300, 400, True)], [], 1, ["disappearance"]),

    # Pred without overlap → appearance
    ([], [make_gene("pred2", 500, 600, False)], 1, ["appearance"]),

    # Several pred overlap with 1 ref → split
    ([make_gene("ref3", 700, 800, True)],
     [make_gene("pred3a", 710, 740, False), make_gene("pred3b", 750, 780, False)],
     1, ["split"]),

    # Several ref overlap with 1 pred → fusion
    ([make_gene("ref4a", 900, 930, True), make_gene("ref4b", 940, 970, True)],
     [make_gene("pred4", 905, 960, False)],
     1, ["fusion"]),

    # Several ref and several pred in cluster → complex
    ([make_gene("ref5a", 1000, 1030, True), make_gene("ref5b", 1040, 1070, True)],
     [make_gene("pred5a", 1005, 1060, False), make_gene("pred5b", 1065, 1090, False)],
     1, ["complex"]),

    # No gene → empty list
    ([], [], 0, []),

    # Genes on different chromosomes → must be in separate groups
    ([make_gene("ref6", 100, 200, True, seqid="chr1")],
     [make_gene("pred6", 150, 250, False, seqid="chr2")],
     2, ["disappearance", "appearance"])
])
def test_overlap_group_types(ref_genes, pred_genes, expected_nb_groups, expected_types):
    groups = OverlapGroup.overlap_groups_from_genes(ref_genes, pred_genes)
    assert len(groups) == expected_nb_groups
    assert sorted(g.get_type() for g in groups) == sorted(expected_types)
