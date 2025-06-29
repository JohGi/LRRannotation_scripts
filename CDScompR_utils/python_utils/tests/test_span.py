import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from CDScompR_lib.gene import Gene

class DummyFeature:
    def __init__(self, start, end, featuretype="gene"):
        self.start = start
        self.end = end
        self.featuretype = featuretype
        self.id = "dummy"

class DummyDB:
    def __init__(self, cds_list):
        self._cds = cds_list

    def children(self, feature, featuretype=None, level=None):
        if featuretype == "CDS":
            return iter(self._cds)
        return iter([])

def test_get_span_gene():
    gene = DummyFeature(100, 200)
    transcript = DummyFeature(150, 180)
    db = DummyDB([])
    assert Gene._get_span(db, gene, transcript, "gene") == (100, 200)

def test_get_span_mrna():
    gene = DummyFeature(100, 200)
    transcript = DummyFeature(150, 180)
    db = DummyDB([])
    assert Gene._get_span(db, gene, transcript, "mRNA") == (150, 180)

def test_get_span_cds():
    gene = DummyFeature(100, 200)
    transcript = DummyFeature(150, 180)
    cds_list = [DummyFeature(160, 165), DummyFeature(170, 175)]
    db = DummyDB(cds_list)
    assert Gene._get_span(db, gene, transcript, "CDS") == (160, 175)

