import gffutils
from attrs import define

@define
class Protein:
    id: str
    db: gffutils.FeatureDB
    feature: gffutils.Feature

    def cds_length(self) -> int:
        """
        Compute the total length of all CDS features in the transcript.
        """
        return sum(len(cds) for cds in self.db.children(self.feature, featuretype='CDS', level=1))

    def cds_count(self) -> int:
        """
        Count the number of CDS features in the transcript.
        """
        return sum(1 for _ in self.db.children(self.feature, featuretype='CDS', level=1))
