import gffutils
from typing import Optional
from attrs import define
from .protein import Protein

@define
class Gene:
    id: str
    span_start: int
    span_end: int
    protein: Protein
    is_ref: bool
    uid: str
    best_hit_id: Optional[str] = None
    identity_score: Optional[float] = None


    def set_identity_scores(self, score_lookup: dict[str, tuple[str, float]]) -> None:
        """
        Set the best hit ID and identity score for this Gene from a lookup dictionary.
        """
        if self.id in score_lookup:
            best_hit_id, identity_score = score_lookup[self.id]
            self.best_hit_id = best_hit_id
            self.identity_score = identity_score


    @staticmethod
    def _get_span(db: gffutils.FeatureDB, gene: gffutils.Feature, transcript: gffutils.Feature, span_type: str) -> tuple[int, int]:
        if span_type == "gene":
            return gene.start, gene.end
        elif span_type == "mRNA":
            return transcript.start, transcript.end
        elif span_type == "CDS":
            cds_coords = [(cds.start, cds.end) for cds in db.children(transcript, featuretype="CDS", level=1)]
            if not cds_coords:
                raise ValueError(f"No CDS found for gene {gene.id}")
            return min(s for s, _ in cds_coords), max(e for _, e in cds_coords)
        else:
            raise ValueError(f"Unsupported span_type: {span_type} (accepted span types are 'gene', 'mRNA' and 'CDS')")


    @classmethod
    def from_gff(cls, db: gffutils.FeatureDB, gene: gffutils.Feature, is_ref: bool, span_type: str = "gene") -> "Gene":
        """
        Create a Gene object from a GFF feature and its associated transcript.
        """
        transcripts = list(db.children(gene, featuretype=('mRNA', 'transcript'), level=1))
        if len(transcripts) != 1:
            raise ValueError(f"Gene {gene.id} has {len(transcripts)} transcripts (expected exactly 1) — One mRNA is required and alternative splicing is currently not supported.")
        transcript = transcripts[0]
        protein = Protein(id=transcript.id, db=db, feature=transcript)
        uid = f"{'ref' if is_ref else 'pred'}:{gene.id}"

        span_start, span_end = Gene._get_span(db, gene, transcript, span_type)

        return cls(id=gene.id, span_start=span_start, span_end=span_end, protein=protein, is_ref=is_ref, uid=uid)
