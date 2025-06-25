import gffutils
import pandas as pd
from typing import Optional
from attrs import define
from .protein import Protein

@define
class Gene:
    id: str
    #start: int
    #end: int
    span_start: int
    span_end: int
    protein: Protein
    is_ref: bool
    uid: str
    best_hit_id: Optional[str] = None
    identity_score: Optional[float] = None

    def set_identity_scores(self, score_df: pd.DataFrame, is_ref: bool) -> None:
        """
        Update this Gene with best hit ID and identity score from a CDScompR score dataframe.
        """
        self_id_col = "ref_id" if is_ref else "alt_id"
        best_hit_id_col = "alt_id" if is_ref else "ref_id"

        score_row = score_df[score_df[self_id_col] == self.id]
        if not score_row.empty:
            self.best_hit_id = score_row.iloc[0][best_hit_id_col]
            score_value = score_row.iloc[0]["identity_score"]
            self.identity_score = float(score_value) if pd.notna(score_value) else None


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
