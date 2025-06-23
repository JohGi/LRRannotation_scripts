import gffutils
import pandas as pd
from typing import Optional
from attrs import define
from .protein import Protein

@define
class Gene:
    id: str
    start: int
    end: int
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

    @classmethod
    def from_gff(cls, db: gffutils.FeatureDB, gene: gffutils.Feature, is_ref: bool) -> "Gene":
        """
        Create a Gene object from a GFF feature and its associated transcript.
        """
        transcripts = list(db.children(gene, featuretype=('mRNA', 'transcript'), level=1))
        if len(transcripts) != 1:
            raise ValueError(...)
        transcript = transcripts[0]
        protein = Protein(id=transcript.id, db=db, feature=transcript)
        uid = f"{'ref' if is_ref else 'pred'}:{gene.id}"
        return cls(id=gene.id, start=gene.start, end=gene.end, protein=protein, is_ref=is_ref, uid=uid)
