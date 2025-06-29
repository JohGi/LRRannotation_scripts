import gffutils
import tempfile

def build_db(gff_path: str) -> gffutils.FeatureDB:
    """
    Create a gffutils database from a GFF file.
    """
    db_path = tempfile.NamedTemporaryFile(delete=True).name
    return gffutils.create_db(
        gff_path,
        dbfn=db_path,
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )

