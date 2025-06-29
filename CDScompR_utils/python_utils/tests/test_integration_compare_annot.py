import subprocess
import difflib
import os
from pathlib import Path

def test_integration_ignore_order(tmp_path: Path) -> None:
    """
    Integration test: run the full pipeline and check that the output TSV matches
    the expected TSV, ignoring line order. Shows a clear diff if they differ.
    """
    
    test_dir: str = os.path.dirname(__file__)
    root_dir: str = os.path.abspath(os.path.join(test_dir, ".."))

    ref_gff = f"{test_dir}/data/test_ref.gff"
    pred_gff = f"{test_dir}/data/test_pred.gff"
    csv_file = f"{test_dir}/data/test_scores.csv"
    expected_tsv = f"{test_dir}/data/expected_output.tsv"

    output_tsv = tmp_path / "observed_output.tsv"

    subprocess.run([
        "python", f"{root_dir}/scripts/compare_annots.py",
        "--ref_gff", ref_gff,
        "--pred_gff", pred_gff,
        "--cdscompr_csv", csv_file,
        "--span_type", "CDS",
        "-o", str(output_tsv)
    ], check=True)


    with open(expected_tsv) as f1, open(output_tsv) as f2:
        expected_lines = f1.readlines()
        observed_lines = f2.readlines()

    expected_header, expected_data = expected_lines[0], expected_lines[1:]
    observed_header, observed_data = observed_lines[0], observed_lines[1:]

    assert expected_header == observed_header, "TSV headers do not match."

    expected_data_sorted = sorted(expected_data)
    observed_data_sorted = sorted(observed_data)

    if expected_data_sorted != observed_data_sorted:
        diff = list(difflib.unified_diff(
            expected_data_sorted,
            observed_data_sorted,
            fromfile="expected_output.tsv",
            tofile="observed_output.tsv",
            lineterm=""
        ))
        diff_text = "\n".join(diff)
        assert False, f"Output TSV data does not match expected (ignoring line order).\nDiff:\n{diff_text}"
