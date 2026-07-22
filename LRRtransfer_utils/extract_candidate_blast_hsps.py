#!/usr/bin/env python3
"""Extract TBLASTN HSPs for one predicted LRRtransfer gene."""

import argparse
import csv
import logging
from pathlib import Path
from typing import Optional


LOGGER = logging.getLogger(__name__)
Region = tuple[str, int, int]
BLAST_COLUMNS = (
    "qseqid sseqid qlen length qstart qend sstart send nident pident "
    "gapopen evalue bitscore positive"
).split()
EXTRA_COLUMNS = [
    "is_round1_query",
    "is_round2_query",
    "is_in_selected_hsp_path",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument(
        "--pred-gene-id",
        required=True,
        help="Predicted gene ID from results/annot_best.gff.",
    )
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    """Parse GFF attributes."""
    return {
        key: value
        for field in text.split(";")
        if (parts := field.partition("="))[1]
        for key, value in [(parts[0], parts[2])]
    }


def read_gene_region(path: Path, gene_id: Optional[str] = None) -> tuple[str, Region]:
    """Return a GFF gene ID and region."""
    with path.open() as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            if len(row) != 9 or row[2] != "gene":
                continue

            current_id = parse_attributes(row[8]).get("ID", "")
            if gene_id is None or current_id == gene_id:
                return current_id, (row[0], int(row[3]), int(row[4]))

    raise ValueError(f"Gene not found in {path}: {gene_id or 'first gene'}")


def find_prediction(
    annotate_one_dir: Path,
    final_region: Region,
) -> tuple[Path, str, str]:
    """Find the individual prediction matching the final gene."""
    matches = []

    for path in annotate_one_dir.glob("annotate_one_*_best.gff"):
        if path.stat().st_size == 0:
            continue
        target_locus, region = read_gene_region(path)
        if region == final_region:
            matches.append((path, target_locus))

    if len(matches) != 1:
        raise ValueError(f"Expected one individual prediction, found {len(matches)}.")

    path, target_locus = matches[0]
    split_id = path.name.removeprefix("annotate_one_").removesuffix("_best.gff")
    return path, target_locus, split_id


def read_round1_query(path: Path, target_locus: str) -> str:
    """Return the initial query associated with a candidate locus."""
    with path.open() as handle:
        for line in handle:
            target, query, _strand = line.rstrip().split("\t")
            if target == target_locus:
                return query
    raise ValueError(f"Candidate locus not found in {path}: {target_locus}")


def read_round2_query(best_path: Path) -> Optional[str]:
    """Return the optional second-round query."""
    prefix = best_path.name.removesuffix("_best.gff")
    pair_id = best_path.with_name(f"{prefix}_pairID")
    if not pair_id.is_file() or pair_id.stat().st_size == 0:
        return None
    return pair_id.read_text().splitlines()[0].split("\t")[1]


def read_candidate_regions(
    path: Path,
    target_locus: str,
) -> tuple[Region, Region, set[tuple[int, int]]]:
    """Return initial, expanded and selected HSP regions."""
    expanded = None
    mrna_id = None
    selected = set()

    with path.open() as handle:
        for line in handle:
            row = line.rstrip().split("\t")
            if len(row) != 9:
                continue

            attributes = parse_attributes(row[8])
            if row[2] == "gene" and attributes.get("ID") == target_locus:
                expanded = (row[0], int(row[3]), int(row[4]))
            elif row[2] == "mRNA" and attributes.get("Parent") == target_locus:
                mrna_id = attributes.get("ID")
            elif row[2] == "CDS" and attributes.get("Parent") == mrna_id:
                selected.add((min(int(row[3]), int(row[4])), max(int(row[3]), int(row[4]))))

    if expanded is None or not selected:
        raise ValueError(f"Incomplete candidate locus in {path}: {target_locus}")

    initial = (expanded[0], min(x[0] for x in selected), max(x[1] for x in selected))
    return initial, expanded, selected


def overlaps(row: list[str], region: Region) -> bool:
    """Return whether a BLAST HSP overlaps a region."""
    start, end = sorted((int(row[6]), int(row[7])))
    return row[1] == region[0] and region[1] <= end and start <= region[2]


def is_selected(
    row: list[str],
    round1_query: str,
    selected: set[tuple[int, int]],
) -> bool:
    """Return whether an HSP is in the selected path."""
    start, end = sorted((int(row[6]), int(row[7])))
    return row[0] == round1_query and (start, end) in selected


def read_hsps(
    path: Path,
    initial: Region,
    expanded: Region,
) -> tuple[list[list[str]], list[list[str]]]:
    """Read HSPs overlapping both candidate regions."""
    initial_rows = []
    expanded_rows = []

    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != len(BLAST_COLUMNS) or not overlaps(row, expanded):
                continue
            expanded_rows.append(row)
            if overlaps(row, initial):
                initial_rows.append(row)

    return initial_rows, expanded_rows


def sort_rows(
    rows: list[list[str]],
    round1_query: str,
    round2_query: Optional[str],
    selected: set[tuple[int, int]],
) -> list[list[str]]:
    """Sort HSPs by query role and genomic position."""
    decorated = []

    for row in rows:
        rank = 0 if row[0] == round1_query else 2
        if round2_query is not None and row[0] == round2_query:
            rank = 1
        start, end = sorted((int(row[6]), int(row[7])))
        key = (rank, row[0], not is_selected(row, round1_query, selected), start, end, int(row[4]))
        decorated.append((key, row))

    return [row for _key, row in sorted(decorated)]


def write_rows(
    path: Path,
    rows: list[list[str]],
    metadata: dict[str, str],
    round1_query: str,
    round2_query: Optional[str],
    selected: set[tuple[int, int]],
    initial: Optional[Region] = None,
) -> None:
    """Write metadata, headers and enriched BLAST rows."""
    columns = BLAST_COLUMNS + EXTRA_COLUMNS + ([] if initial is None else ["is_in_initial_region"])

    with path.open("w", newline="") as handle:
        for key, value in metadata.items():
            handle.write(f"# {key}={value}\n")

        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)

        for row in sort_rows(rows, round1_query, round2_query, selected):
            extra = [
                str(row[0] == round1_query).lower(),
                "NA" if round2_query is None else str(row[0] == round2_query).lower(),
                str(is_selected(row, round1_query, selected)).lower(),
            ]
            if initial is not None:
                extra.append(str(overlaps(row, initial)).lower())
            writer.writerow(row + extra)


def run(results_dir: Path, pred_gene_id: str, output_dir: Path) -> None:
    """Generate initial and expanded HSP reports."""
    _gene_id, final_region = read_gene_region(results_dir / "annot_best.gff", pred_gene_id)
    best_path, target_locus, split_id = find_prediction(results_dir / "annotate_one", final_region)
    round1_query = read_round1_query(results_dir / "list_query_target.txt", target_locus)
    round2_query = read_round2_query(best_path)
    initial, expanded, selected = read_candidate_regions(
        results_dir / "filtered_candidatsLRR.gff",
        target_locus,
    )
    initial_rows, expanded_rows = read_hsps(results_dir / "blast_refProt.tsv", initial, expanded)

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "pred_gene_id": pred_gene_id,
        "target_locus": target_locus,
        "round1_query": round1_query,
        "round2_query": round2_query or "NA",
        "query_target_split_id": split_id,
        "initial_region": f"{initial[0]}:{initial[1]}-{initial[2]}",
        "expanded_region": f"{expanded[0]}:{expanded[1]}-{expanded[2]}",
    }
    initial_path = output_dir / f"{pred_gene_id}_initial_tblastn_hsps.tsv"
    expanded_path = output_dir / f"{pred_gene_id}_expanded_tblastn_hsps.tsv"

    write_rows(
        initial_path,
        initial_rows,
        {**metadata, "selected_region": "initial"},
        round1_query,
        round2_query,
        selected,
    )
    write_rows(
        expanded_path,
        expanded_rows,
        {**metadata, "selected_region": "expanded"},
        round1_query,
        round2_query,
        selected,
        initial,
    )

    LOGGER.info("Initial HSPs: %d -> %s", len(initial_rows), initial_path)
    LOGGER.info("Expanded HSPs: %d -> %s", len(expanded_rows), expanded_path)


def main() -> None:
    """Run the command-line program."""
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    run(
        args.results_dir.resolve(),
        args.pred_gene_id,
        (args.output_dir or args.results_dir).resolve(),
    )


if __name__ == "__main__":
    main()
