#!/usr/bin/env python3
"""Build a per-gene trace table from completed LRRtransfer outputs."""

import argparse
import csv
import logging
from pathlib import Path
from typing import Iterator, Optional

from attrs import define


LOGGER = logging.getLogger(__name__)

OUTPUT_COLUMNS = [
    "predicted_gene_id",
    "target_locus",
    "target_seqid",
    "target_start",
    "target_end",
    "target_strand",
    "round1_query",
    "round1_method",
    "round1_score",
    "round2_run",
    "round2_query",
    "final_query",
    "final_method",
    "final_score",
    "CL_split_id",
]


@define(frozen=True, slots=True)
class Gene:
    """Store the useful fields of one GFF gene feature."""

    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    origin: str = ""
    method: str = ""
    score: str = ""

    @classmethod
    def from_gff_row(cls, row: list[str]) -> "Gene":
        """Create a gene from a nine-column GFF row."""
        attributes = parse_attributes(row[8])
        comment = parse_comment(attributes.get("comment", ""))

        return cls(
            gene_id=attributes.get("ID", ""),
            seqid=row[0],
            start=int(row[3]),
            end=int(row[4]),
            strand=row[6],
            origin=comment.get("Origin", ""),
            method=comment.get("pred", ""),
            score=comment.get("score", ""),
        )

    def coordinate_key(self) -> tuple[str, int, int, str]:
        """Return the coordinates used to match the merged GFF."""
        return self.seqid, self.start, self.end, self.strand


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build a trace table from an LRRtransfer results directory."
    )
    parser.add_argument(
        "results_dir",
        type=Path,
        help="LRRtransfer results directory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV path. Default: RESULTS_DIR/prediction_trace.tsv",
    )
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    """Parse semicolon-separated GFF attributes."""
    attributes = {}

    for field in text.split(";"):
        key, separator, value = field.partition("=")
        if separator:
            attributes[key] = value

    return attributes


def parse_comment(text: str) -> dict[str, str]:
    """Parse slash-separated LRRtransfer comment fields."""
    fields = {}

    for field in text.split(" / "):
        key, separator, value = field.partition(":")
        if separator:
            fields[key] = value.strip()

    return fields


def iter_genes(path: Path) -> Iterator[Gene]:
    """Yield gene features from a GFF file."""
    with path.open() as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")

            if len(row) == 9 and row[2] == "gene":
                yield Gene.from_gff_row(row)


def read_first_gene(path: Path) -> Optional[Gene]:
    """Return the first gene from a GFF file, if present."""
    if not path.is_file() or path.stat().st_size == 0:
        return None

    return next(iter_genes(path), None)


def read_query_targets(path: Path) -> dict[str, str]:
    """Map each candidate locus to its round-one query gene."""
    associations = {}

    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue

            target, query, _strand = line.rstrip("\n").split("\t")
            associations[target] = query

    return associations


def read_round2_query(path: Path) -> str:
    """Read the query selected for the optional second round."""
    if not path.is_file() or path.stat().st_size == 0:
        return ""

    fields = path.read_text().splitlines()[0].split("\t")
    return fields[1]


def validate_results_dir(results_dir: Path) -> None:
    """Check that the required completed-run outputs exist."""
    required_paths = [
        results_dir / "annot_best.gff",
        results_dir / "list_query_target.txt",
        results_dir / "filtered_candidatsLRR.gff",
        results_dir / "annotate_one",
    ]

    missing_paths = [path for path in required_paths if not path.exists()]

    if missing_paths:
        formatted_paths = "\n".join(f"  - {path}" for path in missing_paths)
        raise FileNotFoundError(
            f"Missing required LRRtransfer outputs:\n{formatted_paths}"
        )


def build_trace_rows(results_dir: Path) -> list[list[str]]:
    """Build one trace row for each non-empty final prediction."""
    annotate_one_dir = results_dir / "annotate_one"

    round1_queries = read_query_targets(
        results_dir / "list_query_target.txt"
    )

    candidate_genes = {
        gene.gene_id: gene
        for gene in iter_genes(
            results_dir / "filtered_candidatsLRR.gff"
        )
    }

    final_gene_ids = {
        gene.coordinate_key(): gene.gene_id
        for gene in iter_genes(results_dir / "annot_best.gff")
    }

    rows = []

    for best_path in sorted(
        annotate_one_dir.glob("annotate_one_*_best.gff")
    ):
        final_gene = read_first_gene(best_path)

        if final_gene is None:
            continue

        file_prefix = best_path.name.removesuffix("_best.gff")
        split_id = file_prefix.removeprefix("annotate_one_")

        round1_gene = read_first_gene(
            best_path.with_name(f"{file_prefix}_best1.gff")
        )

        round2_query = read_round2_query(
            best_path.with_name(f"{file_prefix}_pairID")
        )

        target_locus = final_gene.gene_id
        target_gene = candidate_genes.get(target_locus)
        round1_query = round1_queries.get(target_locus, "NA")

        rows.append(
            [
                final_gene_ids.get(final_gene.coordinate_key(), "NA"),
                target_locus,
                target_gene.seqid if target_gene else "NA",
                str(target_gene.start) if target_gene else "NA",
                str(target_gene.end) if target_gene else "NA",
                target_gene.strand if target_gene else "NA",
                round1_query,
                round1_gene.method if round1_gene else "NA",
                round1_gene.score if round1_gene else "NA",
                "true" if round2_query else "false",
                round2_query or "NA",
                round2_query or round1_query,
                final_gene.method or "NA",
                final_gene.score or "NA",
                split_id,
            ]
        )

    return rows


def write_tsv(rows: list[list[str]], output_path: Path) -> None:
    """Write the trace table as a tab-separated file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)


def main() -> None:
    """Run the trace reconstruction."""
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    results_dir = args.results_dir.resolve()
    output_path = (
        args.output
        or results_dir / "prediction_trace.tsv"
    )

    validate_results_dir(results_dir)
    rows = build_trace_rows(results_dir)
    write_tsv(rows, output_path)

    LOGGER.info(
        "Wrote %d predictions to %s",
        len(rows),
        output_path,
    )


if __name__ == "__main__":
    main()
