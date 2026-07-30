#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import sys
from itertools import combinations
from pathlib import Path


MISSING_GENE_VALUES = {"", "~"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize a CDScompare multi-annotation synthesis file in best mode. "
            "The first column contains reference genes. Each following pair of "
            "columns contains a matched gene and its score for one annotation."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input synthesis file in CSV, TSV, or semicolon-separated format.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV summary.",
    )
    parser.add_argument(
        "--prefixes",
        nargs="+",
        help=(
            "Optional annotation labels, in column-pair order. Their number must "
            "equal (number of columns - 1) / 2. When omitted, annot1, annot2, ... "
            "are used."
        ),
    )
    return parser.parse_args()


def sniff_dialect(path: Path) -> csv.Dialect:
    with path.open("r", newline="") as handle:
        sample = handle.read(65536)

    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        return csv.get_dialect("excel")


def is_missing_gene(value: str | None) -> bool:
    return (value or "").strip() in MISSING_GENE_VALUES


def parse_score(
    value: str | None,
    *,
    row_number: int,
    column_name: str,
    gene_found: bool,
) -> float | None:
    if not gene_found:
        return None

    text = (value or "").strip()
    if not text:
        raise ValueError(
            f"Missing score at row {row_number}, column '{column_name}', "
            "although the matched gene is present."
        )

    try:
        score = float(text)
    except ValueError as exc:
        raise ValueError(
            f"Invalid score '{text}' at row {row_number}, "
            f"column '{column_name}'."
        ) from exc

    if not math.isfinite(score):
        raise ValueError(
            f"Non-finite score '{text}' at row {row_number}, "
            f"column '{column_name}'."
        )

    return score


def is_score_100(score: float) -> bool:
    return math.isclose(score, 100.0, rel_tol=0.0, abs_tol=1e-9)


def scores_equal(score1: float, score2: float) -> bool:
    return math.isclose(score1, score2, rel_tol=0.0, abs_tol=1e-9)


def add_row(
    output_rows: list[dict[str, str | int]],
    annot: str,
    overlap: str,
    score: str,
    n: int,
) -> None:
    output_rows.append(
        {
            "Annot": annot,
            "Overlap": overlap,
            "Score": score,
            "N": n,
        }
    )


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    dialect = sniff_dialect(input_path)

    with input_path.open("r", newline="") as handle:
        reader = csv.reader(handle, dialect=dialect)

        try:
            columns = next(reader)
        except StopIteration as exc:
            raise ValueError("Input file is empty.") from exc

        n_columns = len(columns)

        if n_columns < 5:
            raise ValueError(
                f"Expected at least 5 columns, found {n_columns}."
            )

        if n_columns % 2 == 0:
            raise ValueError(
                f"Expected an odd number of columns, found {n_columns}."
            )

        n_annotations = (n_columns - 1) // 2

        if n_columns > 5:
            print(
                "[WARNING] This input contains more than 5 columns "
                f"(ref + {n_annotations} annotations). The script has not yet "
                "been tested with more than two compared annotations "
                "(ref + 2 annotations).",
                file=sys.stderr,
            )

        if args.prefixes is None:
            prefixes = [
                f"annot{i}"
                for i in range(1, n_annotations + 1)
            ]
        else:
            prefixes = args.prefixes

            if len(prefixes) != n_annotations:
                raise ValueError(
                    f"Expected {n_annotations} prefixes for {n_columns} columns, "
                    f"but received {len(prefixes)}: {prefixes}"
                )

            if len(set(prefixes)) != len(prefixes):
                raise ValueError("Prefixes must be unique.")

        annotations = []
        for index, prefix in enumerate(prefixes):
            gene_col_index = 1 + 2 * index
            score_col_index = gene_col_index + 1

            annotations.append(
                {
                    "prefix": prefix,
                    "gene_col_index": gene_col_index,
                    "score_col_index": score_col_index,
                    "gene_col_name": columns[gene_col_index],
                    "score_col_name": columns[score_col_index],
                }
            )

        per_annot_counts = {
            prefix: {
                "score_100": 0,
                "score_lt_100": 0,
                "missed": 0,
            }
            for prefix in prefixes
        }

        pair_counts = {}
        for annot1, annot2 in combinations(annotations, 2):
            prefix1 = annot1["prefix"]
            prefix2 = annot2["prefix"]

            pair_counts[(prefix1, prefix2)] = {
                "annot1_found_annot2_missed_score_100": 0,
                "annot1_found_annot2_missed_score_lt_100": 0,
                "annot2_found_annot1_missed_score_100": 0,
                "annot2_found_annot1_missed_score_lt_100": 0,
                "both_scores_100": 0,
                "same_score_lt_100": 0,
                "annot1_better": 0,
                "annot2_better": 0,
            }

        for row_number, values in enumerate(reader, start=2):
            if not values or all(not value.strip() for value in values):
                continue

            if len(values) != n_columns:
                raise ValueError(
                    f"Row {row_number} has {len(values)} columns; "
                    f"expected {n_columns}."
                )

            row_data = {}

            for annotation in annotations:
                prefix = annotation["prefix"]
                gene_value = values[annotation["gene_col_index"]]
                score_value = values[annotation["score_col_index"]]

                found = not is_missing_gene(gene_value)
                score = parse_score(
                    score_value,
                    row_number=row_number,
                    column_name=annotation["score_col_name"],
                    gene_found=found,
                )

                row_data[prefix] = {
                    "found": found,
                    "score": score,
                }

                if not found:
                    per_annot_counts[prefix]["missed"] += 1
                elif is_score_100(score):
                    per_annot_counts[prefix]["score_100"] += 1
                elif score < 100:
                    per_annot_counts[prefix]["score_lt_100"] += 1
                else:
                    raise ValueError(
                        f"Score above 100 at row {row_number}, "
                        f"column '{annotation['score_col_name']}': {score}"
                    )

            for annot1, annot2 in combinations(annotations, 2):
                prefix1 = annot1["prefix"]
                prefix2 = annot2["prefix"]
                counts = pair_counts[(prefix1, prefix2)]

                found1 = row_data[prefix1]["found"]
                found2 = row_data[prefix2]["found"]
                score1 = row_data[prefix1]["score"]
                score2 = row_data[prefix2]["score"]

                if found1 and not found2:
                    if is_score_100(score1):
                        counts[
                            "annot1_found_annot2_missed_score_100"
                        ] += 1
                    elif score1 < 100:
                        counts[
                            "annot1_found_annot2_missed_score_lt_100"
                        ] += 1

                if found2 and not found1:
                    if is_score_100(score2):
                        counts[
                            "annot2_found_annot1_missed_score_100"
                        ] += 1
                    elif score2 < 100:
                        counts[
                            "annot2_found_annot1_missed_score_lt_100"
                        ] += 1

                if found1 and found2:
                    if is_score_100(score1) and is_score_100(score2):
                        counts["both_scores_100"] += 1
                    elif scores_equal(score1, score2):
                        if score1 < 100:
                            counts["same_score_lt_100"] += 1
                    elif score1 > score2:
                        counts["annot1_better"] += 1
                    else:
                        counts["annot2_better"] += 1

    output_rows: list[dict[str, str | int]] = []

    for prefix in prefixes:
        counts = per_annot_counts[prefix]

        add_row(
            output_rows,
            prefix,
            "with_ref",
            "==100",
            counts["score_100"],
        )
        add_row(
            output_rows,
            prefix,
            "with_ref",
            "<100",
            counts["score_lt_100"],
        )
        add_row(
            output_rows,
            "ref",
            f"missed_by_{prefix}",
            "~",
            counts["missed"],
        )

    for annot1, annot2 in combinations(annotations, 2):
        prefix1 = annot1["prefix"]
        prefix2 = annot2["prefix"]
        counts = pair_counts[(prefix1, prefix2)]

        add_row(
            output_rows,
            prefix1,
            f"with_ref,_{prefix2}_missed",
            "==100",
            counts["annot1_found_annot2_missed_score_100"],
        )
        add_row(
            output_rows,
            prefix1,
            f"with_ref,_{prefix2}_missed",
            "<100",
            counts["annot1_found_annot2_missed_score_lt_100"],
        )
        add_row(
            output_rows,
            prefix2,
            f"with_ref,_{prefix1}_missed",
            "==100",
            counts["annot2_found_annot1_missed_score_100"],
        )
        add_row(
            output_rows,
            prefix2,
            f"with_ref,_{prefix1}_missed",
            "<100",
            counts["annot2_found_annot1_missed_score_lt_100"],
        )
        add_row(
            output_rows,
            prefix1,
            f"with_ref_and_{prefix2}",
            "all_scores==100",
            counts["both_scores_100"],
        )
        add_row(
            output_rows,
            prefix1,
            f"with_ref_and_{prefix2}",
            f"=={prefix2}_score,_both<100",
            counts["same_score_lt_100"],
        )
        add_row(
            output_rows,
            prefix1,
            f"with_ref_and_{prefix2}",
            f">{prefix2}_score",
            counts["annot1_better"],
        )
        add_row(
            output_rows,
            prefix2,
            f"with_ref_and_{prefix1}",
            f">{prefix1}_score",
            counts["annot2_better"],
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Annot", "Overlap", "Score", "N"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Detected annotations: {', '.join(prefixes)}")
    print(f"Summary written to: {output_path}")


if __name__ == "__main__":
    main()
