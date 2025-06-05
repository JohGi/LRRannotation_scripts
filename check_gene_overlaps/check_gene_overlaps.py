"""
Check overlaps between gene features across multiple GFF files.
Only 'gene' features are considered. Prints overlaps with file and line context.
"""

from pathlib import Path
from typing import Dict
from intervaltree import Interval, IntervalTree
import sys

def parse_gff_line(line: str) -> tuple[str, int, int, str] | None:
    """Extract chromosome, start, end, and the full line if it's a 'gene' feature."""
    if line.startswith("#"):
        return None
    parts = line.strip().split("\t")
    if len(parts) < 9 or parts[2].lower() != "gene":
        return None
    chrom = parts[0]
    start = int(parts[3]) - 1  # Convert to 0-based for IntervalTree
    end = int(parts[4])
    return chrom, start, end, line.strip()

def process_gff(file_path: Path, trees_by_chr: Dict[str, IntervalTree]) -> None:
    """Process one GFF file and update interval trees. Print overlaps if found."""
    with file_path.open() as f:
        for line in f:
            parsed = parse_gff_line(line)
            if not parsed:
                continue
            chrom, start, end, original_line = parsed

            tree = trees_by_chr.setdefault(chrom, IntervalTree())
            overlaps = tree.overlap(start, end)

            if overlaps:
                print(f"\n- Overlap found in file: {file_path.name}")
                print(f"\tCurrent gene:\t\t\t{original_line}")
                for ov in sorted(overlaps):
                    print(f"\t-- Overlaps with: {ov.data}")

            tree.add(Interval(start, end, f"[{file_path.name}]\t{original_line}"))

def main(gff_list_file: Path) -> None:
    trees_by_chr: Dict[str, IntervalTree] = {}
    with gff_list_file.open() as f:
        for gff_path in f:
            gff_file = Path(gff_path.strip())
            if gff_file.exists():
                process_gff(gff_file, trees_by_chr)
            else:
                print(f"Error: File not found: {gff_file}", file=sys.stderr)
                sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: check_gene_overlaps.py <gff_list_file>", file=sys.stderr)
        sys.exit(1)

    gff_list_file = Path(sys.argv[1])
    if not gff_list_file.exists():
        print(f"Error: File {gff_list_file} does not exist.", file=sys.stderr)
        sys.exit(1)

    main(gff_list_file)
