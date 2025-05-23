#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd

# Add src/ to the PYTHONPATH to import the module
sys.path.insert(
    0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
)

from CDScompR_lib.plots import plot_identity_hist_from_csv


def main():
    parser = argparse.ArgumentParser(
        description="Plot histogram of identity scores from a CDScompR CSV file"
    )
    parser.add_argument(
        "--csv", required=True, help="Path to the CDScompR CSV file to process"
    )
    parser.add_argument(
        "--ref-name",
        required=True,
        help="Ref name for labeling output and plot",
    )
    parser.add_argument(
        "--alt-name",
        required=True,
        help="Alt name for labeling output and plot",
    )
    parser.add_argument(
        "--output",
        help="Output file path (PNG)",
    )

    args = parser.parse_args()

    plot_identity_hist_from_csv(args.csv, args.ref_name, args.alt_name, args.output)


if __name__ == "__main__":
    main()
