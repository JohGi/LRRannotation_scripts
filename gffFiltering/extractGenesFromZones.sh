#!/bin/bash

set -euo pipefail

if ! command -v AGAT &> /dev/null; then
  module load AGAT/1.2.0-singularity
fi





usage() {
    echo "Usage: $0 <annot.gff> <regions.tsv> <output_folder>"
    echo
    echo "This script creates a gff containing only the genes included or overlapping with a list of regions (chr start end)."
    echo
    echo "Arguments:"
    echo "  <annot.gff>    The input annotation file."
    echo "  <zones.tsv>    The input region file (one region per line)."
    echo -e "  <output_folder>    The output folder name.\n"
}

agat_sp_filter_record_by_coordinates(){
  local gff=$1
  local zones=$2
  mkdir -p AGAT_output
  cd AGAT_output
  agat_sp_filter_record_by_coordinates.pl --gff $gff --tsv $zones
  cd ..
}

concatGFFperChr() {
  local zones=$1
  mkdir -p GFFperChr
  for chr in $(cut -f1 $zones | sort | uniq) ; do
    cat AGAT_output/filter_record_by_coordinates/${chr}*.gff3 > GFFperChr/${chr}_in_zones.gff
  done
}


### Main ###
if [ "$#" -ne 3 ]; then
    usage
    exit 1
fi

gff=$1
zones=$2
output_folder=$3
gff=$(realpath $gff)
zones=$(realpath $zones)

mkdir -p $output_folder
cd $output_folder
agat_sp_filter_record_by_coordinates $gff $zones
concatGFFperChr $zones
cd ..
