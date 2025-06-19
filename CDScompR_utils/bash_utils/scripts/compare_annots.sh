#!/bin/bash
set -euo pipefail

### LAUNCHING:
# ./compare_annots.sh <CA_config.sh> <ref.gff> <alt.gff> <output_suffix> <path/to/output_dir> [--debug]

## With CA_config.sh setting-up the following variables (var=value, one variable per ligne):
# GMT_DIR: path to GeneModelTransfer cloned repo
# GMT_SIF: path to GeneModelTransfer .sif
# PYTHON_LIBS_SIF: path to python libraries .sif (compare_annots.py dependencies)
# CDSCOMPR_DIR: path to CDScompR cloned repo
# CDSCOMPR_UTILS_DIR: path to CDScompR utils directory (https://github.com/JohGi/LRRannotation_scripts/tree/main/CDScompR_utils/python_utils)

module load singularity/3.6.3
module load python/packages/3.8.2

## Functions
check_variables_exist() {
  for varname in "$@"; do
    if [[ -z "${!varname+x}" || -z "${!varname}" ]]; then
      echo "Error: variable '$varname' is not set or is empty." >&2
      exit 1
    fi
  done
}

check_files_exist() {
  for file in "$@"; do
    if [[ ! -f "$file" ]]; then
      echo "Error: file '$file' does not exist." >&2
      exit 1
    fi
  done
}

check_folders_exist() {
  for folder in "$@"; do
    if [[ ! -d "$folder" ]]; then
      echo "Error: folder '$folder' does not exist." >&2
      exit 1
    fi
  done
}

import_and_check_variables() {
  CONFIG_FILE=$(realpath $1)
  REF_GFF=$(realpath $2)
  ALT_GFF=$(realpath $3)
  OUTPUT_SUFFIX=$4
  OUT_DIR=$5

  if [[ "${6:-}" == "--debug" ]]; then
    echo "Debug mode"
    set -xv
  fi

  check_variables_exist CONFIG_FILE
  check_files_exist "$CONFIG_FILE"
  source "${CONFIG_FILE}"

  check_variables_exist GMT_DIR GMT_SIF PYTHON_LIBS_SIF CDSCOMPR_DIR REF_GFF ALT_GFF OUTPUT_SUFFIX OUT_DIR #CDSCOMPR_UTILS_DIR
  check_files_exist "$GMT_SIF" "$PYTHON_LIBS_SIF" "$REF_GFF" "$ALT_GFF"
  check_folders_exist "$GMT_DIR" "$CDSCOMPR_DIR" #"$CDSCOMPR_UTILS_DIR"

  SCRIPT_DIR="$(dirname "$(realpath "$0")")"
  CDSCOMPR_UTILS_DIR=$(realpath ${SCRIPT_DIR}/../..)
}

sort_gffs() {
  local ref_gff=$1
  local alt_gff=$2
  local out_dir=$3

  mkdir -p "${out_dir}/01_sorted_input_gffs"

  local ref_sorted="${out_dir}/01_sorted_input_gffs/ref_"$(basename "$ref_gff" .gff)"_sorted.gff"
  singularity exec "$GMT_SIF" python3 "${GMT_DIR}/SCRIPT/sort_gff.py" -g "$ref_gff" -o "${ref_sorted}"
  local alt_sorted="${out_dir}/01_sorted_input_gffs/alt_"$(basename "$alt_gff" .gff)"_sorted.gff"
  singularity exec "$GMT_SIF" python3 "${GMT_DIR}/SCRIPT/sort_gff.py" -g "$alt_gff" -o "${alt_sorted}"
  echo $(realpath --relative-to="$PWD" "$ref_sorted" "$alt_sorted")
}

run_CDScompR() {
  local sorted_ref_gff=$(realpath $1)
  local sorted_alt_gff=$(realpath $2)
  local suffix=$3
  local out_dir=$(realpath $4)

  cd "${out_dir}" && python3 ${CDSCOMPR_DIR}/CDScompR.py --verbose --reference "$sorted_ref_gff" --alternative "$sorted_alt_gff" >"${out_dir}/CDScompR.log" 2>&1 && cd - >/dev/null
  mv ${out_dir}/results ${out_dir}/02_CDScompR_results
  awk -F ',' '{if (NR > 1 && $3 != "~" && $4 != "~"){print $7}}' ${out_dir}/02_CDScompR_results/*.csv | sort -n | uniq -c >${out_dir}/02_CDScompR_results/${suffix}_overlaping_genes_score_distr.txt
  echo $(realpath --relative-to="$PWD" "${out_dir}/02_CDScompR_results/*.csv")
}

main() {
  import_and_check_variables "$@"
  source ${GMT_DIR}/bin/lib_gff_comment.sh

  mkdir -p "${OUT_DIR}"
  read sorted_ref_gff sorted_alt_gff < <(sort_gffs "$REF_GFF" "$ALT_GFF" "$OUT_DIR")

  CDScompR_output_csv=$(run_CDScompR "${sorted_ref_gff}" "${sorted_alt_gff}" "${OUTPUT_SUFFIX}" "${OUT_DIR}")

  mkdir -p "${OUT_DIR}/03_overlaps"
  singularity run "$PYTHON_LIBS_SIF" "${CDSCOMPR_UTILS_DIR}/python_utils/scripts/compare_annots.py" --ref_gff "${sorted_ref_gff}" --pred_gff "${sorted_alt_gff}" --cdscompr_csv "${CDScompR_output_csv}" -o "${OUT_DIR}/03_overlaps/${OUTPUT_SUFFIX}_overlaps.tsv" 2>&1 | tee "${OUT_DIR}/overlaps.log"

  cut -f11 "${OUT_DIR}/03_overlaps/${OUTPUT_SUFFIX}_overlaps.tsv" | tail -n +2 | sort | uniq -c | awk '{print $2":\t"$1}' >"${OUT_DIR}/03_overlaps/${OUTPUT_SUFFIX}_summary.tsv"
}

main "$@"
