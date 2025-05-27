#!/bin/bash
set -euo pipefail

### LAUNCHING:
# ./build_incremental_LRRome.sh BILRRome_config.sh

## With BILRRome_config.sh setting-up the following variables (var=value, one variable per ligne):
# GMT_DIR: path to GeneModelTransfer cloned repo
# GMT_SIF: path to GeneModelTransfer .sif
# LRRPROFILER_SIF: path to LRRprofiler .sif
# GFF_LIST: path to a txt file listing all expertised gff
# INITIAL_LRROME: path the initial LRRome (built by GMT create_LRRome.sh)
# INITIAL_LRR_GFF: path to the GFF associated with the initial LRRome
# EXP_REF_GENOME: path to the reference genome of the expertised genes
# EXP_PREFIX: prefix for the expertised LRRome (eg: DWSvevo3January)
# INIT_PREFIX: prefix for the initial LRRome (eg: IRGSP)
# SEQ_TYPE: either 'FSprot' (will extract sequences with Extract_sequences_from_genome.py accounting for frameshifts) or 'prot' (will extract sequences with AGAT accounting for CDS phase)


## ----------------------------------- MAIN ------------------------------------------------ ##

main(){

  # importing variables and functions
  CONFIG_FILE=$1

  SCRIPT_DIR="$(realpath "$(dirname "${BASH_SOURCE[0]}")")"
  source ${SCRIPT_DIR}/lib_LRRome.sh

  chmod 755 $CONFIG_FILE && source $CONFIG_FILE
  require_variables GMT_DIR GMT_SIF LRRPROFILER_SIF GFF_LIST INITIAL_LRROME INITIAL_LRR_GFF EXP_REF_GENOME EXP_PREFIX INIT_PREFIX SEQ_TYPE 

  CONCAT_AND_RM_REPEAT_GENES=${SCRIPT_DIR}/concatAndRmRepeatGenes.py

  # building new LRRome
  #mkdir -p 02_build_exp_LRRome 03_LRRome 04_final_GFF

  exp_LRRome_out_dir=01_build_exp_LRRome
  # clean_gff $GFF_LIST $EXP_PREFIX ${GMT_DIR}/SCRIPT $exp_LRRome_out_dir
  
  # build_exp_LRRome_multiGFF $EXP_PREFIX $EXP_REF_GENOME ${exp_LRRome_out_dir}/clean_gff.list $GMT_SIF $SEQ_TYPE ${CONCAT_AND_RM_REPEAT_GENES} ${LRRPROFILER_SIF} ${GMT_DIR}/SCRIPT $exp_LRRome_out_dir
  exp_LRRome=$(realpath ${exp_LRRome_out_dir}/LRRome)
  exp_LRRome_gff=$(realpath ${exp_LRRome_out_dir}/LRR_ANNOT/${EXP_PREFIX}_LRR.gff)

  final_LRRome_out_dir=02_LRRome
  # merge_LRRome $INITIAL_LRROME $exp_LRRome $final_LRRome_out_dir

  final_gff_out_dir=03_final_GFF
  # concat_gff $INITIAL_LRR_GFF $exp_LRRome_gff ${final_gff_out_dir} ${INIT_PREFIX}_${EXP_PREFIX}_LRR.gff
  final_gff=${final_gff_out_dir}/${INIT_PREFIX}_${EXP_PREFIX}_LRR.gff

  create_info_locus ${final_gff} ${final_gff_out_dir}/${INIT_PREFIX}_${EXP_PREFIX}_info_locus.txt

  compute_gene_stats ${final_gff} ${final_gff_out_dir}/${INIT_PREFIX}_${EXP_PREFIX}_gene_stats.tsv

  write_infos $INITIAL_LRROME $GFF_LIST $INIT_PREFIX $EXP_PREFIX $exp_LRRome_out_dir $final_LRRome_out_dir $final_gff_out_dir
}

main "$@"

