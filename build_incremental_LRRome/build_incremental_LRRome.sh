#!/bin/bash

source ./lib_LRRome.sh

## ----------------------------------- MAIN ------------------------------------------------ ##

main(){
  gff_list=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2025_01/01_raw_data/01_gff_EXP/gff_list_order.txt
  initial_LRRome=/storage/replicated/cirad/projects/GE2POP/2023_LRR/IRGSP/LRRome
  initial_LRR_gff=/storage/replicated/cirad/projects/GE2POP/2023_LRR/IRGSP/Oryza_Nipponbare_IRGSP-1.0_LRR-CR__20220209.gff
  exp_ref_genome=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3/Data_Package_01_12_22/DW_Svevo4.fasta
  exp_prefix=DWSvevo3January
  init_prefix=IRGSP
  seq_type=FSprot # FSprot (will extract sequences with Extract_sequences_from_genome.py accounting for frameshifts) or prot (will extract sequences with AGAT accounting for CDS phase)

  mkdir -p 02_build_exp_LRRome 03_LRRome 04_final_GFF

  clean_gff $(realpath 01_gff_EXP) $exp_prefix

  build_exp_LRRome_multiGFF $exp_prefix $exp_ref_genome $gff_list $seq_type

  #build_exp_LRRome $exp_prefix $exp_ref_folder

  exp_LRRome=$(realpath 02_build_exp_LRRome/LRRome)
  merge_LRRome $initial_LRRome $exp_LRRome 03_LRRome

  concat_gff $init_prefix $exp_prefix $initial_LRR_gff

  exp_LRR_gff=$(realpath 02_build_exp_LRRome/LRR_ANNOT/${exp_prefix}_LRR.gff)
  create_info_locus 04_final_GFF/${init_prefix}_${exp_prefix}_LRR.gff $exp_prefix $init_prefix

  compute_gene_stats 04_final_GFF/${init_prefix}_${exp_prefix}_LRR.gff 04_final_GFF/${init_prefix}_${exp_prefix}_gene_stats.tsv

  write_infos $exp_prefix

}

main

# sbatch --partition=agap_normal --mem=20G --wrap="/home/girodollej/scratch/2024_LRR/03_scripts/06_LRRTRANSFER_LAUNCH_PREP/build_incremental_LRRome.sh"
