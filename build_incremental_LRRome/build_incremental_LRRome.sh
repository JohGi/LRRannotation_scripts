#!/bin/bash

set -euo pipefail

if ! command -v python3 &> /dev/null; then
  module load python/packages/3.8.2
fi
if ! command -v singularity &> /dev/null; then
  module load singularity/3.6.3
fi


scripts=/home/girodollej/scratch/GeneModelTransfer/SCRIPT
LRRprofiler_sif=/storage/replicated/cirad/projects/GE2POP/2023_LRR/USEFUL/LRRprofiler.sif

initial_LRRome=$(realpath /storage/replicated/cirad/projects/GE2POP/2023_LRR/IRGSP/LRRome)
exp_ref_folder=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3/Data_Package_01_12_22
prefix=SVEVO_July


## ------------------------------- FUNCTIONS --------------------------------------------- ##


clean_gff () {             # >> Creates 02_build_exp_LRRome/CLEANED_GFF
  echo -e "\n1/ CLEANING EXP GFF FILES...\n"
  local gff_dir=$1
  local cleaner_chr_prefix=$2
  cd 02_build_exp_LRRome
  mkdir -p CLEANED_GFF GFF_chrOK

  ## Change the chromosome names in the gff files if needed
  echo -e "... Fixing chromosome names...\n"
  for chr_num in $(seq 1 7); do
    for genome in A B; do
      chr=${chr_num}${genome}
      for gff in ${gff_dir}/*${chr}*gff; do
        awk -v chr=$chr 'BEGIN{OFS=FS="\t"}{$1="Chr"chr; print $0}' $gff > GFF_chrOK/$(basename $gff)
      done
    done
  done

  ## Clean the gff files with gff_cleaner
  echo -e "... Running gff_cleaner.py...\n"
  for f in GFF_chrOK/*.gff; do
    prefix=$(basename $f)
    python3 ${scripts}/VR/gff_cleaner.py -p $cleaner_chr_prefix -g $f -o CLEANED_GFF/cleaned_${prefix}
  done > gff_cleaner.out

  ## Count the genes in the gff files
  echo -e "... Counting genes in each gff in 02_build_exp_LRRome/gff_gene_count.txt...\n"
  grep -c -w gene CLEANED_GFF/cleaned_*gff > gff_gene_count.txt

  rm -r GFF_chrOK
  cd ..

}




build_exp_LRRome() {        # >> Creates 02_build_exp_LRRome/LRR_ANNOT, 02_build_exp_LRRome/LRRprofile and 02_build_exp_LRRome/LRRome
  echo -e "\n2/ BUILDING EXPERTISED GENES LRROME...\n"
  local new_prefix=$1
  local exp_ref_folder=$2
  cd 02_build_exp_LRRome

  ## Create a unique exp fasta file
  echo -e "... Building exp reference fasta file in 02_build_exp_LRRome/LRR_ANNOT/"${new_prefix}".fasta...\n"
  mkdir -p LRR_ANNOT
  cat ${exp_ref_folder}/Chr*.fasta > LRR_ANNOT/${new_prefix}.fasta

  ## Detect and classify LRR containing genes with LRRprofiler and remove non LRR genes (write the final gff in LRR_ANNOT)
  echo -e "... Running LRRprofiler and removing non-LRR genes >> see final gff in 02_build_exp_LRRome/LRR_ANNOT/"${new_prefix}"_LRR.gff...\n"
  mkdir -p LRRprofile
  cd LRRprofile
  cat ../CLEANED_GFF/cleaned_*.gff > ${new_prefix}.gff
  python3 ${scripts}/Extract_sequences_from_genome.py -g ${new_prefix}.gff -f ../LRR_ANNOT/${new_prefix}.fasta -o ${new_prefix}_proteins.fasta -t FSprot
  singularity run $LRRprofiler_sif --in_proteome ${new_prefix}_proteins.fasta --name ${new_prefix}_LRRprofiler_output
  source $scripts/../bin/lib_gff_comment.sh
  add_family_info ${new_prefix}.gff $PWD/Res_${new_prefix}_LRRprofiler_output/Res_step3/LRR_classification.csv $PWD/${new_prefix}_LRR_info.gff

  grep -w gene ${new_prefix}_LRR_info.gff | grep -v "Fam=" | cut -f2 -d"=" | cut -f1 -d";" > not_LRR_genes.list
  remove_genes_from_gff ${new_prefix}_LRR_info.gff not_LRR_genes.list ../LRR_ANNOT/${new_prefix}_LRR.gff

  cd ..

  ## Build LRRome (in a LRRome folder)
  echo -e "\n... Running create_LRRome.sh in 02_build_exp_LRRome/LRRome...\n"
  ${scripts}/../bin/create_LRRome.sh $PWD/LRR_ANNOT/${new_prefix}.fasta $PWD/LRR_ANNOT/${new_prefix}_LRR.gff $PWD NULL ${scripts}

  cd ..
}


merge_LRRome() {         # >> Fills 03_LRRome
  echo -e "\n3/ MERGING EXPERTISED AND INITIAL LRROMES...\n"
  local LRRome1=$1
  local LRRome2=$2
  local LRRome_new=$3

  cd $LRRome_new
  mkdir -p REF_EXONS REF_PEP REF_cDNA

  for type in EXONS PEP cDNA; do
    cp ${LRRome1}/REF_${type}/* REF_${type}/
    cp ${LRRome2}/REF_${type}/* REF_${type}/
  done

  for type in exons proteins cDNA; do
    cp ${LRRome1}/REF_${type}.fasta . && chmod u+w REF_${type}.fasta
    cat ${LRRome2}/REF_${type}.fasta >> REF_${type}.fasta
  done

  echo -e "The new LRRome was built in: "$(realpath $LRRome_new)"\n"

  cd - > /dev/null
}


write_infos() {
  local new_prefix=$1
  date > LRRome_incremental_build_infos.txt

  echo -e "\n Input files:" >> LRRome_incremental_build_infos.txt
  echo "- Initial LRRome: "$initial_LRRome >> LRRome_incremental_build_infos.txt
  echo "- Expertised gff files in: 01_gff_EXP" >> LRRome_incremental_build_infos.txt

  echo -e "\n Output files:" >> LRRome_incremental_build_infos.txt
  echo "- Exp reference fasta file: 02_build_exp_LRRome/LRR_ANNOT/"${new_prefix}".fasta" >> LRRome_incremental_build_infos.txt
  echo "- Exp gff in 02_build_exp_LRRome/LRR_ANNOT/"${new_prefix}"_LRR.gff" >> LRRome_incremental_build_infos.txt
  echo "- Exp LRRome: 02_build_exp_LRRome/LRRome" >> LRRome_incremental_build_infos.txt
  echo "- New LRRome: 03_LRRome" >> LRRome_incremental_build_infos.txt

  echo -e "\n\nSee LRRome_incremental_build_infos.txt for run information.\n"
}

## ----------------------------------- MAIN ------------------------------------------------ ##


mkdir -p 02_build_exp_LRRome 03_LRRome

clean_gff $(realpath 01_gff_EXP) $prefix

build_exp_LRRome $prefix $exp_ref_folder

exp_LRRome=$(realpath 02_build_exp_LRRome/LRRome)
merge_LRRome $initial_LRRome $exp_LRRome 03_LRRome

write_infos $prefix





# sbatch --partition=agap_normal --mem=20G --wrap="./build_incremental_LRRome.sh"
