#!/bin/bash

set -euo pipefail

if ! command -v python3 &> /dev/null; then
  module load python/packages/3.8.2
fi
if ! command -v singularity &> /dev/null; then
  module load singularity/3.6.3
fi


scripts=/lustre/girodollej/2024_LRR/03_scripts/LRRtransfer/GeneModelTransfer/SCRIPT/
LRRprofiler_sif=/storage/replicated/cirad/projects/GE2POP/2023_LRR/USEFUL/LRRprofiler.sif
concatAndRmRepeatGenes=/home/girodollej/scratch/2024_LRR/03_scripts/LRRannotation_scripts/build_incremental_LRRome/concatAndRmRepeatGenes.py

gff_list=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2025_01/01_raw_data/01_gff_EXP/gff_list_order.txt
initial_LRRome=/storage/replicated/cirad/projects/GE2POP/2023_LRR/IRGSP/LRRome
initial_LRR_gff=/storage/replicated/cirad/projects/GE2POP/2023_LRR/IRGSP/Oryza_Nipponbare_IRGSP-1.0_LRR-CR__20220209.gff
exp_ref_folder=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3/Data_Package_01_12_22
exp_prefix=DWSvevo3January
init_prefix=IRGSP

## ------------------------------- FUNCTIONS --------------------------------------------- ##

rmCR() {
  local file=$1
  if grep -q $'\r' ${file}; then
    sed -i 's/\r$//g' $file
    sed -i 's/\r/\n/g' $file
  fi
}

clean_gff() {             # >> Creates 02_build_exp_LRRome/CLEANED_GFF
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
    python3 ${scripts}/VR/gff_cleaner.py -a -p $cleaner_chr_prefix -g $f -o CLEANED_GFF/cleaned_${prefix}
  done | grep -v "^WARNING: incompatible bounds" > gff_cleaner.out

  ## Count the genes in the gff files
  echo -e "... Counting genes in each gff in 02_build_exp_LRRome/gff_gene_count.txt...\n"
  grep -c -w gene CLEANED_GFF/cleaned_*gff > gff_gene_count.txt

  rm -r GFF_chrOK
  cd ..

}


count_genes_per_chr() {
  local gff=$1
  local out_file=$2
  cut -f1 $gff | sort | uniq -c > $out_file
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
  #cat ../CLEANED_GFF/cleaned_*.gff > ${new_prefix}.gff
  python3 $concatAndRmRepeatGenes --gff_list $gff_list --prefix ${PWD}/../CLEANED_GFF/cleaned_ --output ${new_prefix}.gff
  python3 ${scripts}/Extract_sequences_from_genome.py -g ${new_prefix}.gff -f ../LRR_ANNOT/${new_prefix}.fasta -o ${new_prefix}_proteins.fasta -t FSprot
  singularity run $LRRprofiler_sif --in_proteome ${new_prefix}_proteins.fasta --name ${new_prefix}_LRRprofiler_output
  source $scripts/../bin/lib_gff_comment.sh
  add_family_info ${new_prefix}.gff $PWD/Res_${new_prefix}_LRRprofiler_output/Res_step3/LRR_classification.csv $PWD/${new_prefix}_LRR_info.gff

  grep -w gene ${new_prefix}_LRR_info.gff | grep -v "Fam=" | cut -f2 -d"=" | cut -f1 -d";" > not_LRR_genes.list
  remove_genes_from_gff ${new_prefix}_LRR_info.gff not_LRR_genes.list ../LRR_ANNOT/${new_prefix}_LRR.gff

  ## Compare the number of genes per chromosome before and after removing non-LRR genes
  echo -e "... Counting genes before and after removing non-LRR genes >> see 02_build_exp_LRRome/LRR_ANNOT/genes_per_chr_[before/after]RemovingNonLRR.txt ...\n"
  count_genes_per_chr ${new_prefix}_LRR_info.gff ../LRR_ANNOT/genes_per_chr_beforeRemovingNonLRR.txt
  count_genes_per_chr ../LRR_ANNOT/${new_prefix}_LRR.gff ../LRR_ANNOT/genes_per_chr_afterRemovingNonLRR.txt

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

concat_gff() {         # >> Fills 04_final_GFF
  echo -e "\n4/ MERGING EXPERTISED AND INITIAL GFF FILES...\n"
  local init_prefix=$1
  local exp_prefix=$2
  local initial_gff=$3
  cd 04_final_GFF

  cp $initial_gff ${init_prefix}_${exp_prefix}_LRR.gff
  cat ../02_build_exp_LRRome/LRR_ANNOT/${exp_prefix}_LRR.gff >> ${init_prefix}_${exp_prefix}_LRR.gff
  rmCR ${init_prefix}_${exp_prefix}_LRR.gff

  cd ..
}

create_info_locus() {         # >> Fills 04_final_GFF
  echo -e "\n5/ CREATING THE INFO LOCUS FILE...\n"
  local gff=$1
  local exp_prefix=$2
  local init_prefix=$3


  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^ID=/) {
        Gene_ID=gensub("ID=", "", "g", $i)
      }
      if ($i ~ /^Fam=/) {
        Fam=gensub("Fam=", "", "g", $i)
      }
      if ($i ~ /^Class=/) {
        Class=gensub("Class=", "", "g", $i)
      }
      if ($i ~ /^Gene-Class:/) {
        Class=gensub("Gene-Class:", "", "g", $i)
      }
    }
    print Gene_ID, Fam, Class
  }}' $gff >> 04_final_GFF/${init_prefix}_${exp_prefix}_info_locus.txt

}

compute_gene_stats(){         # >> Fills 04_final_GFF
  echo -e "\n6/ COMPUTING THE FINAL GFF GENE STATS...\n"
  local gff=$1
  local stats_file=$2
  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
    Chr=$1
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^ID=/) {
        Sp=gensub("ID=", "", "g", $i)
        Sp=gensub(/_.*/, "", "g", Sp)
      }
      if ($i ~ /^Fam=/) {
        Fam=gensub("Fam=", "", "g", $i)
      }
      if ($i ~ /^Class=/) {
        Class=gensub("Class=", "", "g", $i)
      }
      if ($i ~ /^Gene-Class:/) {
        Class=gensub("Gene-Class:", "", "g", $i)
      }
    }
    print Sp, Chr, Fam, Class
  }}' $gff | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $5, $1}' > $stats_file
}

write_infos() {
  local new_prefix=$1
  date > LRRome_incremental_build_infos.txt

  echo -e "\n Input files:" >> LRRome_incremental_build_infos.txt
  echo "- Initial LRRome: "$initial_LRRome >> LRRome_incremental_build_infos.txt
  echo "- Expertised gff files in: 01_gff_EXP" >> LRRome_incremental_build_infos.txt

  echo -e "\n Output files:" >> LRRome_incremental_build_infos.txt
  echo "- Exp reference fasta file: 02_build_exp_LRRome/LRR_ANNOT/${new_prefix}.fasta" >> LRRome_incremental_build_infos.txt
  echo "- Exp gff in 02_build_exp_LRRome/LRR_ANNOT/${new_prefix}_LRR.gff" >> LRRome_incremental_build_infos.txt
  echo "- List of removed non-LRR genes: 02_build_exp_LRRome/LRRprofile/not_LRR_genes.list" >> LRRome_incremental_build_infos.txt
  echo "- Exp LRRome: 02_build_exp_LRRome/LRRome" >> LRRome_incremental_build_infos.txt
  echo "- New LRRome: 03_LRRome" >> LRRome_incremental_build_infos.txt
  echo "- New gff file: 04_final_GFF/${init_prefix}_${exp_prefix}_LRR.gff" >> LRRome_incremental_build_infos.txt
  echo "- Locus info file: 04_final_GFF/${init_prefix}_${exp_prefix}_info_locus.txt" >> LRRome_incremental_build_infos.txt

  echo -e "\n Stats files:" >> LRRome_incremental_build_infos.txt
  echo "- Gene numbers per chromosome before removing non-LRR genes: 02_build_exp_LRRome/LRR_ANNOT/genes_per_chr_beforeRemovingNonLRR.txt" >> LRRome_incremental_build_infos.txt
  echo "- Gene numbers per chromosome after removing non-LRR genes: 02_build_exp_LRRome/LRR_ANNOT/genes_per_chr_afterRemovingNonLRR.txt" >> LRRome_incremental_build_infos.txt
  echo "- Final LRRome gene stats in tidy format: 04_final_GFF/${init_prefix}_${exp_prefix}_gene_stats.tsv" >> LRRome_incremental_build_infos.txt

  echo -e "\n\nSee LRRome_incremental_build_infos.txt for run information.\n"
}

## ----------------------------------- MAIN ------------------------------------------------ ##


mkdir -p 02_build_exp_LRRome 03_LRRome 04_final_GFF

clean_gff $(realpath 01_gff_EXP) $exp_prefix

build_exp_LRRome $exp_prefix $exp_ref_folder

exp_LRRome=$(realpath 02_build_exp_LRRome/LRRome)
merge_LRRome $initial_LRRome $exp_LRRome 03_LRRome

concat_gff $init_prefix $exp_prefix $initial_LRR_gff

exp_LRR_gff=$(realpath 02_build_exp_LRRome/LRR_ANNOT/${exp_prefix}_LRR.gff)
create_info_locus 04_final_GFF/${init_prefix}_${exp_prefix}_LRR.gff $exp_prefix $init_prefix

compute_gene_stats 04_final_GFF/${init_prefix}_${exp_prefix}_LRR.gff 04_final_GFF/${init_prefix}_${exp_prefix}_gene_stats.tsv

write_infos $exp_prefix

# sbatch --partition=agap_normal --mem=20G --wrap="/home/girodollej/scratch/2024_LRR/03_scripts/06_LRRTRANSFER_LAUNCH_PREP/build_incremental_LRRome.sh"
