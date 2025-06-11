#!/bin/bash

set -euo pipefail

if ! command -v python3 &>/dev/null; then
  module load python/packages/3.8.2
fi
if ! command -v singularity &>/dev/null; then
  module load singularity/3.6.3
fi
if ! command -v agat &>/dev/null; then
  module load AGAT/1.2.0-singularity
fi

## ------------------------------- FUNCTIONS --------------------------------------------- ##

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

rmCR() {
  local file=$1
  if grep -q $'\r' ${file}; then
    sed -i 's/\r$//g' $file
    sed -i 's/\r/\n/g' $file
  fi
}

get_gene_ids() {
  local gff=$1
  
  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^ID=/) {
        gene_ID=gensub("ID=", "", "g", $i)
        break
      }
    }
    print gene_ID
  }}' $gff
}

get_gene_ids_from_gff_list(){
  local gff_list=$1

  for gff in $(cat $gff_list) ; do
    get_gene_ids $gff
  done
}

clean_gff() {
  local gff_list=$1
  local cleaner_chr_prefix=$2
  local GMT_dir=$3
  local out_dir=$4

  mkdir -p ${out_dir}/GFF_chrOK

  ## Change the chromosome names in the gff files if needed
  echo -e "... Fixing chromosome names...\n"
  for chr_num in $(seq 1 7); do
    for genome in A B; do
      chr=${chr_num}${genome}
      for gff in $(cat ${gff_list}); do
        if [[ $gff == *hr${chr}* ]]; then
          awk -v chr=$chr 'BEGIN{OFS=FS="\t"}{$1="Chr"chr; print $0}' $gff >${out_dir}/GFF_chrOK/$(basename $gff)
        fi
      done
    done
  done

  ## Clean the gff files with gff_cleaner
  echo -e "... Running gff_cleaner.py...\n"
  rm -f ${out_dir}/clean_gff.list
  for init_gff in $(cat ${gff_list}); do
    gff_name=$(basename $init_gff)
    gff=${out_dir}/GFF_chrOK/${gff_name}
    echo -e python3 ${GMT_dir}/SCRIPT/VR/gff_cleaner.py -a -p $cleaner_chr_prefix -g $gff -o ${out_dir}/cleaned_${gff_name}"\n"
    python3 ${GMT_dir}/SCRIPT/VR/gff_cleaner.py -a -p $cleaner_chr_prefix -g $gff -o ${out_dir}/cleaned_${gff_name}
    echo $(realpath ${out_dir}/cleaned_${gff_name}) >>${out_dir}/clean_gff.list
  done | grep -v "^WARNING: incompatible bounds" >${out_dir}/gff_cleaner.out

  ## Count the genes in the gff files
  echo -e "... Counting genes in each gff in ${out_dir}/gff_gene_count.txt...\n"
  grep -c -w gene ${out_dir}/cleaned_*gff >${out_dir}/gff_gene_count.txt

  rm -r ${out_dir}/GFF_chrOK
}

clean_gff_for_LRRprofiler(){
  local cleaner_chr_prefix=$1
  local GMT_dir=$2
  local new_gff_list=$3
  local out_dir=$4
  local extra_gff_list="${5:-}"

  echo -e "\n1/ CLEANING INPUT GFF FILES...\n"

  clean_gff $new_gff_list $cleaner_chr_prefix ${GMT_dir} ${out_dir}/NEW_GFF
  cp ${out_dir}/NEW_GFF/clean_gff.list ${out_dir}/
  if [[ -n "${extra_gff_list}" ]] ; then 
    clean_gff $extra_gff_list $cleaner_chr_prefix ${GMT_dir} ${out_dir}/EXTRA_GFF
    cat ${out_dir}/clean_gff.list ${out_dir}/EXTRA_GFF/clean_gff.list > ${out_dir}/clean_gff_tmp.list && mv ${out_dir}/clean_gff_tmp.list ${out_dir}/clean_gff.list
    get_gene_ids_from_gff_list ${out_dir}/EXTRA_GFF/clean_gff.list > ${out_dir}/extra_gene_IDs.list
  fi
}

count_genes_per_chr() {
  local gff=$1
  local out_file=$2
  awk '{if ($3 == "gene"){print $1}}' $gff | sort | uniq -c >$out_file
}

extract_prot_sequences(){
    local input_gff=$1
    local input_fasta=$2
    local seq_type=$3
    local GMT_sif=$4
    local GMT_dir=$5
    local prot_sequences_output=$6

    
    if [[ ${seq_type} == "FSprot" ]]; then
      echo -e singularity exec ${GMT_sif} python3 ${GMT_dir}/SCRIPT/Extract_sequences_from_genome.py -g ${input_gff} -f ${input_fasta} -o ${prot_sequences_output} -t FSprot --no_FS_codon"\n"
      singularity exec ${GMT_sif} python3 ${GMT_dir}/SCRIPT/Extract_sequences_from_genome.py -g ${input_gff} -f ${input_fasta} -o ${prot_sequences_output} -t FSprot --no_FS_codon
    elif [[ ${seq_type} == "prot" ]]; then
      echo -e agat_sp_extract_sequences.pl -g ${input_gff} -f ${input_fasta} -t cds -p -o ${prot_sequences_output}"\n"
      agat_sp_extract_sequences.pl -g ${input_gff} -f ${input_fasta} -t cds -p -o ${prot_sequences_output}
      sed -i -e '/^>/s/>.*gene=/>/' -e '/^>/s/ seq_id=.*//' ${prot_sequences_output}
    else
      echo "ERROR: Unknow sequence type (${seq_type})"
      exit 1
    fi
  }

run_LRRprofiler(){
  local input_prot_sequences=$(realpath $1)
  local LRRprofiler_sif=$(realpath $2)
  local out_name=$3
  local out_dir=$4

  cd ${out_dir}
  echo -e singularity run ${LRRprofiler_sif} --in_proteome ${input_prot_sequences} --name ${out_name} --dev"\n"
  singularity run ${LRRprofiler_sif} --in_proteome ${input_prot_sequences} --name ${out_name} --dev
  cd - > /dev/null
}

remove_non_LRR_genes(){
  local init_gff=$1
  local LRRprofiler_res_dir=$2
  local GMT_dir=$3
  local out_gff_prefix=$4
  local LRR_info_out_dir=$5
  local LRR_gff_out_dir=$6
  local extra_genes_to_rm_list="${7:-}"

  local LRR_info_gff=${LRR_info_out_dir}/${out_gff_prefix}_LRR_info.gff

  source ${GMT_dir}/bin/lib_gff_comment.sh
  echo -e add_family_info ${init_gff} ${LRRprofiler_res_dir}/Res_step3/LRR_classification.csv ${LRR_info_gff}"\n"
  add_family_info ${init_gff} ${LRRprofiler_res_dir}/Res_step3/LRR_classification.csv ${LRR_info_gff} && rm -f __LRR_family.tmp

  if [[ -n "${extra_genes_to_rm_list}" ]] ; then
    local LRR_info_with_extra_genes_gff=${LRR_info_out_dir}/${out_gff_prefix}_with_extra_genes_LRR_info.gff
    mv ${LRR_info_gff} ${LRR_info_with_extra_genes_gff}
    echo "remove_genes_from_gff ${LRR_info_with_extra_genes_gff} $extra_genes_to_rm_list ${LRR_info_gff}"
    remove_genes_from_gff ${LRR_info_with_extra_genes_gff} $extra_genes_to_rm_list ${LRR_info_gff}
  fi

  local non_LRR_list=${LRR_info_out_dir}/not_LRR_genes.list 
  awk '{if ($3 == "gene"){print $0}}' ${LRR_info_gff} | grep -v "Fam=" | cut -f2 -d"=" | cut -f1 -d";" > ${non_LRR_list} || true
  awk '{if ($3 == "gene"){print $0}}' ${LRR_info_gff} | grep -E "Fam=UC|Fam=F-box" | cut -f2 -d"=" | cut -f1 -d";" >> ${non_LRR_list} || true

  echo -e remove_genes_from_gff ${LRR_info_gff} ${non_LRR_list} ${LRR_gff_out_dir}/${out_gff_prefix}_LRR.gff"\n"
  remove_genes_from_gff ${LRR_info_gff} ${non_LRR_list} ${LRR_gff_out_dir}/${out_gff_prefix}_LRR.gff
}

detect_and_rm_non_LRR_genes(){
  local input_gff=$1
  local input_fasta=$2
  local seq_type=$3
  local GMT_sif=$4
  local GMT_dir=$5
  local LRRprofiler_sif=$6
  local new_prefix=$7
  local LRRprofiler_out_dir=$8
  local LRR_gff_out_dir=$9
  local extra_genes_to_rm_list="${10:-}"

  echo -e "... Running LRRprofiler and removing non-LRR genes >> see final gff in ${LRR_gff_out_dir}/${new_prefix}_LRR.gff...\n"

  if [[ -n "$extra_genes_to_rm_list" ]] ; then
    local prot_sequences=${LRRprofiler_out_dir}/${new_prefix}_with_extra_genes_proteins.fasta
  else
    local prot_sequences=${LRRprofiler_out_dir}/${new_prefix}_proteins.fasta
  fi
  
  extract_prot_sequences ${input_gff} ${input_fasta} ${seq_type} ${GMT_sif} ${GMT_dir} ${prot_sequences}

  run_LRRprofiler ${prot_sequences} ${LRRprofiler_sif} ${new_prefix}_LRRprofiler_output ${LRRprofiler_out_dir}

  remove_non_LRR_genes ${input_gff} ${LRRprofiler_out_dir}/Res_${new_prefix}_LRRprofiler_output ${GMT_dir} ${new_prefix} ${LRRprofiler_out_dir} ${LRR_gff_out_dir} ${extra_genes_to_rm_list}
}

compute_stats_before_after_removing_non_LRR(){
  local gff_before=$1
  local gff_after=$2
  local out_dir=$3

  echo -e "... Counting genes before and after removing non-LRR genes >> see ${out_dir}/genes_per_chr_[before/after]RemovingNonLRR.txt ...\n"
  count_genes_per_chr ${gff_before} ${out_dir}/genes_per_chr_beforeRemovingNonLRR.txt
  count_genes_per_chr ${gff_after} ${out_dir}/genes_per_chr_afterRemovingNonLRR.txt
}

run_create_LRRome(){
  local input_gff=$(realpath $1)
  local input_fasta=$(realpath $2)
  local GMT_sif=$3
  local GMT_dir=$4
  local LRRome_out_dir=$(realpath $5)

  local out_dir=$(dirname "${LRRome_out_dir%/}")
  echo -e "\n... Running create_LRRome.sh in ${LRRome_out_dir}...\n"
  echo -e singularity exec ${GMT_sif} ${GMT_dir}/bin/create_LRRome.sh ${input_fasta} ${input_gff} ${out_dir} NULL ${GMT_dir}/SCRIPT"\n"
  singularity exec ${GMT_sif} ${GMT_dir}/bin/create_LRRome.sh ${input_fasta} ${input_gff} ${out_dir} NULL ${GMT_dir}/SCRIPT
  mv ${out_dir}/LRRome ${LRRome_out_dir}
}

build_exp_LRRome() {
  echo -e "\n2/ BUILDING EXPERTISED GENES LRROME...\n"
  local new_prefix=$1
  local input_fasta=$2
  local input_gff=$3
  local GMT_sif=$4
  local seq_type=$5
  local LRRprofiler_sif=$6
  local GMT_dir=$7
  local out_dir=$8
  local extra_gene_list=${9:-}

  local LRRprofiler_out_dir=${out_dir}/01_LRRprofiler
  local LRR_gff_out_dir=${out_dir}/02_LRR_gff
  mkdir -p ${LRRprofiler_out_dir} ${LRR_gff_out_dir}

  detect_and_rm_non_LRR_genes ${input_gff} ${input_fasta} ${seq_type} ${GMT_sif} ${GMT_dir} ${LRRprofiler_sif} ${new_prefix} ${LRRprofiler_out_dir} ${LRR_gff_out_dir} ${extra_gene_list}
  final_gff=${LRR_gff_out_dir}/${new_prefix}_LRR.gff

  compute_stats_before_after_removing_non_LRR ${LRRprofiler_out_dir}/${new_prefix}_LRR_info.gff ${final_gff} ${LRR_gff_out_dir}
  
  run_create_LRRome ${final_gff} ${input_fasta} ${GMT_sif} ${GMT_dir} ${out_dir}/03_LRRome 
}

build_exp_LRRome_multiGFF() {
  local new_prefix=$1
  local input_fasta=$2
  local gff_list=$3
  local GMT_sif=$4
  local seq_type=$5
  local concatAndRmRepeatGenes=$6
  local LRRprofiler_sif=$7
  local GMT_dir=$8
  local out_dir=$9
  local extra_gene_list=${10:-}

  local LRRprofiler_out_dir=${out_dir}/01_LRRprofiler
  mkdir -p ${LRRprofiler_out_dir}

  if [[ -n "$extra_gene_list" ]] ; then
    local concat_gff=${LRRprofiler_out_dir}/${new_prefix}_with_extra_genes.gff
  else
    local concat_gff=${LRRprofiler_out_dir}/${new_prefix}.gff
  fi
  python3 $concatAndRmRepeatGenes --gff_list $gff_list --output ${concat_gff}

  build_exp_LRRome $new_prefix $input_fasta ${concat_gff} ${GMT_sif} ${seq_type} ${LRRprofiler_sif} ${GMT_dir} ${out_dir} ${extra_gene_list}
}

merge_LRRome() {
  echo -e "\n3/ MERGING EXPERTISED AND INITIAL LRROMES...\n"
  local LRRome1=$(realpath $1)
  local LRRome2=$(realpath $2)
  local LRRome_new=$(realpath $3)

  mkdir -p ${LRRome_new}/REF_EXONS ${LRRome_new}/REF_PEP ${LRRome_new}/REF_cDNA

  local type
  for type in EXONS PEP cDNA; do
    cp ${LRRome1}/REF_${type}/* ${LRRome_new}/REF_${type}/
    cp ${LRRome2}/REF_${type}/* ${LRRome_new}/REF_${type}/
  done

  for type in exons proteins cDNA; do
    cp ${LRRome1}/REF_${type}.fasta ${LRRome_new}/ && chmod u+w ${LRRome_new}/REF_${type}.fasta
    cat ${LRRome2}/REF_${type}.fasta >>${LRRome_new}/REF_${type}.fasta
  done

  echo -e "The new LRRome was built in: "$(realpath ${LRRome_new})"\n"
}

concat_gff() {
  echo -e "\n4/ MERGING EXPERTISED AND INITIAL GFF FILES...\n"
  local initial_gff=$1
  local exp_LRRome_gff=$2
  local final_gff_out_dir=$3
  local final_gff_name=$4

  mkdir -p ${final_gff_out_dir}
  cp $initial_gff ${final_gff_out_dir}/${final_gff_name}
  chmod u+w ${final_gff_out_dir}/${final_gff_name}
  cat ${exp_LRRome_gff} >>${final_gff_out_dir}/${final_gff_name}
  rmCR ${final_gff_out_dir}/${final_gff_name}

}

create_info_locus() {
  echo -e "\n5/ CREATING THE INFO LOCUS FILE...\n"
  local gff=$1
  local output=$2

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
  }}' $gff >$output

}

compute_gene_stats() {
  echo -e "\n6/ COMPUTING THE FINAL GFF GENE STATS...\n"
  local gff=$1
  local stats_file=$2
  awk -F '[ \t;]|comment=' 'BEGIN{
    OFS="\t"
    fam_keys = "^(Fam:|Fam=|family:)"
    class_keys = "^(Class=|Gene-Class:|class:)"
  }
  {if ($3 == "gene"){
    Fam="Non-LRR"
    Chr=$1
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^ID=/) {
        Sp=gensub("ID=", "", "g", $i)
        Sp=gensub(/_.*/, "", "g", Sp)
      }
      if (match($i, fam_keys)) {
        Fam = gensub(fam_keys, "", 1, $i)
      }
      if (match($i, class_keys)) {
        Class = gensub(class_keys, "", 1, $i)
      }
    }
    print Sp, Chr, Fam, Class
  }}' $gff | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $5, $1}' >$stats_file
}

compute_gene_stats_generic_gff() {
  local gff=$1
  local stats_file=$2
  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
    Fam="Non-LRR"
    Chr=$1
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^Fam=/) {
        Fam=gensub("Fam=", "", "g", $i)
      }
    }
    print Chr, Fam
  }}' $gff | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' >$stats_file
}

write_infos() {
  local initial_LRRome=$1
  local exp_gff_list=$2
  local init_prefix=$3
  local exp_prefix=$4
  local build_exp_LRRome_dir=$5
  local LRRome_dir=$6
  local final_gff_dir=$7
  local extra_gene_list=${8:-}

  local LRRprofiler_dir=${build_exp_LRRome_dir}/01_LRRprofiler
  local LRR_gff_dir=${build_exp_LRRome_dir}/02_LRR_gff
  local exp_LRRome_dir=${build_exp_LRRome_dir}/03_LRRome

  date >LRRome_incremental_build_infos.txt

  echo -e "\n Input files:" >>LRRome_incremental_build_infos.txt
  echo "- Initial LRRome: "$initial_LRRome >>LRRome_incremental_build_infos.txt
  echo "- Expertised gff files listed in: ${exp_gff_list}" >>LRRome_incremental_build_infos.txt

  echo -e "\n Output files:" >>LRRome_incremental_build_infos.txt
  echo "- Exp gff in ${LRR_gff_dir}/${exp_prefix}_LRR.gff" >>LRRome_incremental_build_infos.txt
  echo "- List of removed non-LRR genes: ${LRRprofiler_dir}/not_LRR_genes.list" >>LRRome_incremental_build_infos.txt
  echo "- Exp LRRome: ${exp_LRRome_dir}" >>LRRome_incremental_build_infos.txt
  echo "- New final LRRome: ${LRRome_dir}" >>LRRome_incremental_build_infos.txt
  echo "- New final gff file: ${final_gff_dir}/${init_prefix}_${exp_prefix}_LRR.gff" >>LRRome_incremental_build_infos.txt
  echo "- Associated locus info file: ${final_gff_dir}/${init_prefix}_${exp_prefix}_info_locus.txt" >>LRRome_incremental_build_infos.txt

  echo -e "\n Stats files:" >>LRRome_incremental_build_infos.txt
  echo "- Gene numbers per chromosome before removing non-LRR genes: ${LRR_gff_dir}/genes_per_chr_beforeRemovingNonLRR.txt" >>LRRome_incremental_build_infos.txt
  echo "- Gene numbers per chromosome after removing non-LRR genes: ${LRR_gff_dir}/genes_per_chr_afterRemovingNonLRR.txt" >>LRRome_incremental_build_infos.txt
  echo "- Final LRRome gene stats in tidy format: ${final_gff_dir}/${init_prefix}_${exp_prefix}_gene_stats.tsv" >>LRRome_incremental_build_infos.txt

  echo -e "\n\nSee LRRome_incremental_build_infos.txt for run information.\n"
}
