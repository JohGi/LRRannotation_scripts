#!/bin/bash

set -euo pipefail

if ! command -v python3 &> /dev/null; then
  module load python/packages/3.8.2
fi
if ! command -v singularity &> /dev/null; then
  module load singularity/3.6.3
fi
if ! command -v agat &> /dev/null; then
  module load AGAT/1.2.0-singularity
fi


## ------------------------------- FUNCTIONS --------------------------------------------- ##

require_variables() {
  local missing=0
  for varname in "$@"; do
    if [[ -z "${!varname+x}" || -z "${!varname}" ]]; then
      echo "Error: variable '$varname' is not set or is empty." >&2
      missing=1
    fi
  done
  return $missing
}

rmCR() {
  local file=$1
  if grep -q $'\r' ${file}; then
    sed -i 's/\r$//g' $file
    sed -i 's/\r/\n/g' $file
  fi
}

clean_gff() {             # >> Creates 02_build_exp_LRRome/CLEANED_GFF
  echo -e "\n1/ CLEANING EXP GFF FILES...\n"
  local gff_list=$1
  local cleaner_chr_prefix=$2
  local scripts=$3
  local out_dir=$4
  mkdir -p $out_dir && cd $out_dir
  mkdir -p CLEANED_GFF GFF_chrOK

  ## Change the chromosome names in the gff files if needed
  echo -e "... Fixing chromosome names...\n"
  for chr_num in $(seq 1 7); do
    for genome in A B; do
      chr=${chr_num}${genome}
      for gff in $(cat ${gff_list}); do
        if [[ $gff == *hr${chr}* ]]; then
          awk -v chr=$chr 'BEGIN{OFS=FS="\t"}{$1="Chr"chr; print $0}' $gff > GFF_chrOK/$(basename $gff)
        fi 
      done
    done
  done

  ## Clean the gff files with gff_cleaner
  echo -e "... Running gff_cleaner.py...\n"
  rm -f clean_gff.list
  for init_gff in $(cat ${gff_list}); do
    gff_name=$(basename $init_gff)
    gff=GFF_chrOK/${gff_name}
    echo -e python3 ${scripts}/VR/gff_cleaner.py -a -p $cleaner_chr_prefix -g $gff -o CLEANED_GFF/cleaned_${gff_name}"\n"
    python3 ${scripts}/VR/gff_cleaner.py -a -p $cleaner_chr_prefix -g $gff -o CLEANED_GFF/cleaned_${gff_name}
    echo $(realpath CLEANED_GFF/cleaned_${gff_name}) >> clean_gff.list
  done | grep -v "^WARNING: incompatible bounds" > gff_cleaner.out

  ## Count the genes in the gff files
  echo -e "... Counting genes in each gff in ${out_dir}/gff_gene_count.txt...\n"
  grep -c -w gene CLEANED_GFF/cleaned_*gff > gff_gene_count.txt

  rm -r GFF_chrOK
  cd ..

}


count_genes_per_chr() {
  local gff=$1
  local out_file=$2
  awk '{if ($3 == "gene"){print $1}}' $gff | sort | uniq -c > $out_file
}

build_exp_LRRome_multiGFF(){
  local new_prefix=$1
  local input_fasta=$2
  local gff_list=$(realpath $3)
  local GMT_sif=$4
  local seq_type=$5
  local concatAndRmRepeatGenes=$6
  local LRRprofiler_sif=$7
  local scripts=$8
  local out_dir=$9

  mkdir -p ${out_dir}/LRRprofile
  cd ${out_dir}/LRRprofile
  python3 $concatAndRmRepeatGenes --gff_list $gff_list --output ${new_prefix}.gff

  cd ../..

  build_exp_LRRome $new_prefix $input_fasta ${new_prefix}.gff ${GMT_sif} ${seq_type} ${LRRprofiler_sif} ${scripts} ${out_dir}
}

build_exp_LRRome() {        # >> Creates 02_build_exp_LRRome/LRR_ANNOT, 02_build_exp_LRRome/LRRprofile and 02_build_exp_LRRome/LRRome
  echo -e "\n2/ BUILDING EXPERTISED GENES LRROME...\n"
  local new_prefix=$1
  local input_fasta=$2
  local input_gff=$3
  local GMT_sif=$4
  local seq_type=$5
  local LRRprofiler_sif=$6
  local scripts=$7
  local out_dir=$8
  mkdir -p ${out_dir}
  cd ${out_dir}

  ## Detect and classify LRR containing genes with LRRprofiler and remove non LRR genes (write the final gff in LRR_ANNOT)
  echo -e "... Running LRRprofiler and removing non-LRR genes >> see final gff in ${out_dir}/LRR_ANNOT/"${new_prefix}"_LRR.gff...\n"

  mkdir -p LRRprofile LRR_ANNOT
  cd LRRprofile

  if [[ ${seq_type} == "FSprot" ]] ; then
    echo -e singularity exec $GMT_sif python3 ${scripts}/Extract_sequences_from_genome.py -g ${input_gff} -f ${input_fasta} -o ${new_prefix}_proteins.fasta -t FSprot"\n"
    singularity exec $GMT_sif python3 ${scripts}/Extract_sequences_from_genome.py -g ${input_gff} -f ${input_fasta} -o ${new_prefix}_proteins.fasta -t FSprot
  elif [[ ${seq_type} == "prot" ]] ; then
    echo -e agat_sp_extract_sequences.pl -g ${input_gff} -f ${input_fasta} -t cds -p -o ${new_prefix}_proteins.fasta"\n"
    agat_sp_extract_sequences.pl -g ${input_gff} -f ${input_fasta} -t cds -p -o ${new_prefix}_proteins.fasta
    sed -i -e '/^>/s/>.*gene=/>/' -e '/^>/s/ seq_id=.*//' ${new_prefix}_proteins.fasta
  else
    echo "ERROR: Unknow sequence type (${seq_type})"
    exit 1
  fi

  echo -e singularity run $LRRprofiler_sif --in_proteome ${new_prefix}_proteins.fasta --name ${new_prefix}_LRRprofiler_output --dev"\n"
  singularity run $LRRprofiler_sif --in_proteome ${new_prefix}_proteins.fasta --name ${new_prefix}_LRRprofiler_output --dev

  source $scripts/../bin/lib_gff_comment.sh
  echo -e add_family_info ${input_gff} $PWD/Res_${new_prefix}_LRRprofiler_output/Res_step3/LRR_classification.csv $PWD/${new_prefix}_LRR_info.gff"\n"
  add_family_info ${input_gff} $PWD/Res_${new_prefix}_LRRprofiler_output/Res_step3/LRR_classification.csv $PWD/${new_prefix}_LRR_info.gff

  awk '{if ($3 == "gene"){print $0}}' ${new_prefix}_LRR_info.gff | grep -v "Fam=" | cut -f2 -d"=" | cut -f1 -d";" > not_LRR_genes.list
  awk '{if ($3 == "gene"){print $0}}' ${new_prefix}_LRR_info.gff | grep -E "Fam=UC|Fam=F-box" | cut -f2 -d"=" | cut -f1 -d";" >> not_LRR_genes.list

  echo -e remove_genes_from_gff ${new_prefix}_LRR_info.gff not_LRR_genes.list ../LRR_ANNOT/${new_prefix}_LRR.gff"\n"
  remove_genes_from_gff ${new_prefix}_LRR_info.gff not_LRR_genes.list ../LRR_ANNOT/${new_prefix}_LRR.gff

  ## Compare the number of genes per chromosome before and after removing non-LRR genes
  echo -e "... Counting genes before and after removing non-LRR genes >> see ${out_dir}/LRR_ANNOT/genes_per_chr_[before/after]RemovingNonLRR.txt ...\n"
  count_genes_per_chr ${new_prefix}_LRR_info.gff ../LRR_ANNOT/genes_per_chr_beforeRemovingNonLRR.txt
  count_genes_per_chr ../LRR_ANNOT/${new_prefix}_LRR.gff ../LRR_ANNOT/genes_per_chr_afterRemovingNonLRR.txt

  cd ..

  ## Build LRRome (in a LRRome folder)
  echo -e "\n... Running create_LRRome.sh in ${out_dir}/LRRome...\n"
  echo -e singularity exec $GMT_sif ${scripts}/../bin/create_LRRome.sh ${input_fasta} $PWD/LRR_ANNOT/${new_prefix}_LRR.gff $PWD NULL ${scripts}"\n"
  singularity exec $GMT_sif ${scripts}/../bin/create_LRRome.sh ${input_fasta} $PWD/LRR_ANNOT/${new_prefix}_LRR.gff $PWD NULL ${scripts}

  cd ..
}


merge_LRRome() {         # >> Fills 03_LRRome
  echo -e "\n3/ MERGING EXPERTISED AND INITIAL LRROMES...\n"
  local LRRome1=$1
  local LRRome2=$2
  local LRRome_new=$3

  mkdir -p $LRRome_new && cd $LRRome_new
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
  local initial_gff=$1
  local exp_LRRome_gff=$2
  local final_gff_out_dir=$3
  local final_gff_name=$4
  
  mkdir -p ${final_gff_out_dir}
  cp $initial_gff ${final_gff_out_dir}/${final_gff_name}
  cat ${exp_LRRome_gff} >> ${final_gff_out_dir}/${final_gff_name}
  rmCR ${final_gff_out_dir}/${final_gff_name}

}

create_info_locus() {         # >> Fills 04_final_GFF
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
  }}' $gff > $output

}

compute_gene_stats(){         # >> Fills 04_final_GFF
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
  }}' $gff | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $5, $1}' > $stats_file
}

compute_gene_stats_generic_gff(){
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
  }}' $gff | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2, $3, $1}' > $stats_file
}


write_infos() {
  local initial_LRRome=$1
  local exp_gff_list=$2
  local init_prefix=$3
  local exp_prefix=$4
  local exp_LRRome_dir=$5
  local LRRome_dir=$6
  local final_gff_dir=$7

  date > LRRome_incremental_build_infos.txt

  echo -e "\n Input files:" >> LRRome_incremental_build_infos.txt
  echo "- Initial LRRome: "$initial_LRRome >> LRRome_incremental_build_infos.txt
  echo "- Expertised gff files in: ${exp_gff_list}" >> LRRome_incremental_build_infos.txt

  echo -e "\n Output files:" >> LRRome_incremental_build_infos.txt
  echo "- Exp reference fasta file: ${exp_LRRome_dir}/LRR_ANNOT/${exp_prefix}.fasta" >> LRRome_incremental_build_infos.txt
  echo "- Exp gff in ${exp_LRRome_dir}/LRR_ANNOT/${exp_prefix}_LRR.gff" >> LRRome_incremental_build_infos.txt
  echo "- List of removed non-LRR genes: ${exp_LRRome_dir}/LRRprofile/not_LRR_genes.list" >> LRRome_incremental_build_infos.txt
  echo "- Exp LRRome: ${exp_LRRome_dir}/LRRome" >> LRRome_incremental_build_infos.txt
  echo "- New LRRome: ${LRRome_dir}" >> LRRome_incremental_build_infos.txt
  echo "- New gff file: ${final_gff_dir}/${init_prefix}_${exp_prefix}_LRR.gff" >> LRRome_incremental_build_infos.txt
  echo "- Locus info file: ${final_gff_dir}/${init_prefix}_${exp_prefix}_info_locus.txt" >> LRRome_incremental_build_infos.txt

  echo -e "\n Stats files:" >> LRRome_incremental_build_infos.txt
  echo "- Gene numbers per chromosome before removing non-LRR genes: ${exp_LRRome_dir}/LRR_ANNOT/genes_per_chr_beforeRemovingNonLRR.txt" >> LRRome_incremental_build_infos.txt
  echo "- Gene numbers per chromosome after removing non-LRR genes: ${exp_LRRome_dir}/LRR_ANNOT/genes_per_chr_afterRemovingNonLRR.txt" >> LRRome_incremental_build_infos.txt
  echo "- Final LRRome gene stats in tidy format: ${final_gff_dir}/${init_prefix}_${exp_prefix}_gene_stats.tsv" >> LRRome_incremental_build_infos.txt

  echo -e "\n\nSee LRRome_incremental_build_infos.txt for run information.\n"
}
