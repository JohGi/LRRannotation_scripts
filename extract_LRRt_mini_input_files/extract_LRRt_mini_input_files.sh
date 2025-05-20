#extract_LRRt_mini_input_files.sh <UTILS_DIR> <GENE_LIST> <BLAST> <TARGET_GENOME> <REF_GENOME> <REF_GFF> <REF_LOCUS_INFO> <MAX_TARGET_ZONES> <OUTDIR> <REF_FASTA_MARGIN> <TARGET_FASTA_MARGIN>
# E. g. when LRRtransfer has been run and you want to quickly re-run the transfer on a small subset of genes.
# Requires :
# a gene list file,
# the blast produced by LRRt,
# the target genome/ref genome/ref gff/ref locus info used as input for LRRt,
# the max number of target zones you want to retain,
# and the margins (bp) to add around sequences in the output ref and target fasta files

set -euo pipefail

UTILS_DIR=$1
GENE_LIST=$2
BLAST=$3
TARGET_GENOME=$4
REF_GENOME=$5
REF_GFF=$6
REF_LOCUS_INFO=$7
MAX_TARGET_ZONES=$8
OUTDIR=$9
REF_FASTA_MARGIN=${10}
TARGET_FASTA_MARGIN=${11}

source ${UTILS_DIR}/bash_utils/init_utils.sh

## Functions

# create_ref_fasta() {
#   local ref_gff=$1
#   local out_ref_fasta=$2

#   module_load bedtools
#   tmp_dir=$(create_and_return_tmp_dir)
#   awk '{if ($3 == "gene") {print $1,$4,$5}}' $ref_gff | sort -k1,1 -k2,2n >${tmp_dir}/ref_genome_zones.bed
#   chrs=$(awk '{print $1}' ${tmp_dir}/ref_genome_zones.bed | sort | uniq)
#   margin=30
#   for chr in $chrs; do
#     last_pos=$(grep -w $chr ${tmp_dir}/ref_genome_zones.bed | tail -1 | cut -d ' ' -f3)
#     last_pos_with_margin=$(($last_pos + $margin))
#     echo -e "${chr}\t0\t${last_pos_with_margin}" >>${tmp_dir}/final_ref_genome_zones.bed
#   done

#   bedtools getfasta -fi ${REF_GENOME} -bed ${tmp_dir}/final_ref_genome_zones.bed -fo ${out_ref_fasta}
#   sed -i '/^>/ s/:.*//' ${out_ref_fasta}

#   rm -r $tmp_dir

# }

check_gff_fasta_consistency() {
  local old_gff=$1
  local new_gff=$2
  local old_fasta=$3
  local new_fasta=$4

  log "Checking for gff and fasta inconsistency"
  module_load bedtools __
  tmp_dir=$(create_and_return_tmp_dir)

  gff_fasta_corresp=(
    "$old_gff $old_fasta old_ref"
    "$new_gff $new_fasta new_ref"
  )
  for elements in "${gff_fasta_corresp[@]}"; do
    read -r gff fasta name <<<"$elements"
    awk -F '[\t=;]' 'BEGIN{OFS="\t"}{print $1, $4-1, $5, $10}' $gff >${tmp_dir}/${name}.bed #Convert GFF to BED: GFF is 1-based with inclusive start and end, while BED is 0-based with inclusive start and exclusive end >> we substract 1 from start and leave end unchanged
    bedtools getfasta -fi $fasta -bed ${tmp_dir}/${name}.bed -name | sed 's/::.*//' >${tmp_dir}/${name}.fasta
  done

  if ! diff ${tmp_dir}/old_ref.fasta ${tmp_dir}/new_ref.fasta >/dev/null; then
    exit_with_error "Extracted sequences from fasta files don't match. Maybe there was an error when computing the new coordinates or extracting the corresponding sequences"
  fi
  rm -r $tmp_dir
}

create_ref_gff_and_fasta() {
  local extracted_gff=$1
  local out_ref_gff=$2
  local out_ref_genome=$3

  log "Creating the new reference gff and fasta files with a margin of ${REF_FASTA_MARGIN} bp of flanking sequence around each gene"
  module_load singularity __
  singularity exec $GMT_SINGULARITY_IMAGE python ${UTILS_DIR}/python_utils/scripts/gff_utils/extract_seqs_and_recalculate_coords.py -g ${extracted_gff} -f ${REF_GENOME} -m ${REF_FASTA_MARGIN} -og ${out_ref_gff} -of ${out_ref_genome}
  check_gff_fasta_consistency ${extracted_gff} ${out_ref_gff} ${REF_GENOME} ${out_ref_genome}
}

create_target_fasta() {
  local out_target_fasta=$1

  log "Creating the new target fasta file by sampling a maximum of ${MAX_TARGET_ZONES} zones, with a margin of ${TARGET_FASTA_MARGIN} bp of flanking sequence around each zone"
  module_load bedtools __
  tmp_dir=$(create_and_return_tmp_dir)
  grep -w -f ${GENE_LIST} ${BLAST} | cut -f2,7,8 | awk -v margin=$TARGET_FASTA_MARGIN 'BEGIN{OFS="\t"}{if ($2>$3){print $1, $3-margin, $2+margin} else {print $1,$2-margin,$3+margin}}' >${tmp_dir}/blast.bed
  bedtools sort -i ${tmp_dir}/blast.bed | bedtools merge | shuf -n $MAX_TARGET_ZONES >${tmp_dir}/blast_merged_sample.bed
  bedtools getfasta -fi ${TARGET_GENOME} -bed ${tmp_dir}/blast_merged_sample.bed -fo ${out_target_fasta}
  sed -i '/^>/ s/:/_/' ${out_target_fasta}

  rm -r $tmp_dir
}

create_dummy_infoLocus() {
  local gff=$1
  local out_file=$2
  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^ID=/) {
          Gene_ID=gensub("ID=", "", "g", $i)
        }
      }
      print Gene_ID, "UC", "Canonical"
    }}' $gff >$out_file
}

create_info_locus() {
  out_info_locus=$1
  out_ref_gff=$2

  if [[ ${REF_LOCUS_INFO} == "__" ]]; then
    log "Creating a dummy info locus file"
    create_dummy_infoLocus ${out_ref_gff} ${out_info_locus}
  else
    log "Creating the info locus file"
    grep -w -f ${GENE_LIST} ${REF_LOCUS_INFO} >${out_info_locus}
  fi
}

print_config() {
  local out_target_genome=$1
  local out_ref_genome=$2
  local out_ref_gff=$3
  local out_info_locus=$4

  echo -e "\nInformation for LRRtransfer config file:"
  echo -e "target_genome: \"$(realpath ${out_target_genome})\""
  echo -e "ref_genome: \"$(realpath ${out_ref_genome})\""
  echo -e "ref_gff: \"$(realpath ${out_ref_gff})\""
  echo -e "ref_locus_info: \"$(realpath ${out_info_locus})\""
}

## Main
mkdir -p $OUTDIR

# extract relevant genes from ref gff
module_load singularity __
singularity exec $GMT_SINGULARITY_IMAGE python $GMT/SCRIPT/CANDIDATE_LOCI/filter_gff_by_geneID.py -g ${REF_GFF} -l ${GENE_LIST} -o ${OUTDIR}/tmp_ref.gff

# make ref gff and ref fasta
create_ref_gff_and_fasta ${OUTDIR}/tmp_ref.gff ${OUTDIR}/ref.gff ${OUTDIR}/ref.fasta
rm ${OUTDIR}/tmp_ref.gff

# make target fasta
create_target_fasta ${OUTDIR}/target.fasta

# make info_locus
create_info_locus ${OUTDIR}/info_locus.tsv ${OUTDIR}/ref.gff

# print config
print_config ${OUTDIR}/target.fasta ${OUTDIR}/ref.fasta ${OUTDIR}/ref.gff ${OUTDIR}/info_locus.tsv
