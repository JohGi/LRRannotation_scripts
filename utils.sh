function create_fake_infoLocus() {
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

get_gene_class() {
  local gff=$1
  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
    gsub(/comment=/, "", $0)
    ID = ""
    Class = ""
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^ID=/) {
        ID=gensub("ID=", "", "g", $i)
      }
      if ($i ~ /^Class=/) {
        Class=gensub("Class=", "", "g", $i)
      }
      if ($i ~ /^Gene-Class:/) {
        Class=gensub("Gene-Class:", "", "g", $i)
      }
    }
    print ID, Class
  }}' $gff
}