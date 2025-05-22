function create_fake_infoLocus(){
  local gff=$1
  local out_file=$2
  awk -F '[ \t;]' 'BEGIN{OFS="\t"}{if ($3 == "gene"){
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^ID=/) {
          Gene_ID=gensub("ID=", "", "g", $i)
        }
      }
      print Gene_ID, "UC", "Canonical"
    }}' $gff > $out_file
}