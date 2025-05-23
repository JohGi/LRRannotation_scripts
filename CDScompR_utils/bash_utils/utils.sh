print_common_genes() {
  local csv=$1
  awk -F ',' '{
    if (NR == 1){
      print $0
    } else if ($3 != "~" && $4 != "~") {
      print $0
    }
  }' $csv
}

print_specific_genes() {
  local csv=$1
  local dataset=$2
  if [[ "$dataset" == "REF" ]]; then
    i=3
    j=4
  elif [[ "$dataset" == "ALT" ]]; then
    i=4
    j=3
  else
    echo "Error: unknown dataset value '$dataset'" >&2
    exit 1
  fi

  awk -F ',' -v i=$i -v j=$j '{
    if (NR == 1){
      print $0
    } else if ($i != "~" && $j == "~") {
      print $0
    }
  }' $csv
}
