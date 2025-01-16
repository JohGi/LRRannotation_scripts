#commit: https://github.com/JohGi/LRRannotation_scripts/commit/182c6d2067ff95300412ac05c47858ffaecf4133

gff=/storage/replicated/cirad/projects/GE2POP/2023_LRR/SVEVO3_LRR_ANNOT_2024_07_12/02_OUTPUT_2024_08_03/annot_best.gff
zones=zones.tsv

./extractGenesFromZones.sh $gff $zones OUTPUTS