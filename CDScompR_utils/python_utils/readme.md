## compare_annots.py

Run with:  
``` 
sif=/mnt/c/Users/girodolle/Documents/Apptainer/Apptainer_python3.11/v0.2/env_python3.11_libs.sif
python_utils_dir=/mnt/c/Users/girodolle/Documents/GitHub/LRRannotation_scripts/CDScompR_utils/python_utils

apptainer run --bind /mnt/c/Users/girodolle/Documents $sif ${python_utils_dir}/scripts/compare_annots.py --ref_gff ${python_utils_dir}/example_data/ref_test.gff --pred_gff ${python_utils_dir}/example_data/pred_test.gff --cdscompr_csv ${python_utils_dir}/example_data/pred_test_clean.csv -o overlaps_test.tsv

```

Run tests with:  
```
apptainer exec --bind /mnt/c/Users/girodolle/Documents $sif pytest ${python_utils_dir}/tests/test_overlap_group.py -v
```
