### Lancement Muse sur fichiers tests
```srun --partition=agap_normal --pty bash```  
```git clone git@github.com:JohGi/LRRannotation_scripts.git```  
```cd LRRannotation_scripts/findOverlaps```  
```module load python/packages/3.8.2```  

```mkdir results && python ./findOverlaps.py --ref_gff test_files/ref.gff --alt_gff test_files/alt.gff --overlaps_output results/results_overlaps.txt --groups_output results/results_groups.txt --overlap_thr 0.2 --overreach_thr 0.5 --verbose --show_all_genes --show_all_groups```

&nbsp;
### Lancement tests unitaires
```pytest -vv```  
