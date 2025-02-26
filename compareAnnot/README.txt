#Lancement sur Muse :

profile="CONFIG/CLUSTER_PROFILE_SLURM"
config_file="CONFIG/config.yaml"
smk_Pipeline="/home/girodollej/scratch/2024_LRR/03_scripts/13_LAUNCH_LRRt_PRIMATES/compare_annotations.smk"
sbatch --partition=agap_normal --wrap="module load snakemake/7.32.4-conda ; module load singularity/3.6.3 ; snakemake --profile $profile --snakefile $smk_Pipeline --configfile $config_file --keep-going --singularity-args \"--bind ${HOME}\""
