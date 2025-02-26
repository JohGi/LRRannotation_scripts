#!/usr/bin/env python

# Workflow params
default_threads = 1 #this will be erased by the user's specifications for each rule in the profile yaml
singularity_image = "/storage/simple/users/girodollej/projects/GE2POP/2023_LRR/LRRtransfer_image/LRRtransfer.sif"

path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"

# Output directories
outDir = "OUTPUTS"
CDScompR_outDir = outDir+"/01_CDScompR"
mergeCompR_outDir = outDir+"/02_merge_compR"
plots_outDir = outDir+"/03_plots"


# Functions
def tsv2dict(tsv):
    d = {}
    with open(tsv, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                columns = line.strip().split("\t")
                key = columns[0]
                values = columns[1:]
                d[key] = values
    return d



# Arguments from config file
ref_gff = config["ref_gff"]
ref_name = config["ref_name"]

CDScompR = config["CDScompR_path"]
sort_gff = config["sortGFF_path"]
merge_compR = config["merge_path"]

gff_dict = tsv2dict(config["alt_gff_list"])




## WORKFLOW

# --------------------------------------------------------
rule all:
  input:
    plots_outDir+"/id_score_plot.png",
    plots_outDir+"/id_score_distribution.png"
# --------------------------------------------------------

rule sort_refGFF:
    input:
        ref_gff=ref_gff
    output:
        sorted_gff=temp(outDir+"/"+ref_name+"_sorted.gff")
    singularity:
        singularity_image
    threads:
        default_threads
    shell:
        "python {sort_gff} -g {input} -o {output}"


rule sort_altGFF:
    input:
        alt_gff=lambda wildcards: gff_dict[wildcards.base]
    output:
        sorted_gff=temp(outDir+"/{base}_sorted.gff")
    singularity:
        singularity_image
    threads:
        default_threads
    shell:
        "python {sort_gff} -g {input} -o {output}"


rule CDScompR:
    input:
        ref_gff=outDir+"/"+ref_name+"_sorted.gff",
        alt_gff=outDir+"/{base}_sorted.gff"
    output:
        CDScompR_outDir+"/"+ref_name+"_{base}/output.csv"
    envmodules:
        "python/packages/3.7.2"
    threads:
        default_threads
    shell:
        "python {CDScompR} --verbose --reference {input.ref_gff} --alternative {input.alt_gff} --out_dir {CDScompR_outDir}/{ref_name}_{wildcards.base};"
        "mv {CDScompR_outDir}/{ref_name}_{wildcards.base}/*.csv {output}"



rule create_csvNameCorresp:
    input:
        CDScompR_outDir+"/"+ref_name+"_{base}/output.csv"
    output:
        temp(mergeCompR_outDir+"/{base}.list")
    threads: 
        default_threads
    shell:
        """
        echo "{input}\t{wildcards.base}" > {output}
        """


rule merge_CDSCRoutputs:
    input:
        ref_gff=outDir+"/"+ref_name+"_sorted.gff",
        csv_files=expand(mergeCompR_outDir+"/{base}.list", base=gff_dict.keys())
    output:
        merge_output = mergeCompR_outDir+"/CDScompR_outputs_merge.tsv",
        csv_file=mergeCompR_outDir+"/csv_files.list"
    envmodules:
        "python/packages/3.7.2"
    threads: 
        default_threads
    shell:
        """
        cat {input.csv_files} > {output.csv_file}
        python {merge_compR} --comp_csv_list {output.csv_file} --ref_gff {input.ref_gff} --output {output.merge_output}
        """

rule plot_scores:
    input:
        mergeCompR_outDir+"/CDScompR_outputs_merge.tsv"
    output:
        plots_outDir+"/id_score_plot.png",
        plots_outDir+"/id_score_distribution.png"
    envmodules:
        "R/packages/4.3.1"
    threads:
        default_threads
    params:
        alt1_name=list(gff_dict.keys())[0],
        alt2_name=list(gff_dict.keys())[1]
    shell:
        "Rscript {scripts_dir}/plot_scores.R {input} {params.alt1_name} {params.alt2_name} {plots_outDir}"



