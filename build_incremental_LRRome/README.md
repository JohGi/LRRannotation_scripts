# `build_incremental_LRRome.sh`

Incrementally builds a new expertised LRRome by integrating curated GFF annotations into an existing LRRome reference.  
This pipeline includes cleaning GFF files, running LRRprofiler, removing non-LRR genes, merging annotation sets, and generating updated statistics and metadata.

---

## Overview

This script allows you to:

- Integrate new expertised LRR genes into a reference LRRome.
- Optionally provide extra GFFs to enrich LRRprofiler's domain detection (without adding these genes to the final LRRome).
- Run LRRprofiler to detect and validate LRR domains based on curated genes.
- Merge initial and new annotations after filtering non-LRR genes.
- Generate updated GFFs, statistics, and sequence files.
- Trace the build with comprehensive log and metadata output.

---

## Usage
```
./build_incremental_LRRome.sh path/to/config.sh

```
Use debug mode for verbose output:
```
./build_incremental_LRRome.sh path/to/config.sh --debug
```
---

## Input Requirements

The script expects a configuration bash file defining the following variables:

```bash
GMT_DIR=...               # Path to GeneModelTransfer repo
GMT_SIF=...               # Path to GMT Singularity image
LRRPROFILER_SIF=...       # Path to LRRprofiler Singularity image
GFF_LIST=...              # TXT file listing expertised GFFs (most to least recent)
EXTRA_GFF_LIST=...        # Optional GFFs for profiling (excluded from new LRRome)
INITIAL_LRROME=...        # Path to the original LRRome directory
INITIAL_LRR_GFF=...       # Path to GFF file associated with initial LRRome
EXP_REF_GENOME=...        # Path to reference genome
EXP_PREFIX=...            # Prefix for expertised output
INIT_PREFIX=...           # Prefix for initial LRRome
SEQ_TYPE=prot|FSprot      # Type of sequence to extract (CDS or frameshift-aware)
```
> See BILRRome_config.sh for an example.

#### About `EXTRA_GFF_LIST`

If provided, `EXTRA_GFF_LIST` should point to a text file listing additional expertised GFF files.

These GFFs are:
- **Used only for LRRprofiler detection**, to increase the diversity of the LRR profile database.
- **Not included in the final LRRome**: genes from these files are excluded from the output sequences and GFF.

Set this variable to `""` (empty string) or omit it if not needed.

> **Example use case**  
> You have an initial LRRome built from **two different species** (e.g. *Oryza sativa* and *Triticum turgidum*), and you now want to:
> - Add new manually expertised LRR genes **only** for *T. turgidum* (via `GFF_LIST`),
> - But still want LRRprofiler to use **all previously validated *T. turgidum*** genes (from the initial LRRome) to build better HMM profiles.
>
> In this case:
> - Provide the new *T. turgidum* GFFs in `GFF_LIST`,
> - Extract *only the *T. turgidum* GFFs from the initial LRRome* and list them in `EXTRA_GFF_LIST`.
>
> This ensures that:
> - LRRprofiler uses a richer training set for that species,
> - Only the newly validated expertised genes are actually added to the updated LRRome.

---

## Outputs
The script will generate:
```
📂 01_cleaned_gff/
 └── Cleaned expertised input GFFs (NEW and EXTRA)

📂 02_build_exp_LRRome/
 ├── 01_LRRprofiler/
 |   └── LRRprofiler output
 |   └── Intermediary GFFs files
 |   └── Intermediary protein FASTA file
 |   └── Non-LRR genes list
 ├── 02_LRR_gff/
 |   └── GFF of the expertised genes (without potential extra genes) after removing non-LRR genes
 |   └── Gene count files
 └── 03_LRRome/
     └── Expertised LRRome sequences and annotation (without potential extra genes)

📂 03_LRRome/
 └── Final merged LRRome (initial + expertised genes, without potential extra genes)

📂 04_final_GFF/
 ├── Final merged GFF file
 ├── Locus info file
 └── Gene stats (tidy table)

📄 LRRome_incremental_build_infos.txt (log file summarizing output file paths)
```

---
## Project Structure

```
.
├── build_incremental_LRRome.sh   # Main launcher script  
├── lib_LRRome.sh                 # Library with all modular bash functions  
├── concatAndRmRepeatGenes.py     # Python script to concatenate GFF files  
├── BILRRome_config.sh            # Example config file (user provided)  
└── ...
```
---
## 🛠 Requirements

- Bash
- Singularity
- Python 3
- [GeneModelTransfer](https://github.com/ranwez/GeneModelTransfer/tree/master) (via `.sif`)
- [LRRprofiler](https://github.com/ranwez/LRRprofiler) (via `.sif`)
- AGAT (`agat_sp_extract_sequences.pl`)

> Modules will be automatically loaded if not available in the environment (only works on Muse HPC).

---

