Script pour préparer un nouveau LRRome (dans 03_LRRome) à partir de :
- une série de gff expertisés dans 01_gff_EXP
- le génome de référence des gff (un dossier contenant les fasta des chromosomes individuels nommés Chr*.fasta)
- un LRRome déjà fait (ex : riz)



--------------------------- LANCEMENT ---------------------------
sbatch --partition=agap_normal --mem=20G --wrap="./build_incremental_LRRome.sh"


# seff du dernier lancement :
Job ID: 25929535
Cluster: cluster
User/Group: girodollej/agap
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:30:41
CPU Efficiency: 92.61% of 00:33:08 core-walltime
Job Wall-clock time: 00:33:08
Memory Utilized: 9.57 GB
Memory Efficiency: 47.83% of 20.00 GB

NB: 
- Les warnings affichés au moment du create_LRRome.sh concernant les tailles non multiples de 3 vont dans le stdout (slurm...out si run dans un job).
Ils viennent de biopython ("Partial codon, len(sequence) not a multiple of three.") et de Extract_sequences_from_genome.py (donne le nom du gène + la séquence).
- Les warnings/errors de gff_cleaner (notamment sur les CDS/mRNA sans parent) sont écrits dans 02_build_exp_LRRome/gff_cleaner.out
-------------------------------------------------------------------



------------------------- IDEE GENERALE ---------------------------
Si on a déjà un LRRome (par ex riz + blé dur) et qu'on veut rajouter des gènes expertisés de blé, on peut repartir du LRRome du riz tel quel, mais il faut refaire celui du blé, puis fusionner les deux. C'est à cause de l'étape de LRRprofiler qui identifie les gènes (LRR ou non et quel type de LRR) en prenant en compte des caractéristiques du jeu de données (et donc du génome). Il faut donc l'appeler séparément sur tous les gènes du riz et sur tous les gènes du blé.

Les étapes pour reconstruire le LRRome de blé (build_exp_LRRome dans le script) :
- concaténer les fastas des chromosomes en 1 seul fasta (1)
- corriger les seqid (noms de chromosomes) dans les gff exp + faire tourner gff_cleaner dessus
- concaténer les gff en 1 seul
- extraction des séquences protéiques correspondant au gff (Extract_sequences_from_genome.py -FSprot)
- identification et classification des gènes LRR à partir des séquences protéiques (LRRprofiler)
- filtre des gènes non-LRR > gff final du blé (2)
- construction du LRRome avec create_LRRome.sh à partir du fasta (1) et du gff (2)
-------------------------------------------------------------------



---------------------------- SORTIES ------------------------------
NB : Le fichier REF_proteins.fasta sert d'entrée à LRRtransfert.

1/ le fasta du génome des gènes exp [02_build_exp_LRRome/LRR_ANNOT/${prefix}.fasta]

2/ le gff des gènes exp [02_build_exp_LRRome/LRR_ANNOT/${prefix}_LRR.gff]

3/ Un LRRome (créé par create_LRRome.sh) [03_LRRome]
= Un dossier qui contient :
- les fastas entiers des cDNA, exons et protéines à la racine
- un dossier par fasta entier qui contient les fastas individuels de chaque séquence du fasta

03_LRRome
|-- REF_EXONS
|   |-- DWSvevo3July_Chr1A_0006928950_mrna_1_CDS_1
|   |-- DWSvevo3July_Chr1A_0006928950_mrna_1_CDS_2
|   |-- OSJnip_Chr12_27344060:cds_1
|   |-- OSJnip_Chr12_27344060:cds_2
|   |-- ...
|-- REF_PEP
|   |-- DWSvevo3July_Chr1A_0006928950
|   |-- DWSvevo3July_Chr1A_0007534996
|   |-- OSJnip_Chr12_27122673
|   |-- OSJnip_Chr12_27344060
|   |-- ...
|-- REF_cDNA
|   |-- DWSvevo3July_Chr1A_0006928950
|   |-- DWSvevo3July_Chr1A_0007534996
|   |-- OSJnip_Chr12_27122673
|   |-- OSJnip_Chr12_27344060
|   |-- ...
|-- REF_cDNA.fasta
|-- REF_exons.fasta
`-- REF_proteins.fasta
-------------------------------------------------------------------

