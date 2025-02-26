# Usage: Rscript plot_scoreComp.R <CDScompR_outputs_merge.tsv> <HsapChr1LRRt> <HsapChr1LO> <OUTPUTS/03_plots>


library(ggplot2)
library(reshape)
library(scales)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
merged_CDScompR_output<-args[1]
alt1_name<-args[2]
alt2_name<-args[3]
outdir<-args[4]

# Output files
id_score_plot<-paste0(outdir,"/id_score_plot.png")
id_score_distribution<-paste0(outdir,"/id_score_distribution.png")

#alt1_name<-"HsapChr1LRRt"
#alt2_name<-"HsapChr1LO"
#setwd("~/2024_LRR/03_scripts/13_LAUNCH_LRRt_PRIMATES/SCRIPTS")
#merged_CDScompR_output<-"CDScompR_outputs_merge.tsv"

# Expected column names
alt1_locus_colname<-paste0("Alt_locus_", alt1_name)
alt2_locus_colname<-paste0("Alt_locus_", alt2_name)
idscore_alt1_colname<-paste0("Identity_Score_", alt1_name)
idscore_alt2_colname<-paste0("Identity_Score_", alt2_name)

# Import data
data<-read.table(merged_CDScompR_output, header=T)
data_subset<-subset(data, data[[alt1_locus_colname]] != "~" & data[[alt2_locus_colname]] != "~")[,c("gene_id", idscore_alt1_colname, idscore_alt2_colname)]
nb_genes<-nrow(data_subset)


# Identity scores plot
png(id_score_plot, width=750, height=750)
plot(data_subset[[idscore_alt1_colname]]~data_subset[[idscore_alt2_colname]], data=data_subset, pch=16, col=alpha("red", 0.5), cex=2, xlab=idscore_alt2_colname, ylab=idscore_alt1_colname)
abline(a=0, b=1)
dev.off()

# Identity scores distribution
colnames(data_subset)[2:3]<-c(alt1_name, alt2_name)
data_subset_for_ggplot<-melt(data_subset, id="gene_id")
colnames(data_subset_for_ggplot)[2:3]<-c("Tool","Identity_Score")
data_subset_for_ggplot$Identity_Score<-as.numeric(data_subset_for_ggplot$Identity_Score)
png(id_score_distribution, width=750) 
ggplot(data_subset_for_ggplot, aes(x = Identity_Score, fill = Tool)) +
geom_histogram(position = "dodge", bins = 20, alpha = 0.6) +
labs(title = paste0("Distributions of CDScompR identity scores (%) when compared to the reference annotation (",nb_genes," genes)"),
     x = "Identity Score",
     y = "Frequency") +
theme_minimal() +
scale_fill_manual(values = c("#FF5733", "#33C1FF"))
dev.off()
