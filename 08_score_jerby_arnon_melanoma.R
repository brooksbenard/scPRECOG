# Armon Azizi
#
# Analyze jerby arnon melanoma single cell data using prognostic scoring code

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(NMF)

setwd("~/Bioinformatics/sc_prognosis/")

source("scripts/sc_prognostic_scoring.R")
source("scripts/sc_precog_functions.R")

precog_data<-read.table("PRECOG/PRECOG_metaZ_table.txt",sep="\t",header=TRUE,row.names=1)

precog_z<-precog_data[,"Melanoma_Metastasis",drop=FALSE]

marker_genes<-readRDS("outputs/melanoma/jerby_arnon_celltype_marker_genes.rds")
malignant_genes<-marker_genes$gene[marker_genes$p_val_adj<0.05&marker_genes$cluster=="Mal"&marker_genes$avg_logFC>0]
precog_z_malignant<-precog_z[intersect(malignant_genes, rownames(precog_z)),,drop=FALSE]


counts<-as.data.frame(fread("standardized_data/jerby_arnon_melanoma/normalized_counts.tsv"))
counts<-column_to_rownames(counts, colnames(counts)[1])
phenotypes<-as.data.frame(fread("standardized_data/jerby_arnon_melanoma/phenotypes.tsv"))

sample_info<-as.data.frame(fread("standardized_data/jerby_arnon_melanoma/sample_info.txt"))
sample_info<-sample_info[sample_info$`Lesion type`=="metastasis",]

# add celltype count info
# celltype_counts<-as.data.frame.matrix(table(phenotypes[phenotypes$sample%in%sample_info$Sample,c(2,3)]))
# sample_info<-cbind(sample_info, celltype_counts[match(sample_info$Sample, rownames(celltype_counts)),])

#### Celltype prognostic score comparison ####

# metastatic samples
met_samples<-sample_info$Sample
met_cells<-phenotypes$ID[phenotypes$sample%in%met_samples]
met_counts<-counts[,met_cells]
met_pheno<-phenotypes[phenotypes$ID%in%met_cells,]

met_scores<-run_prognostic_scoring(expression = met_counts, prognostic_data = precog_z, cor_method = "pearson")

colors<-get_palette(palette="ucscgb", k=20)

plot_input<-data.frame(score=met_scores$prognostic_score,
                       celltype=met_pheno$CellType[match(rownames(met_scores),met_pheno$ID)], stringsAsFactors = FALSE)
plot_input<-plot_input[plot_input$celltype!="?",]
plot_input<-plot_input[complete.cases(plot_input),]
plot_input$celltype<-as.character(plot_input$celltype)
plot_input$celltype[plot_input$celltype=="Cancer-Associated Fibroblast"]<-"CAF"
order<-aggregate(plot_input[,"score"], list(plot_input$celltype), mean)
order<-rev(order[order(order[,2]),][,1])
ggplot(plot_input, aes(x=celltype, y=score, fill=celltype)) +
  geom_jitter(position=position_jitter(0.1), size=0.1) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.size=-1) +
  scale_fill_manual(values = colors[1:length(unique(plot_input$celltype))]) +
  xlab(NULL)+
  ylab("Prognostic Score (PRECOG)") +
  theme_bw() +
  scale_x_discrete(limits=order) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=10),
        legend.position = "none") 
