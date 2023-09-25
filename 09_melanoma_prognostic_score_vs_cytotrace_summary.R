# Armon Azizi
#
#
# analyze malignant cells in melanoma and summarize cytotrace simialrity

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
source("scripts/cytotrace.R")

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
samples<-sample_info$Sample

cyto_similarity<-data.frame(sample=samples,
                            correlation=rep(NA,length(samples)),
                            p=rep(NA,length(samples)),
                            cells=rep(NA,length(samples)))

for(s in samples){
  
  message(paste("analyzing:", s))
  
  sample_cells<-phenotypes$ID[phenotypes$CellType=="Mal"&phenotypes$sample==s]
  sample_counts<-counts[,sample_cells]
  
  if(length(sample_cells)<10){next}
  
  prognostic_scores<-run_prognostic_scoring(expression = sample_counts, prognostic_data = precog_z, cor_method = "pearson")
  
  subsample_size<-ncol(sample_counts)-1
  cytotrace_scores<-cytoTRACE(sample_counts, enableFast = TRUE, subsamplesize = subsample_size, rank_space = FALSE)
  cytotrace_scores<-as.data.frame(cytotrace_scores$CytoTRACE)
  colnames(cytotrace_scores)<-c("cytotrace")
  
  scores<-data.frame(prognostic=prognostic_scores$prognostic_score[match(sample_cells, rownames(prognostic_scores))],
                     cytotrace=cytotrace_scores$cytotrace[match(sample_cells, rownames(cytotrace_scores))])
  
  scores<-scores[scores$cytotrace>0,]
  
  # plot<-ggplot(scores, aes(x=prognostic,y=cytotrace)) +
  #   geom_point()
  # print(plot)
  
  res<-corr.test(scores$prognostic,scores$cytotrace,method = "spearman")
  
  cyto_similarity[cyto_similarity$sample==s,"correlation"]<-res$r
  cyto_similarity[cyto_similarity$sample==s,"p"]<-res$p
  cyto_similarity[cyto_similarity$sample==s,"cells"]<-length(sample_cells)
}


plot_input<-cbind(sample_info[,c("Age","Sex","Lesion type","lymph_node_status")],cyto_similarity)
plot_input$significant<-rep(T, nrow(plot_input))
plot_input$significant[plot_input$p>0.05]<-F

ggplot(plot_input, aes(x=lymph_node_status, y=correlation, color=significant)) +
  geom_jitter(size=3, width=0.2) +
  xlab("") +
  ylab("Correlation Between Prognostic Score\nAnd cytoTRACE in Malignant cells") +
  theme_bw()

ggplot(plot_input, aes(x=correlation, y=-log10(p), color=significant, size=cells)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("black","darkred")) +
  xlab("Correlation Between Prognostic Score\nAnd cytoTRACE in Malignant cells") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text= element_text(color="black", size=10))

stem_samples<-cyto_similarity$sample[cyto_similarity$correlation>0&cyto_similarity$p<0.05]
stem_samples<-as.character(stem_samples[!is.na(stem_samples)])
saveRDS(stem_samples, "outputs/melanoma/jerby_arnon_samples_with_stem_heterogeneity.rds")
