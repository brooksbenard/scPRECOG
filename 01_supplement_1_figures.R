# Supplement 1 figures: comparing no exclusion of dropout to excluding dropout in correlation

library(openxlsx)
library(ggplot2)
library(ggpubr)

setwd("~/Bioinformatics/sc_prognosis/")

precog_data<-read.table("PRECOG/PRECOG_metaZ_table.txt",sep="\t",header=TRUE,row.names=1)

dataset_info<-read.xlsx("dataset_info/dataset_info.xlsx")

scored_data<-read.xlsx("dataset_info/celltype_assignments.xlsx")

all_counts<-readRDS("combined_datasets/normalized_counts.rds")
all_phenotypes<-readRDS("combined_datasets/phenotypes.rds")


cells_of_interest<-c("1_AACCGCGGTACCGTAT_15", "1_TGGGCGTTCGAGAGCA_19")


cell_of_interest<-all_counts[,cells_of_interest[2],drop=FALSE]

prognostic_data<-precog_data[,"Lung_cancer_ADENO",drop=FALSE]

expression_data<-cell_of_interest[intersect(rownames(cell_of_interest), rownames(prognostic_data)),,drop=FALSE]

detected_genes<-rownames(expression_data)[expression_data[,1]>0]

nonzeroed_scores<-data.frame(matrix(nrow=length(detected_genes),ncol=3))
colnames(nonzeroed_scores)<-c("prognostic_score","num_genes_detected","fraction_genes_detected")

for(i in 1:(length(detected_genes)-10)){
  
  expression_zscores<-expression_data
  zeroed_genes<-sample(detected_genes,i,replace=FALSE)
  expression_zscores[zeroed_genes,1]<-0
  #expression_zscores[expression_zscores==0]<-NA
  expression_zscores[,1]<-scale(expression_zscores[,1])
  prognostic_data_sub<-prognostic_data[rownames(expression_zscores),,drop=FALSE]
  
  exp_vector<-expression_zscores[,cells_of_interest[2],drop=FALSE]
  exp_vector<-exp_vector[!is.na(exp_vector[,cells_of_interest[2]]),,drop=FALSE]
  zscores<-prognostic_data_sub[rownames(exp_vector),,drop=FALSE]
  nonzeroed_scores[i,"prognostic_score"]<-stats::cor(zscores[,1],exp_vector[,1], method = "pearson")
  nonzeroed_scores[i,"num_genes_detected"]<-length(detected_genes)-i
  nonzeroed_scores[i,"fraction_genes_detected"]<-(length(detected_genes)-i)/nrow(expression_data)
}


ggplot(nonzeroed_scores, aes(x=num_genes_detected,y=prognostic_score)) +
  geom_point() +
  theme_bw() +
  xlab("Number of Genes Detected") +
  ylab("Prognostic Score") +
  stat_cor(method = "pearson", size=5) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(color="black", size=14))


zeroed_scores<-data.frame(matrix(nrow=length(detected_genes),ncol=3))
colnames(zeroed_scores)<-c("prognostic_score","num_genes_detected","fraction_genes_detected")

for(i in 1:(length(detected_genes)-10)){
  
  expression_zscores<-expression_data
  zeroed_genes<-sample(detected_genes,i,replace=FALSE)
  expression_zscores[zeroed_genes,1]<-0
  expression_zscores[expression_zscores==0]<-NA
  expression_zscores[,1]<-scale(expression_zscores[,1])
  prognostic_data_sub<-prognostic_data[rownames(expression_zscores),,drop=FALSE]
  
  exp_vector<-expression_zscores[,cells_of_interest[2],drop=FALSE]
  exp_vector<-exp_vector[!is.na(exp_vector[,cells_of_interest[2]]),,drop=FALSE]
  zscores<-prognostic_data_sub[rownames(exp_vector),,drop=FALSE]
  zeroed_scores[i,"prognostic_score"]<-stats::cor(zscores[,1],exp_vector[,1], method = "pearson")
  zeroed_scores[i,"num_genes_detected"]<-length(detected_genes)-i
  zeroed_scores[i,"fraction_genes_detected"]<-(length(detected_genes)-i)/nrow(expression_data)
}


ggplot(zeroed_scores, aes(x=num_genes_detected,y=prognostic_score)) +
  geom_point() +
  theme_bw() +
  xlab("Number of Genes Detected") +
  ylab("Prognostic Score") +
  stat_cor(method = "pearson", size=5) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color="black", size=12),
        axis.title = element_text(color="black", size=14))
