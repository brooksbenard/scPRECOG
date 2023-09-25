# Armon Azizi
# combine all data

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)

setwd("~/Bioinformatics/sc_prognosis/")

source("scripts/sc_prognostic_scoring.R")

#source("~/Bioinformatics/single_cell/scripts/sc_precog_functions.R")

dataset_info<-read.xlsx("dataset_info/dataset_info.xlsx")

all_counts<-data.frame()
all_phenotypes<-data.frame()

for(i in 1:nrow(dataset_info)){
  
  dataset<-dataset_info[i,"dataset"]
  precog_label<-dataset_info[i,"precog_label"]
  tcga_label<-dataset_info[i,"tcga_label"]
  counts_file<-dataset_info[i,"counts_file"]
  cell_column<-dataset_info[i,"cell_column"]
  celltype_column<-dataset_info[i,"celltype_column"]
  sample_column<-dataset_info[i,"sample_column"]
  
  message(paste("analyzing:", dataset))
  
  counts<-as.data.frame(fread(paste0("standardized_data/",dataset,"/",counts_file)))
  counts<-column_to_rownames(counts, colnames(counts)[1])
  phenotypes<-as.data.frame(fread(paste0("standardized_data/",dataset,"/phenotypes.tsv")))
  
  colnames(counts)<-paste(as.character(i),colnames(counts),sep="_")
  phenotypes[,cell_column]<-paste(as.character(i),phenotypes[,cell_column],sep="_")
  
  phenotypes_mod<-phenotypes[,c(cell_column,celltype_column)]
  if(!is.na(sample_column)){
    phenotypes_mod$sample<-phenotypes[,sample_column]
  }else{
    phenotypes_mod$sample<-rep(NA,nrow(phenotypes_mod))
  }
  
  colnames(phenotypes_mod)<-c("ID","CellType","Sample")
  phenotypes_mod$dataset<-rep(dataset,nrow(phenotypes_mod))
  
  if(i==1){
    all_counts<-counts
  }else{
    shared_genes<-intersect(rownames(all_counts),rownames(counts))
    all_counts<-all_counts[shared_genes,]
    counts<-counts[shared_genes,]
    all_counts<-cbind(all_counts,counts)
  }
  
  all_phenotypes<-rbind(all_phenotypes,phenotypes_mod)
}

all_cells<-intersect(colnames(all_counts),all_phenotypes$ID)
all_counts<-all_counts[,all_cells]
all_phenotypes<-all_phenotypes[match(all_cells,all_phenotypes$ID),]

saveRDS(all_counts,"combined_datasets/normalized_counts.rds")
saveRDS(all_phenotypes,"combined_datasets/phenotypes.rds")
