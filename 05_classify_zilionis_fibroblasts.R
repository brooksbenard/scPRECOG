# Armon Azizi
#
#
# classify zilionis fibroblasts using the lambrechts differentially expressed genes.

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)

setwd("~/Bioinformatics/sc_prognosis/")

source("scripts/sc_prognostic_scoring.R")
source("scripts/sc_precog_functions.R")

dataset_info<-read.xlsx("dataset_info/dataset_info.xlsx")
scored_data<-read.xlsx("dataset_info/celltype_assignments.xlsx")

i<-14

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

fibroblast_cells<-phenotypes[,cell_column][phenotypes[,celltype_column]=="Fibroblasts"]

fibroblast_counts<-counts[,fibroblast_cells]

diff_genes<-as.data.frame(fread("outputs/LUAD_fibroblasts/lambrechts_luad_adv_vs_fav_fibroblasts_nmf_diff_exp.tsv"))
adv_genes<-diff_genes$Gene[rev(order(diff_genes$avg_lfc))][1:100]
fav_genes<-diff_genes$Gene[order(diff_genes$avg_lfc)][1:100]

#### Identify bad fibroblasts using lambrechts signature ####
samples<-unique(phenotypes[phenotypes$CellType=="Fibroblasts",][,sample_column])

bad_fibroblasts_from_gene_sig<-c()

for(s in samples){
  
  message(paste("analyzing:", s))
  
  fibroblast_sample_cells<-phenotypes[,cell_column][phenotypes[,celltype_column]=="Fibroblasts"&phenotypes[,sample_column]==s]
  fibroblast_sample_counts<-fibroblast_counts[,fibroblast_sample_cells]
  
  scores<-colSums(fibroblast_sample_counts[intersect(adv_genes,rownames(fibroblast_sample_counts)),])-colSums(fibroblast_sample_counts[intersect(fav_genes,rownames(fibroblast_sample_counts)),])
  
  #hist(scores)
  
  adv_fibroblasts<-names(scores)[scores>=0]
  
  bad_fibroblasts_from_gene_sig<-c(bad_fibroblasts_from_gene_sig,adv_fibroblasts)
}

#### Subtype specific expression signatures on lambrechts genes ####
fibroblast_1<-bad_fibroblasts_from_gene_sig
fibroblast_2<-setdiff(fibroblast_cells,bad_fibroblasts_from_gene_sig)


scaled_fibroblast_counts<-log2(fibroblast_counts[,phenotypes$ID[phenotypes$patientID%in%samples&phenotypes$CellType=="Fibroblasts"]]+1)
for(s in samples){
  message(s)
  sample_cells<-intersect(colnames(scaled_fibroblast_counts),phenotypes$ID[phenotypes$PatientNumber==s])
  scaled_data<-as.data.frame(t(scale(t(scaled_fibroblast_counts[,sample_cells,drop=FALSE]))))
  scaled_fibroblast_counts<-cbind(scaled_fibroblast_counts[,!colnames(scaled_fibroblast_counts)%in%sample_cells],scaled_data)
}
scaled_fibroblast_counts[is.na(scaled_fibroblast_counts)]<-0

fibroblast_1_genes<-adv_genes
fibroblast_2_genes<-fav_genes

prognostic_fibroblast_counts<-scaled_fibroblast_counts[c(fibroblast_1_genes,fibroblast_2_genes),c(fibroblast_1,fibroblast_2)]

hmcs<-prognostic_fibroblast_counts
hmcs[is.na(hmcs)] <- 0
for(i in 1:10){
  print(i)
  hmcs <- t(apply(hmcs, 1, function(x) x/sqrt(sum(x^2))))
  hmcs[is.na(hmcs)] <- 0
  hmcs <- hmcs - apply(hmcs,1,mean)
  hmcs[which(hmcs<0)] <- 0
}
hmcs <- hmcs - apply(hmcs,1,mean)
thres <- .025
hmcs[(hmcs>thres)] <- thres
hmcs[(hmcs<= -thres)] <- -thres

#sliding window
if(1==1){
  hmcs[(hmcs<= 0)] <- 0
  w <- 1
  for(i in 1:nrow(hmcs))
  {
    d <- sapply(c(1:(ncol(hmcs)-w)), function(j) mean(hmcs[i,c(j:(j+w))]))
    if(i==1) smooth <- d
    else{smooth <- rbind(smooth,d)}
  }
  rownames(smooth) = rownames(hmcs)
  colnames(smooth) = colnames(hmcs)[1:(ncol(hmcs)-w)]
  hmcs <- smooth
  hmcs <- (hmcs - sapply(1:nrow(hmcs),function(i) mean(hmcs[i,])))
  thres <- .005
  hmcs[which(hmcs>thres)] <- thres
  hmcs[which(hmcs<= -thres)] <- -thres
}
heatmap_input<-hmcs

if(TRUE){
  rownames(phenotypes)<-NULL
  annotations<-column_to_rownames(phenotypes[,c("ID","CellType","patientID")],"ID")
  colnames(annotations)<-c("CellType","sample")
  annotations<-annotations[annotations$sample%in%samples,]
  annotations<-annotations[annotations$CellType=="Fibroblasts",]
  annotations<-annotations[,"sample",drop=FALSE]
  annotations$CellType<-rep("Favorable Fibroblasts",nrow(annotations))
  annotations$CellType[rownames(annotations)%in%fibroblast_1]<-"Adverse Fibroblasts"
  annotations$sample<-as.factor(annotations$sample)
  annotations$CellType<-as.factor(annotations$CellType)
  annotations<-annotations[colnames(heatmap_input),,drop=FALSE]
  celltype_colors<-get_palette("ucscgb",k=length(levels(annotations$CellType)))
  names(celltype_colors)<-levels(annotations$CellType)
  sample_colors<-get_palette("Set1",k=length(levels(annotations$sample)))
  names(sample_colors)<-levels(annotations$sample)
  anno_colors<-list(CellType=celltype_colors,sample=sample_colors) 
}

h<-Heatmap(heatmap_input, 
           cluster_rows = FALSE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           show_column_names = FALSE,
           top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_legend_param = list(direction = "horizontal", nrow=2)),
           col = viridis(10),
           #col = colorRampPalette(c("lightblue","white","darkred"))(10),
           heatmap_legend_param = list(title = "Scaled Expression", nrow=1))
ComplexHeatmap::draw(h, annotation_legend_side = "bottom")
