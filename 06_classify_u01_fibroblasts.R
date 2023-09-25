# Armon Azizi
#
# Cluster U01 Fibroblast samples on fibroblast genes from lambrechts.

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)

setwd("~/Bioinformatics/sc_prognosis/")

U01<-fread("U01_RNAseq/GSE111907_processed_data.txt")
U01<-column_to_rownames(U01, "Name")

sample_info<-as.data.frame(fread("U01_RNAseq/U01_sample_info.txt"))
adeno_samples<-sample_info$Sample[sample_info$Pathology=="ADENO"]

U01<-U01[,adeno_samples]

fibroblast_sample_names<-grep("CD10",colnames(U01), value=TRUE)

U01_fibroblasts<-U01[,fibroblast_sample_names]
duplicate_samples<-c("T28_CD10","T30_CD10")
U01_fibroblasts<-U01_fibroblasts[,!colnames(U01_fibroblasts)%in%duplicate_samples]
#U01_fibroblasts<-log2(U01_fibroblasts+1)


diff_genes<-as.data.frame(fread("outputs/LUAD_fibroblasts/lambrechts_luad_adv_vs_fav_fibroblasts_nmf_diff_exp.tsv"))
adv_genes<-diff_genes$Gene[rev(order(diff_genes$avg_lfc))][1:100]
fav_genes<-diff_genes$Gene[order(diff_genes$avg_lfc)][1:100]

scores<-colSums(U01_fibroblasts[adv_genes,])-colSums(U01_fibroblasts[fav_genes,])
scores<-scores[rev(order(scores))]

# monte carlo to find significance of classifications
p_vals<-data.frame()
for(s in colnames(U01_fibroblasts)){
  message(s)
  score<-sum(U01_fibroblasts[adv_genes,s])-sum(U01_fibroblasts[fav_genes,s])
  simulations<-c()
  for(i in 1:100){
    values<-sample(U01_fibroblasts[c(adv_genes, fav_genes),s])
    simulations<-c(simulations, sum(values[1:100])-sum(values[101:200]))
  }
  
  p <- min(sum(simulations<score)/100,sum(simulations>score)/100)
  
  if(p==0){
    p<-1/100
  }
  
  p<-p*2
  
  p_vals[s,"p_val"]<-p
  p_vals[s,"score"]<-score
  
}

adv_fibroblast_samples<-intersect(names(scores[scores>0]),rownames(p_vals)[p_vals$p_val<0.05])
fav_fibroblast_samples<-intersect(names(scores[scores<0]),rownames(p_vals)[p_vals$p_val<0.05])

# barplot(scores[rev(order(scores))])

heatmap_input<-log2(U01_fibroblasts[c(adv_genes,fav_genes),c(adv_fibroblast_samples,fav_fibroblast_samples)]+1)
heatmap_input<-log2(U01_fibroblasts[c(adv_genes,fav_genes),rownames(p_vals)[rev(order(p_vals$score))]]+1)

# gene_order<-c(adv_genes[rev(order(log2(rowMeans(heatmap_input[adv_genes,adv_fibroblast_samples])/rowMeans(heatmap_input[adv_genes,fav_fibroblast_samples]))))],
#               fav_genes[order(log2(rowMeans(heatmap_input[fav_genes,adv_fibroblast_samples])/rowMeans(heatmap_input[fav_genes,fav_fibroblast_samples])))])
# 
# heatmap_input<-heatmap_input[gene_order,]

#heatmap_input<-t(scale(t(heatmap_input)))

hmcs<-heatmap_input
hmcs[is.na(hmcs)] <- 0
for(i in 1:1){
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

heatmap_input<-hmcs

Heatmap(heatmap_input, 
        cluster_rows = FALSE,
        cluster_columns=FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        col = viridis(10),
        heatmap_legend_param = list(title = "Expression"))
