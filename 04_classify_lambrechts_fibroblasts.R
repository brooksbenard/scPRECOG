# Armon Azizi
# aazizi@stanford.edu
#
# Code to analyze fibroblast populations from lambrechts Lung adenocarcinoma samples
# Score each cell using PRECOG prognostic scores
# Cluster cells and do differential expression
# Do heatmap of differential genes and volcano plot
# Do pca plot of different features
# Generate cibersort signature matrix input
#
# downstream analysis not included: make cibersort signature matrix input, analyze deconvolved cibersort fractions.
#
# NOTE: Lambrechts normalized data is already normalized by seurat NormalizeData. So DO NOT re-normalize with seurat after creating seurat object.

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(Seurat)
library(NMF)
library(viridis)

setwd("~/Bioinformatics/sc_prognosis/")
source("scripts/sc_prognostic_scoring.R")
source("scripts/sc_precog_functions.R")
precog_data<-read.table("PRECOG/PRECOG_metaZ_table.txt",sep="\t",header=TRUE,row.names=1)

# preset info about dataset and phenotypes file
dataset<-"lambrechts_lung_adeno"
precog_label<-"Lung_cancer_ADENO"
cell_column<-"ID"
celltype_column<-"CellType"
sample_column<-"PatientNumber"

# read normalized counts file and phenotypes
counts<-as.data.frame(fread(paste0("standardized_data/",dataset,"/normalized_counts.tsv")))
counts<-column_to_rownames(counts, colnames(counts)[1])
phenotypes<-as.data.frame(fread(paste0("standardized_data/",dataset,"/phenotypes.tsv")))

fibroblast_cells<-phenotypes[,cell_column][phenotypes[,celltype_column]=="Fibroblasts"]

fibroblast_counts<-counts[,fibroblast_cells]

samples<-unique(phenotypes[,sample_column])

precog_z<-precog_data[,precog_label,drop=FALSE]
fibroblast_genes<-readRDS("outputs/LUAD_fibroblasts/lambrechts_fibroblast_marker_genes.rds")
precog_z<-precog_z[intersect(fibroblast_genes, rownames(precog_z)),,drop=FALSE]

# Classify fibroblasts using genes correlated with prognostic score and NMF.
bad_fibroblasts<-c()

for(s in samples){
  
  message(paste("analyzing:", s))
  
  # subset normalized data to 
  sample_cells<-phenotypes[,cell_column][phenotypes[,celltype_column]=="Fibroblasts"&phenotypes[,sample_column]==s]
  sample_counts<-fibroblast_counts[,sample_cells]
  
  # run prognostic scoring algorithm to score each cell
  scores<-run_prognostic_scoring(expression = sample_counts, prognostic_data = precog_z, cor_method = "pearson")
  
  # find top 100 positive and negative correlated genes with prognostic score
  gene_correlations<-cor(t(sample_counts), scores$prognostic_score, method = "spearman")
  gene_correlations<-gene_correlations[complete.cases(gene_correlations),,drop=FALSE]
  cluster_genes<-c(rownames(gene_correlations)[order(gene_correlations[,1])][1:100],
                   rownames(gene_correlations)[rev(order(gene_correlations[,1]))][1:100])
  
  # run nmf on data subset to correlated genes
  # using 5 runs and getting consensus, but can do more. Results shouldnt change much.
  nmf_input<-sample_counts[cluster_genes,]
  nmf_res <- nmf(nmf_input, 2, nrun=5, .options='v')
  nmf_clusters<-predict(nmf_res)
  
  # plot consensus heatmap
  message(paste("cophcor:",cophcor(nmf_res)))
  heatmap_input<-nmf_res@consensus
  h<-Heatmap(heatmap_input,
             col=viridis(10),
             column_title = paste("Correlated Genes Clustering\ncophcor:",round(cophcor(nmf_res),4)),
             show_row_names=FALSE,
             show_column_names=FALSE,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8),
             use_raster = TRUE, raster_quality = 1)
  ComplexHeatmap::draw(h)
  
  
  # use seurat to run pca
  seurat_data<- CreateSeuratObject(sample_counts, min.cells = 5, min.features = 0)
  seurat_data<-FindVariableFeatures(object=seurat_data, verbose = FALSE)
  seurat_data<-ScaleData(object=seurat_data, verbose = FALSE)
  seurat_data<-RunPCA(object = seurat_data, verbose=FALSE)
  pca<-seurat_data@reductions$pca@cell.embeddings[,1:5]
  
  metadata<-data.frame(cell=sample_cells,
                       nmf_cluster=as.factor(nmf_clusters[sample_cells]), 
                       score=scores$prognostic_score,
                       score_rank=rank(scores$prognostic_score))
  metadata<-cbind(metadata, pca)
  
  # plot pca colored by prognostic score
  plot<-ggplot(metadata, aes(x=PC_1,y=PC_2,color=score_rank)) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("Prognostic Score Rank") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5))
  print(plot)
  
  # plot violin plot of nmf clusters and score.
  plot1<-ggplot(metadata, aes(x=nmf_cluster,y=score)) +
    geom_violin() +
    geom_jitter(width=0.1) +
    theme_bw() +
    ggtitle("Clustered On Correlated Genes") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(method = "t.test", comparisons=list(c("1","2")))
  print(plot1)
  
  # plot nmf clusters on PCA plot
  plot5<-ggplot(metadata, aes(x=PC_1,y=PC_2,color=nmf_cluster)) +
    geom_point() +
    scale_color_manual(values=c("red","blue")) +
    ggtitle("Correlated Genes Clusters") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5))
  print(plot5)
  
  # determine if NMF clusters have significantly different prognistic scores.
  # If so, save cells in bad cluster.
  p<-t.test(metadata$score[metadata$nmf_cluster=="1"], metadata$score[metadata$nmf_cluster=="2"])$p.value
  diff<-mean(metadata$score[metadata$nmf_cluster=="1"])-mean(metadata$score[metadata$nmf_cluster=="2"])
  
  if(p<0.05){
    if(diff>0){
      bad_fibroblasts<-c(bad_fibroblasts, names(nmf_clusters)[nmf_clusters=="1"])
    }else{
      bad_fibroblasts<-c(bad_fibroblasts, names(nmf_clusters)[nmf_clusters=="2"])
    }
  }
}

write.table(bad_fibroblasts, "outputs/LUAD_fibroblasts/lambrechts_luad_bad_fibroblast_IDs_nmf.txt", row.names=FALSE, col.names = FALSE, quote=FALSE)



# differential expression
bad_fibroblasts<-fread("outputs/LUAD_fibroblasts/lambrechts_luad_bad_fibroblast_IDs_nmf.txt", header = FALSE)$V1
fibroblast_1<-bad_fibroblasts
fibroblast_2<-c(setdiff(fibroblast_cells[fibroblast_cells%in%phenotypes$ID[phenotypes$PatientNumber=="4"]],bad_fibroblasts),
                setdiff(fibroblast_cells[fibroblast_cells%in%phenotypes$ID[phenotypes$PatientNumber=="3"]],bad_fibroblasts))

# make seurat object, dont normalize, find markers.
f_seurat<-CreateSeuratObject(fibroblast_counts, min.cells = 5, min.features = 200)
f_seurat<-AddMetaData(object = f_seurat,
                      metadata = column_to_rownames(data.frame(ID=c(fibroblast_1,fibroblast_2),
                                                               label=c(rep("fib_1",length(fibroblast_1)),rep("fib_2",length(fibroblast_2)))),"ID"),
                      col.name = "type")
f_seurat<-SetIdent(f_seurat,value = "type")
diff_exp_seurat<-FindMarkers(f_seurat, ident.1 = "fib_1",ident.2 = "fib_2",logfc.threshold = 0)

colnames(diff_exp_seurat)<-c("p_val","avg_lfc","pct.1","pct.2","p_adj")

write.table(rownames_to_column(diff_exp_seurat,"Gene"),"outputs/LUAD_fibroblasts/lambrechts_luad_adv_vs_fav_fibroblasts_nmf_diff_exp.tsv",row.names = FALSE,sep="\t",quote=FALSE)


# heatmap of differential genes
diff_exp<-as.data.frame(fread("outputs/LUAD_fibroblasts/lambrechts_luad_adv_vs_fav_fibroblasts_nmf_diff_exp.tsv"))
diff_exp<-column_to_rownames(diff_exp, "Gene")

# scale each patient separately to gene z scores.
scaled_fibroblast_counts<-fibroblast_counts
for(s in samples){
  message(s)
  sample_cells<-intersect(colnames(scaled_fibroblast_counts),phenotypes$ID[phenotypes$PatientNumber==s])
  scaled_data<-as.data.frame(t(scale(t(scaled_fibroblast_counts[,sample_cells,drop=FALSE]))))
  scaled_fibroblast_counts<-cbind(scaled_fibroblast_counts[,!colnames(scaled_fibroblast_counts)%in%sample_cells],scaled_data)
}
scaled_fibroblast_counts[is.na(scaled_fibroblast_counts)]<-0

# get top genes by fold change.
fibroblast_1_genes<-rownames(diff_exp[rev(order(diff_exp$avg_lfc)),,drop=FALSE])[1:100]
fibroblast_2_genes<-rownames(diff_exp[order(diff_exp$avg_lfc),,drop=FALSE])[1:100]


prognostic_fibroblast_counts<-fibroblast_counts[c(fibroblast_1_genes,fibroblast_2_genes),c(fibroblast_1,fibroblast_2)]

# normalize matrix by iterative sum of squares. (code from Aaron Newman/Bogdan)
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

#sliding window smoothing
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

# make annotations
if(TRUE){
  rownames(phenotypes)<-NULL
  annotations<-column_to_rownames(phenotypes[,c("ID","CellType","PatientNumber")],"ID")
  colnames(annotations)<-c("CellType","sample")
  annotations<-annotations[annotations$CellType=="Fibroblasts",]
  annotations<-annotations[,"sample",drop=FALSE]
  annotations$CellType<-rep("Favorable Fibroblasts",nrow(annotations))
  annotations$CellType[rownames(annotations)%in%fibroblast_1]<-"Adverse Fibroblasts"
  annotations$sample<-as.factor(annotations$sample)
  annotations$CellType<-as.factor(annotations$CellType)
  annotations<-annotations[colnames(heatmap_input),,drop=FALSE]
  celltype_colors<-get_palette("ucscgb",k=length(levels(annotations$CellType)))
  names(celltype_colors)<-levels(annotations$CellType)
  sample_colors<-get_palette("Set1",k=length(levels(annotations$CellType)))
  names(sample_colors)<-levels(annotations$sample)
  anno_colors<-list(CellType=celltype_colors,sample=sample_colors) 
}

# plot with complexheatmap
h<-Heatmap(heatmap_input, 
           cluster_rows = FALSE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           show_column_names = FALSE,
           top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_legend_param = list(direction = "horizontal")),
           col = viridis(10),
           #col = colorRampPalette(c("lightblue","white","darkred"))(10),
           heatmap_legend_param = list(title = "Scaled Expression", nrow=1))
ComplexHeatmap::draw(h, annotation_legend_side = "bottom")





#### volcano plot of differential genes ####
library(ggrepel)

diff_genes<-as.data.frame(fread("outputs/LUAD_fibroblasts/lambrechts_luad_adv_vs_fav_fibroblasts_nmf_diff_exp.tsv"))

volcano_input<-diff_genes
volcano_input$color<-rep("black",nrow(volcano_input))
volcano_input$color[volcano_input$p_adj<0.05&(volcano_input$avg_lfc>0.5|volcano_input$avg_lfc<(-0.5))]<-"red"
interesting_genes<-c(volcano_input$Gene[order(volcano_input$avg_lfc)][1:10],volcano_input$Gene[rev(order(volcano_input$avg_lfc))][1:10])
volcano_input$annotation<-volcano_input$Gene
volcano_input$annotation[!volcano_input$annotation%in%interesting_genes]<-""
ggplot(volcano_input, aes(x=avg_lfc,y=-log10(p_adj), color=color, label=annotation)) +
  geom_point(size=1) +
  scale_color_manual(values=c("black","darkred")) +
  theme_bw() +
  xlab("Avg Log Fold Change") +
  ylab("-log10(Adjusted p-val)") +
  geom_label_repel(size=3) +
  theme(panel.border = element_rect(color="black", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(color="black"))



#### Fibroblast PCA plots ####
# PCA of patient, fibroblast classification, and prognostic scire labels
bad_fibroblasts<-fread("outputs/LUAD_fibroblasts/lambrechts_luad_bad_fibroblast_IDs_nmf.txt", header = FALSE)$V1
fibroblast_1<-bad_fibroblasts
fibroblast_2<-c(setdiff(fibroblast_cells[fibroblast_cells%in%phenotypes$ID[phenotypes$PatientNumber=="4"]],bad_fibroblasts),
                setdiff(fibroblast_cells[fibroblast_cells%in%phenotypes$ID[phenotypes$PatientNumber=="3"]],bad_fibroblasts))

f_seurat<-CreateSeuratObject(fibroblast_counts, min.cells = 5, min.features = 200)
f_seurat<-AddMetaData(object = f_seurat,
                      metadata = column_to_rownames(data.frame(ID=c(fibroblast_1,fibroblast_2),
                                                               label=c(rep("fib_1",length(fibroblast_1)),rep("fib_2",length(fibroblast_2)))),"ID"),
                      col.name = "type")
md<-phenotypes[match(rownames(f_seurat@meta.data), phenotypes$ID),c("ID","PatientNumber")]
md$PatientNumber<-as.character(md$PatientNumber)
rownames(md)<-NULL
f_seurat<-AddMetaData(object = f_seurat,
                      metadata = column_to_rownames(md,"ID"),
                      col.name = "patient")

f_seurat<-FindVariableFeatures(object=f_seurat)
f_seurat<-ScaleData(object=f_seurat)
f_seurat<-RunPCA(object = f_seurat, pc.genes = f_seurat@var.genes)

plot_metadata_pca(f_seurat, "patient",c("red","blue"))
plot_metadata_pca(f_seurat, "type",c("red","blue"))

fibroblast_scores<-scored_data[match(rownames(f_seurat@meta.data),scored_data$ID),]
fibroblast_scores$precog_score_pearson_rank<-rank(fibroblast_scores$precog_score_pearson)
rownames(fibroblast_scores)<-NULL
f_seurat<-AddMetaData(object = f_seurat, metadata = column_to_rownames(fibroblast_scores[,c("ID","precog_score_pearson_rank")],"ID"), col.name = "precog_score_rank")
f_seurat<-AddMetaData(object = f_seurat, metadata = column_to_rownames(fibroblast_scores[,c("ID","precog_score_pearson")],"ID"), col.name = "precog_score")
colors<-colorRampPalette(c("#000089","#0106FF","#00F5FF","#73FF86","#F7FF08","#FF1E00","#750000"))
p<-FeaturePlot(f_seurat, "precog_score_rank",reduction = "pca") +
  scale_color_viridis() +
  #scale_color_gradientn(colours = colors(20)) +
  ggtitle("Lambrechts All Fibroblasts")
print(p)



#### Generate Cibersort Input ####
bad_fibroblasts<-fread("outputs/LUAD_fibroblasts/lambrechts_luad_bad_fibroblast_IDs_nmf.txt", header = FALSE)$V1
fibroblast_1<-bad_fibroblasts
fibroblast_2<-c(setdiff(fibroblast_cells[fibroblast_cells%in%phenotypes$ID[phenotypes$PatientNumber=="4"]],bad_fibroblasts),
                setdiff(fibroblast_cells[fibroblast_cells%in%phenotypes$ID[phenotypes$PatientNumber=="3"]],bad_fibroblasts))

celltype_reference<-phenotypes[,c("ID","CellType"),drop=FALSE]
celltype_reference$CellType[celltype_reference$CellType%in%c("CD8 T cells", "Monocytes and Macrophages", "CD4 T cells", "B cells", "Tregs", "NK cells", "Granulocytes", "Dendritic cells")]<-"Immune"
celltype_reference$CellType[celltype_reference$CellType%in%c("PCs", "", "Normal epithelia", "Erythroblasts", "Mast cells")]<-"Other"
celltype_reference$CellType[celltype_reference$CellType=="Cancer epithelia"]<-"Tumor"
celltype_reference$CellType[celltype_reference$CellType=="Endothelial cells"]<-"Endothelial"
celltype_reference<-celltype_reference[which(celltype_reference$CellType!="Other"),]

celltype_reference$CellType[celltype_reference$ID%in%fibroblast_1]<-"Fibroblast_adv"
celltype_reference$CellType[celltype_reference$ID%in%fibroblast_2]<-"Fibroblast_fav"

rownames(celltype_reference)<-NULL
celltype_reference<-column_to_rownames(celltype_reference,"ID")

coi<-sample(rownames(celltype_reference),5000)
celltype_reference<-celltype_reference[coi,,drop=FALSE]

cibersort_input<-counts
cibersort_input<-cibersort_input[,coi]
cibersort_input<-(2^cibersort_input)-1
cibersort_input<-sweep(cibersort_input,2,colSums(cibersort_input),'/')
cibersort_input<-cibersort_input*1000000

cibersort_input<-rbind(celltype_reference$CellType,cibersort_input)

write.table(cibersort_input,"cibersort_inputs/luad_rename_fibroblast_split_csx_input.tsv",sep="\t",col.names = FALSE, row.names = TRUE, quote = FALSE)

cibersort_input[cibersort_input=="Fibroblast_adv"]<-"Fibroblasts"
cibersort_input[cibersort_input=="Fibroblast_fav"]<-"Fibroblasts"

write.table(cibersort_input,"cibersort_inputs/luad_rename_fibroblast_whole_csx_input.tsv",sep="\t",col.names = FALSE, row.names = TRUE, quote = FALSE)

