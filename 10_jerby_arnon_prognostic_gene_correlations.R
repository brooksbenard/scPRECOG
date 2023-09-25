# Armon Azizi
#
# jerby arnon melanoma malignant cell analysis
# Find genes that are correlated with prognostic score in samples that have differentiation heterogeneity
# Generate correlation heatmap
# generate correlation comparisons to other cell types
# generate gsea on correlated genes.

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(NMF)
library(circlize)

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

celltypes<-unique(phenotypes$CellType)

celltype_prognostic_corr<-data.frame(matrix(nrow=nrow(counts), ncol=0))
rownames(celltype_prognostic_corr)<-rownames(counts)


# Malignant correlations
for(s in samples){
  if(length(phenotypes$ID[phenotypes$CellType=="Mal"&phenotypes$sample==s])<10){next}
  
  sample_celltype_counts<-counts[,phenotypes$ID[phenotypes$sample==s&phenotypes$CellType=="Mal"]]
  
  scores<-run_prognostic_scoring(expression = sample_celltype_counts, prognostic_data = precog_z_malignant, cor_method = "pearson")
  
  gene_correlations<-cor(t(log2(sample_celltype_counts+1)), scores$prognostic_score, method = "spearman")
  
  gene_correlations<-as.data.frame(gene_correlations)
  colnames(gene_correlations)<-c(paste("Mal",s,sep="_"))
  
  celltype_prognostic_corr<-cbind(celltype_prognostic_corr, gene_correlations)
}


malignant_corr<-celltype_prognostic_corr[,grep("Mal",colnames(celltype_prognostic_corr))]
malignant_corr[is.na(malignant_corr)]<-0
saveRDS(malignant_corr, "outputs/melanoma/melanoma_malignant_gene_correlations_to_precog_score.rds")
malignant_corr$mean_correlation<-rowMeans(malignant_corr)
write.table(rownames_to_column(malignant_corr,"Gene"),"outputs/melanoma/melanoma_malignant_gene_correlations_to_precog_score.tsv",row.names = F,quote=F,sep="\t")
malignant_corr<-malignant_corr[,colnames(malignant_corr)!="mean_correlation"]
top_genes<-c(rownames(malignant_corr[order(rowSums(malignant_corr)),])[1:100],
             rownames(malignant_corr[rev(order(rowSums(malignant_corr))),])[1:100])
plot_input<-malignant_corr[top_genes,]

sample_order<-sample_info[paste("Mal",sample_info$Sample,sep="_")%in%colnames(plot_input),]
sample_order<-sample_order$Sample[order(sample_order$lymph_node_status)]
plot_input<-plot_input[,paste("Mal",sample_order,sep="_")]

plot_genes<-rownames(plot_input)[c(1:20,101:125)]
genes_idx<-match(plot_genes, rownames(plot_input))

if(TRUE){
  annotations<-sample_info[,c("Sample","lymph_node_status")]
  rownames(annotations)<-NULL
  annotations<-column_to_rownames(annotations,"Sample")
  colnames(annotations)<-c("Metastasis")
  annotations$Metastasis<-as.factor(annotations$Metastasis)
  annotations<-annotations[sapply(strsplit(colnames(plot_input),"_"),'[',2),,drop=FALSE]
  sample_colors<-get_palette("ucscgb",k=length(levels(annotations$Metastasis)))
  names(sample_colors)<-levels(annotations$Metastasis)
  anno_colors<-list(Metastasis=sample_colors) 
}

Heatmap(plot_input, 
        cluster_rows = FALSE,
        cluster_columns=FALSE,
        show_row_names = FALSE,
        show_column_names = F,
        #col = viridis(10),
        top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors),
        left_annotation = rowAnnotation(genes = anno_mark(at = genes_idx, 
                                                          labels = plot_genes, 
                                                          side="left",
                                                          labels_gp = gpar(fontsize = 8),
                                                          padding = 2,link_width = unit(15, "mm"))),
        col = colorRamp2(seq(-0.1,0.1,length.out = 10), viridis(10)),
        heatmap_legend_param = list(title = "Corelation\nWith Precog Score", nrow=1))


# only stem samples
stem_samples<-readRDS("outputs/melanoma/jerby_arnon_samples_with_stem_heterogeneity.rds")

malignant_corr<-readRDS("outputs/melanoma/melanoma_malignant_gene_correlations_to_precog_score.rds")
malignant_corr<-malignant_corr[,paste("Mal",stem_samples,sep="_")]
top_genes<-c(rownames(malignant_corr[order(rowSums(malignant_corr)),])[1:100],
             rownames(malignant_corr[rev(order(rowSums(malignant_corr))),])[1:100])
plot_input<-malignant_corr[top_genes,]
#plot_input<-scale(plot_input)

sample_order<-sample_info[paste("Mal",sample_info$Sample,sep="_")%in%colnames(plot_input),]
sample_order<-sample_order$Sample[order(sample_order$lymph_node_status)]
plot_input<-plot_input[,paste("Mal",sample_order,sep="_")]

plot_genes<-rownames(plot_input)[c(1:20,101:125)]
genes_idx<-match(plot_genes, rownames(plot_input))

if(TRUE){
  annotations<-sample_info[,c("Sample","lymph_node_status")]
  rownames(annotations)<-NULL
  annotations<-column_to_rownames(annotations,"Sample")
  colnames(annotations)<-c("Metastasis")
  annotations$Metastasis<-as.factor(annotations$Metastasis)
  annotations<-annotations[sapply(strsplit(colnames(plot_input),"_"),'[',2),,drop=FALSE]
  sample_colors<-get_palette("ucscgb",k=length(levels(annotations$Metastasis)))
  names(sample_colors)<-levels(annotations$Metastasis)
  anno_colors<-list(Metastasis=sample_colors) 
}

Heatmap(plot_input, 
        cluster_rows = FALSE,
        cluster_columns=FALSE,
        show_row_names = FALSE,
        show_column_names = F,
        #col = viridis(10),
        top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors),
        left_annotation = rowAnnotation(genes = anno_mark(at = genes_idx, 
                                                          labels = plot_genes, 
                                                          side="left",
                                                          labels_gp = gpar(fontsize = 8),
                                                          padding = 2,link_width = unit(15, "mm"))),
        col = colorRamp2(seq(-0.1,0.1,length.out = 10), viridis(10)),
        heatmap_legend_param = list(title = "Corelation\nWith Precog Score", nrow=1))


# compare celltype correlations

celltype_prognostic_corr<-data.frame(matrix(nrow=nrow(counts), ncol=0))
rownames(celltype_prognostic_corr)<-rownames(counts)

for(c in c("Mal","Macrophage","Endo.","T.CD4","CAF","T.CD8","NK","B.cell")){
  
  for(s in samples){
    if(length(phenotypes$ID[phenotypes$CellType==c&phenotypes$sample==s])<10){next}
    
    sample_celltype_counts<-counts[,phenotypes$ID[phenotypes$sample==s&phenotypes$CellType==c]]
    
    celltype_markers<-readRDS("outputs/melanoma/jerby_arnon_celltype_marker_genes.rds")
    celltype_genes<-rownames(celltype_markers)[celltype_markers$cluster==c&celltype_markers$p_val_adj<0.05&celltype_markers$avg_logFC>0]
    precog_z_celltype<-precog_z[intersect(celltype_genes, rownames(precog_z)),,drop=FALSE]
    
    scores<-run_prognostic_scoring(expression = sample_celltype_counts, prognostic_data = precog_z_malignant, cor_method = "pearson")
    
    gene_correlations<-cor(t(log2(sample_celltype_counts+1)), scores$prognostic_score, method = "spearman")
    
    gene_correlations<-as.data.frame(gene_correlations)
    colnames(gene_correlations)<-c(paste(c,s,sep="_"))
    
    celltype_prognostic_corr<-cbind(celltype_prognostic_corr, gene_correlations)
  }
}


top_genes<-c()
for(c in c("Mal","Macrophage","Endo.","T.CD4","CAF","T.CD8","NK","B.cell")){
  celltype_corr<-celltype_prognostic_corr[,grep(c,colnames(celltype_prognostic_corr)),drop=FALSE]
  celltype_corr[is.na(celltype_corr)]<-0
  top_genes<-c(top_genes,
               rownames(celltype_corr[order(rowSums(celltype_corr)),,drop=FALSE])[1:100],
               rownames(celltype_corr[rev(order(rowSums(celltype_corr))),,drop=FALSE])[1:100])
}

plot_input<-celltype_prognostic_corr[top_genes,]
plot_input[is.na(plot_input)]<-0


Heatmap(plot_input, 
        cluster_rows = FALSE,
        cluster_columns=FALSE,
        show_row_names = FALSE,
        show_column_names = T,
        column_names_gp = gpar(fontsize = 5),
        #col = viridis(10),
        col = colorRamp2(seq(-0.1,0.1,length.out = 10), viridis(10)),
        heatmap_legend_param = list(title = "Scaled Corelation\nWith Precog Score", nrow=1))


# top 100 positively correlated genes

malignant_corr<-celltype_prognostic_corr[,grep("Mal",colnames(celltype_prognostic_corr))]
malignant_corr[is.na(malignant_corr)]<-0
top_genes<-c(rownames(malignant_corr[rev(order(rowSums(malignant_corr))),])[1:100])

celltype_prognostic_corr_zeroed<-celltype_prognostic_corr
celltype_prognostic_corr_zeroed[is.na(celltype_prognostic_corr_zeroed)]<-0

plot_input<-data.frame(matrix(nrow=length(top_genes), ncol=0))

temp<-data.frame(a=rowMeans(celltype_prognostic_corr_zeroed[top_genes,grep("Mal",colnames(celltype_prognostic_corr_zeroed)),drop=FALSE]))
colnames(temp)<-c("Malignant")
plot_input<-cbind(plot_input, temp)
temp<-data.frame(a=rowMeans(celltype_prognostic_corr_zeroed[top_genes,grep("Mal",colnames(celltype_prognostic_corr_zeroed), invert = T),drop=FALSE]))
colnames(temp)<-c("Non-Malignant")
plot_input<-cbind(plot_input, temp)

plot_input<-reshape2::melt(as.matrix(plot_input))

ggplot(plot_input, aes(x=Var2, y=value, fill=Var2)) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.size=-1) +
  scale_x_discrete(labels=c("Malignant","Non-Malignant\nCelltypes Average")) +
  scale_fill_manual(values = get_palette("ucscgb", length(unique(plot_input$Var2)))) +
  xlab(NULL)+
  ylab("Average correlation of gene\nwith score across all patients") +
  theme_bw() +
  #scale_x_discrete(limits=order) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=10),
        legend.position = "none") +
  stat_compare_means(method="t.test", paired = F, comparisons = list(c("Malignant","Non-Malignant")))


# top 100 negatively correlated genes

malignant_corr<-celltype_prognostic_corr[,grep("Mal",colnames(celltype_prognostic_corr))]
malignant_corr[is.na(malignant_corr)]<-0
top_genes<-c(rownames(malignant_corr[order(rowSums(malignant_corr)),])[1:100])

celltype_prognostic_corr_zeroed<-celltype_prognostic_corr
celltype_prognostic_corr_zeroed[is.na(celltype_prognostic_corr_zeroed)]<-0

plot_input<-data.frame(matrix(nrow=length(top_genes), ncol=0))

temp<-data.frame(a=rowMeans(celltype_prognostic_corr_zeroed[top_genes,grep("Mal",colnames(celltype_prognostic_corr_zeroed)),drop=FALSE]))
colnames(temp)<-c("Malignant")
plot_input<-cbind(plot_input, temp)
temp<-data.frame(a=rowMeans(celltype_prognostic_corr_zeroed[top_genes,grep("Mal",colnames(celltype_prognostic_corr_zeroed), invert = T),drop=FALSE]))
colnames(temp)<-c("Non-Malignant")
plot_input<-cbind(plot_input, temp)

plot_input<-reshape2::melt(as.matrix(plot_input))

ggplot(plot_input, aes(x=Var2, y=value, fill=Var2)) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.size=-1) +
  scale_x_discrete(labels=c("Malignant","Non-Malignant\nCelltypes Average")) +
  scale_fill_manual(values = get_palette("ucscgb", length(unique(plot_input$Var2)))) +
  xlab(NULL)+
  ylab("Average correlation of gene\nwith score across all patients") +
  theme_bw() +
  #scale_x_discrete(limits=order) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=10),
        legend.position = "none") +
  stat_compare_means(method="t.test", paired = F, comparisons = list(c("Malignant","Non-Malignant")))




# gsea on correlated genes
library(fgsea)
library(qusage)
stem_samples<-readRDS("outputs/melanoma/jerby_arnon_samples_with_stem_heterogeneity.rds")

malignant_corr<-readRDS("outputs/melanoma/melanoma_malignant_gene_correlations_to_precog_score.rds")
malignant_corr<-malignant_corr[,paste("Mal",stem_samples,sep="_")]

ranked_genes<-rowMeans(malignant_corr)
ranked_genes<-ranked_genes[ranked_genes!=0]

pathways1<-read.gmt("gene_sets/msigdb_curated_c2.gmt")
pathways2<-read.gmt("gene_sets/msigdb_go_bioprocess.gmt")
pathways<-c(pathways1,pathways2)

fgsea_res<-fgsea(pathways, ranked_genes)
write.xlsx(fgsea_res, "outputs/melanoma/prognostic_gene_correlations_gsea.xlsx")
