# Armon Azizi
#
# heatmap of favorable and unfavorable genes for each celltype

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(fgsea)
library(qusage)

setwd("~/Bioinformatics/sc_prognosis/")

source("scripts/sc_precog_functions.R")

precog_data<-read.table("PRECOG/PRECOG_metaZ_table.txt",sep="\t",header=TRUE,row.names=1)

celltype_order<-c("Fibroblasts","Endothelial","Epithelial/Malignant","Other","Macrophage w/ (Mono/Myeloid)","Other Immune","T Cells","Dendritic","NK Cells","B Cells","Mast Cells")

celltype_mapping<-c("Fibroblasts"="Fibroblasts",
                    "Endothelial"="Endothelial",
                    "Epithelial/Malignant"="Epithelial.Malignant",
                    "Other"="Other",
                    "Macrophage w/ (Mono/Myeloid)"="mac.mono.myeloid",
                    "Other Immune"="other.immune",
                    "T Cells"="t.cells",
                    "Dendritic"="dendritic",
                    "NK Cells"="NK",
                    "B Cells"="b.cells",
                    "Mast Cells"="mast.cells")

dataset_info<-read.xlsx("dataset_info/dataset_info.xlsx")

scored_data<-read.xlsx("outputs/all_datasets_prognostic_scores.xlsx")

# Load data
all_counts<-readRDS("combined_datasets/normalized_counts.rds")
all_phenotypes<-readRDS("combined_datasets/phenotypes.rds")



#### Do differential expression between celltypes for all datasets ####
for(i in 1:nrow(dataset_info)){
  
  dataset<-dataset_info[i,"dataset"]
  precog_label<-dataset_info[i,"precog_label"]
  tcga_label<-dataset_info[i,"tcga_label"]
  counts_file<-dataset_info[i,"counts_file"]
  cell_column<-dataset_info[i,"cell_column"]
  celltype_column<-dataset_info[i,"celltype_column"]
  sample_column<-dataset_info[i,"sample_column"]
  
  message(dataset)
  
  cells_subset<-all_phenotypes$ID[all_phenotypes$dataset==dataset]
  
  counts<-all_counts[,cells_subset]
  phenotypes<-all_phenotypes[match(cells_subset,all_phenotypes$ID),]
  
  phenotypes<-rename_celltypes(phenotypes, "CellType")
  
  # generate seurat object
  seurat_data<-CreateSeuratObject(counts, min.cells = 5, min.features = 200)
  if(!dataset%in%c("lambrechts_lung_adeno","lambrechts_lung_scc")){
    seurat_data<-NormalizeData(seurat_data)
  }
  seurat_data <- ScaleData(object=seurat_data, features=rownames(seurat_data))
  
  # label celltypes
  celltype_labels<-phenotypes[match(rownames(seurat_data@meta.data),phenotypes$ID),c("ID","CellType")]
  celltype_labels[,"CellType"]<-as.character(celltype_labels[,"CellType"])
  rownames(celltype_labels)<-NULL
  celltype_labels<-column_to_rownames(celltype_labels,"ID")
  seurat_data<-AddMetaData(object = seurat_data, metadata = celltype_labels, col.name = "celltype")
  seurat_data<-SetIdent(seurat_data,value="celltype")
  
  # find differential genes
  markers<-FindAllMarkers(seurat_data, max.cells.per.ident = 500)
  
  write.table(rownames_to_column(markers,"Gene"),paste0("outputs/celltype_differential_expression/",dataset,"_findmarkers_output.tsv"),sep="\t", row.names = FALSE)
}


#### Combine datasets and find significant celltype markers ####

set.seed(1)
cells_subset<-sample(all_phenotypes$ID,5000,replace = FALSE)
counts<-all_counts[,cells_subset]
phenotypes<-all_phenotypes[match(cells_subset,all_phenotypes$ID),]

phenotypes<-rename_celltypes(phenotypes, "CellType")


favorable_genes<-rownames(precog_data)[precog_data$Unweighted_meta.Z_of_all_cancers<(-2)]
adverse_genes<-rownames(precog_data)[precog_data$Unweighted_meta.Z_of_all_cancers>2]

prognostic_counts<-counts[intersect(rownames(counts),c(adverse_genes,favorable_genes)),]

# normalize each dataset by gene z-scores
for(d in unique(phenotypes$dataset)){
  message(d)
  dataset_cells<-intersect(colnames(prognostic_counts),phenotypes$ID[phenotypes$dataset==d])
  scaled_data<-as.data.frame(t(scale(t(prognostic_counts[,dataset_cells,drop=FALSE]))))
  prognostic_counts<-cbind(prognostic_counts[,!colnames(prognostic_counts)%in%dataset_cells],scaled_data)
}

prognostic_counts[is.na(prognostic_counts)]<-0

rownames(phenotypes)<-NULL
cell_labels<-column_to_rownames(phenotypes[,c("ID","CellType","dataset")],"ID")
cell_labels$group<-as.factor(paste(phenotypes$CellType,phenotypes$dataset,sep = "_"))
cell_labels$CellType<-as.factor(phenotypes$CellType)

# find genes that are associated with each celltype
gene_assignments<-reorder_dataframe(prognostic_counts, 
                                    cell_labels[,"CellType",drop=FALSE],
                                    sample_order=sort(unique(cell_labels$CellType)), 
                                    return_assignments = TRUE)

gene_assignments_1<-gene_assignments
gene_assignments_1$meta_precog_prognosticity<-rep("Favorable",nrow(gene_assignments_1))
gene_assignments_1$meta_precog_prognosticity[gene_assignments_1$gene%in%adverse_genes]<-"Adverse"
write.table(gene_assignments_1, "outputs/celltype_prognostic_gene_assignments.tsv",sep="\t",row.names=FALSE,quote=FALSE)

gene_assignments$score<-rep(0,nrow(gene_assignments))

for(d in unique(phenotypes$dataset)){
  message(d)
  differential_data<-as.data.frame(fread(paste0("outputs/celltype_differential_expression/",d,"_findmarkers_output.tsv")))
  differential_data<-differential_data[differential_data$p_val_adj<0.05,]
  for(i in 1:nrow(gene_assignments)){
    if(sum(differential_data$cluster==gene_assignments[i,"CellType"]&differential_data$gene==gene_assignments[i,"gene"])>0){
      gene_assignments[i,"score"]<-gene_assignments[i,"score"]+1
    }
  }
}


final_genes<-gene_assignments$gene[gene_assignments$score>1]

prognostic_counts<-prognostic_counts[final_genes,]

prognostic_counts<-reorder_dataframe(prognostic_counts, cell_labels[,"CellType",drop=FALSE],sample_order=sort(unique(cell_labels$CellType)), priority_features = intersect(adverse_genes,rownames(prognostic_counts)))

new_order<-sort(unique(cell_labels$group))
cell_order<-rownames(cell_labels)[order(match(cell_labels$group,new_order))]
prognostic_counts<-prognostic_counts[,cell_order]

# heatmap_input<-prognostic_counts
# heatmap_input<-apply(heatmap_input, 1, function(x) x - mean(x))
# heatmap_input<-apply(heatmap_input, 1, function(x) x/sd(x))


hmcs<-prognostic_counts
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
  w <- 2
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

heatmap_input<-heatmap_input-min(heatmap_input)
heatmap_input[intersect(rownames(heatmap_input),favorable_genes),]<-heatmap_input[intersect(rownames(heatmap_input),favorable_genes),]*-1

if(TRUE){
  rownames(phenotypes)<-NULL
  annotations<-column_to_rownames(phenotypes[,c("ID","CellType","dataset")],"ID")
  colnames(annotations)<-c("CellType","dataset")
  annotations$CellType<-as.factor(annotations$CellType)
  annotations$dataset<-as.factor(annotations$dataset)
  annotations<-annotations[colnames(heatmap_input),,drop=FALSE]
  celltype_colors<-get_palette("ucscgb",k=length(levels(annotations$CellType)))
  names(celltype_colors)<-levels(annotations$CellType)
  dataset_colors<-get_palette("Set1",k=length(levels(annotations$dataset)))
  names(dataset_colors)<-levels(annotations$dataset)
  anno_colors<-list(CellType=celltype_colors,dataset=dataset_colors) 
}

Heatmap(heatmap_input, 
        cluster_rows = FALSE,
        cluster_columns=FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors),
        col = colorRampPalette(c("#FDE725FF","#21908CFF","#440154FF","#21908CFF","#FDE725FF"))(10),
        #col = colorRampPalette(c("green","black","black","red"))(10),
        heatmap_legend_param = list(title = "Expression"))

Heatmap(heatmap_input, 
        cluster_rows = FALSE,
        cluster_columns=FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors),
        #col = colorRampPalette(c("#FDE725FF","#21908CFF","#440154FF","#21908CFF","#FDE725FF"))(10),
        col = colorRampPalette(c("green","black","black","red"))(10),
        heatmap_legend_param = list(title = "Expression"), use_raster = TRUE, raster_quality = 10, raster_device = "png")


## gene set enrichment
dir.create("outputs/pathway_enrichment/meta_enrichment", recursive = TRUE)

gene_assignments<-as.data.frame(fread("outputs/celltype_prognostic_gene_assignments.tsv"))[,c(1,2)]

gene_assignments$prognosis<-rep(NA,nrow(gene_assignments))
gene_assignments$prognosis[gene_assignments$gene%in%favorable_genes]<-"favorable"
gene_assignments$prognosis[gene_assignments$gene%in%adverse_genes]<-"adverse"

write.xlsx(gene_assignments, paste0("outputs/pathway_enrichment/meta_enrichment/gene_assignments.xlsx"))

pathways1<-read.gmt("gene_sets/msigdb_curated_c2.gmt")
pathways2<-read.gmt("gene_sets/msigdb_go_bioprocess.gmt")

pathways<-c(pathways1,pathways2)

all_genes<-rownames(all_counts)

for(c in unique(gene_assignments$CellType)){
  for(p in c("adverse","favorable")){
    
    message(paste(c,p))
    
    goi<-gene_assignments$gene[gene_assignments$CellType==c&gene_assignments$prognosis==p]
    result<-data.frame()
    
    for(i in 1:length(pathways)){
      pathway<-pathways[[i]]
      
      pathway_genes<-intersect(all_genes,pathway)
      
      p_val<-dhyper(length(intersect(goi,pathway_genes)),length(pathway_genes),length(all_genes)-length(pathway_genes),length(goi))
      
      lfc<-log2(length(intersect(goi,pathway_genes))/(length(goi)*(length(pathway_genes)/length(all_genes))))
      
      result[i,"pathway"]<-names(pathways)[i]
      result[i,"lfc"]<-lfc
      result[i,"p"]<-p_val
      result[i,"num_pathway_genes"]<-length(pathway_genes)
      result[i,"num_celltype_genes"]<-length(goi)
      result[i,"num_intersect"]<-length(intersect(goi,pathway_genes))
      result[i,"intersect_genes"]<-paste(intersect(goi,pathway_genes), collapse = ";")
      
      result<-result[result$lfc>(-Inf),]
      result<-result[result$p<0.05,]
      result<-result[complete.cases(result),]
    }
    
    write.xlsx(result, paste0("pathway_enrichment/meta_enrichment/",celltype_mapping[c],"_",p,"_genes_pathway_enrichment.xlsx"))
  }
}


# plot curated pathways
pathway_data<-read.xlsx("outputs/pathway_enrichment/curated_pathways.xlsx")

ggplot(pathway_data, aes(x=group, y=-log10(p_val), label=pathway)) +
  geom_bar(stat="identity", fill="darkblue", color="black", width = 0.15) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_text(aes(label=pathway), vjust=-1.5, hjust=1, size=3) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))



# plot number of genes assigned per celltype
gene_assignments<-read.xlsx("outputs/pathway_enrichment/meta_enrichment/gene_assignments.xlsx")
gene_assignment_counts<-as.data.frame(table(gene_assignments[,c(2,3)]))

ggplot(gene_assignment_counts, aes(x=CellType,y=Freq,group=prognosis, fill=prognosis)) +
  geom_bar(position="dodge", stat="identity", color="black") +
  scale_fill_manual(values=c("red","green")) +
  ylab("Number of Genes") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=12, angle = 45, hjust=1),
        axis.text.y = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        legend.text = element_text(color="black", size=12),
        legend.title = element_text(color="black", size=14))
