# Armon Azizi
#
# Find lambrechts fibroblast markers

# find fibroblast specific genes in lambrechts


library(data.table)
library(Seurat)

setwd("~/Bioinformatics/sc_prognosis/")

counts<-as.data.frame(fread("standardized_data/lambrechts_lung_adeno/normalized_counts.tsv"))
counts<-column_to_rownames(counts, colnames(counts)[1])
phenotypes<-as.data.frame(fread("standardized_data/lambrechts_lung_adeno/phenotypes.tsv"))

seurat_obj<-CreateSeuratObject(counts, min.cells = 5, min.features = 200)

celltype_ref<-phenotypes[match(colnames(seurat_obj),phenotypes$ID),c("ID","CellType")]
rownames(celltype_ref)<-NULL
seurat_obj<-AddMetaData(object = seurat_obj,
                        metadata = column_to_rownames(celltype_ref, "ID"),
                        col.name = "celltype")
seurat_obj<-SetIdent(seurat_obj,value = "celltype")
diff_exp_seurat<-FindMarkers(seurat_obj, ident.1 = "Fibroblasts")
fibroblast_markers<-rownames(diff_exp_seurat)[diff_exp_seurat$p_val_adj<0.05&diff_exp_seurat$avg_logFC>0]

saveRDS(fibroblast_markers, "outputs/LUAD_fibroblasts/lambrechts_fibroblast_marker_genes.rds")
