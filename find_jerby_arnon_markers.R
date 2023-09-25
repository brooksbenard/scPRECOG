# find celltype specific genes in jerby-arnon

setwd("~/Bioinformatics/sc_prognosis/")

source("scripts/sc_prognostic_scoring.R")
source("scripts/sc_precog_functions.R")

precog_data<-read.table("/Users/majetilab/Bioinformatics/PRECOG/PRECOG_metaZ_table.txt",sep="\t",header=TRUE,row.names=1)

precog_z<-precog_data[,"Melanoma_Metastasis",drop=FALSE]

counts<-as.data.frame(fread("standardized_data/jerby_arnon_melanoma/normalized_counts.tsv"))
counts<-column_to_rownames(counts, colnames(counts)[1])
phenotypes<-as.data.frame(fread("standardized_data/jerby_arnon_melanoma/phenotypes.tsv"))

seurat_obj<-CreateSeuratObject(counts, min.cells = 5, min.features = 200)
seurat_obj<-NormalizeData(seurat_obj)

celltype_ref<-phenotypes[match(colnames(seurat_obj),phenotypes$ID),c("ID","CellType")]
rownames(celltype_ref)<-NULL
seurat_obj<-AddMetaData(object = seurat_obj,
                        metadata = column_to_rownames(celltype_ref, "ID"),
                        col.name = "celltype")
seurat_obj<-SetIdent(seurat_obj,value = "celltype")
diff_exp_seurat<-FindAllMarkers(seurat_obj)

saveRDS(diff_exp_seurat, "outputs/melanoma/jerby_arnon_celltype_marker_genes.rds")
