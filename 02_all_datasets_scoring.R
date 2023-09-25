# Armon Azizi
#
# Do prognostic scoring on all standardized datasets

setwd("~/Bioinformatics/sc_prognosis/")

library(data.table)
library(tibble)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)

source("scripts/sc_prognostic_scoring.R")

precog_data<-read.table("PRECOG/PRECOG_metaZ_table.txt",sep="\t",header=TRUE,row.names=1)

dataset_info<-read.xlsx("dataset_info/dataset_info.xlsx")

final_data<-data.frame()

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
  
  precog_z<-precog_data[,precog_label,drop=FALSE]
  
  precog_scores_pearson<-run_prognostic_scoring(expression = counts, prognostic_data = precog_z, cor_method = "pearson")
  
  if(!is.na(sample_column)){
    sample_info<-phenotypes[match(rownames(precog_scores_pearson),phenotypes[,cell_column]),sample_column]
  }else{
    sample_info<-rep(NA, nrow(precog_scores_pearson))
  }
  
  result<-data.frame(ID=rownames(precog_scores_pearson),
                     dataset=rep(dataset,nrow(precog_scores_pearson)),
                     sample_info=sample_info,
                     celltype=phenotypes[match(rownames(precog_scores_pearson),phenotypes[,cell_column]),celltype_column],
                     precog_score_pearson=precog_scores_pearson$prognostic_score)
  
  final_data<-rbind(final_data,result)
}

write.xlsx(final_data, "outputs/all_datasets_prognostic_scores.xlsx")


# number of cells in each celltype
scored_data<-read.xlsx("outputs/all_datasets_prognostic_scores.xlsx")
res<-as.data.frame(table(paste(scored_data$dataset,scored_data$celltype,sep=",")))
res<-data.frame(dataset=sapply(strsplit(as.character(res$Var1),","),'[',1),
                celltype=sapply(strsplit(as.character(res$Var1),","),'[',2),
                count=res$Freq)
write.xlsx(res,"outputs/celltype_counts.xlsx")


## Heatmap of datasets vs celltypes
scored_data<-read.xlsx("outputs/all_datasets_prognostic_scores.xlsx")

tcells<-c("CD4 T cells","CD8 T cells","T-Cell","Tregs","T Cells","Tcell","T cell","T and NK cells","T cells","NKT cells","gd T cells",
          "T:CD4+NAIVE","T:CD8+NAIVE","T:CD8+EM","T:CD4+CM","T:Reg","T:CD8+CM","T:CD4+EM","CD4 T-cells","Tm","Treg","Th","CD8 T-cells")
bcells<-c("B cells","B-Cell","Bcell","B cell","B:","B-cells","Breg")
mono_mac<-c("Monocytes and Macrophages","Macrophage","TAM","Myeloid & Macrophages","MACROPHAGE:","MONOCYTE:precursor","MONOCYTE:","MACROPHAGE","MONOCYTE","MICROGLIA/MACROPHAGE")
mast<-c("Mast cells","Mast","MastCell","MAST:","MAST")
nk<-c("NK cells","NK","NK:CD56+16+3+NKT","NK:CD56-16+3-","NK:CD56+16+3-","NKT")
other_immune<-c("PCs","Granulocytes","Immune","pDCs","NEUTROPHIL:","NEUTROPHIL","PMNs")
fibroblasts<-c("Fibroblasts","Cancer-Associated Fibroblast","CAF","Fibroblast","normal fibroblasts","metastatic fibroblasts","primary fibroblasts","FIBROBLAST")
dendritic<-c("Dendritic cells","Dendritic","Activated DC","mDC:","pDC:","DENDRITIC","DENDRITIC (ACTIVATED)")
epithelial_malignant<-c("Cancer epithelia","Normal epithelia","Malignant","Tumor","Malignant cell","Epithelial","Ductal cells 1","Ductal cells 2",
                        "Ductal cells 3","Malignant cells","HGSOC epithelial","HGSOC-F epithelial","HG3 epithelial","EPITHELIAL","LGSOC epithelial", "Epithelial cells")
endothelial<-c("Endothelial","TEC","Endothelial cells","endothelial progenitor","ENDOTHELIAL")
other<-c("","Erythroblasts","Unknown","Stromal","Myeloid","unclassified","HPC-like","Plasma cells",
         "Fenestrated EC","RBCs","Perivascular cells","Acinar cells","myocyte","Myeloid cells","metastatic myeloid lineage",
         "primary myeloid lineage","cancer stromal cells","PERICYTE","MDSC","PROLIFERATING MESENCHYMAL PROGENITOR","IG","normal stromal cells","benign mesothelial")

celltype_assignments<-data.frame(tcells=c(tcells,rep(NA,length(tcells)-length(tcells))),
                                 bcells=c(bcells,rep(NA,length(tcells)-length(bcells))),
                                 mono_mac=c(mono_mac,rep(NA,length(tcells)-length(mono_mac))),
                                 mast=c(mast,rep(NA,length(tcells)-length(mast))),
                                 other_immune=c(other_immune,rep(NA,length(tcells)-length(other_immune))),
                                 fibroblasts=c(fibroblasts,rep(NA,length(tcells)-length(fibroblasts))),
                                 dendritic=c(dendritic,rep(NA,length(tcells)-length(dendritic))),
                                 epithelial_malignanat=c(epithelial_malignant,rep(NA,length(tcells)-length(epithelial_malignant))),
                                 endothelial=c(endothelial,rep(NA,length(tcells)-length(endothelial))),
                                 other=c(other,rep(NA,length(tcells)-length(other))))


write.xlsx(celltype_assignments, "dataset_info/celltype_assignments.xlsx")

heatmap_data<-scored_data[,c("dataset","celltype","precog_score_pearson")]

colnames(heatmap_data)<-c("dataset","celltype","score")
heatmap_data$celltype[heatmap_data$celltype%in%tcells]<-"T Cells"
heatmap_data$celltype[heatmap_data$celltype%in%bcells]<-"B Cells"
heatmap_data$celltype[heatmap_data$celltype%in%mono_mac]<-"Macrophage w/ (Mono/Myeloid)"
heatmap_data$celltype[heatmap_data$celltype%in%mast]<-"Mast Cells"
heatmap_data$celltype[heatmap_data$celltype%in%nk]<-"NK Cells"
heatmap_data$celltype[heatmap_data$celltype%in%other_immune]<-"Other Immune"
heatmap_data$celltype[heatmap_data$celltype%in%fibroblasts]<-"Fibroblasts"
heatmap_data$celltype[heatmap_data$celltype%in%dendritic]<-"Dendritic"
heatmap_data$celltype[heatmap_data$celltype%in%epithelial_malignant]<-"Epithelial/Malignant"
heatmap_data$celltype[heatmap_data$celltype%in%endothelial]<-"Endothelial"
heatmap_data$celltype[heatmap_data$celltype%in%other]<-"Other"

heatmap_data<-heatmap_data[complete.cases(heatmap_data),]

collapsed<-aggregate(score~dataset+celltype,data=heatmap_data,FUN=mean)
heatmap_input<-reshape2::dcast(collapsed,dataset~celltype,mean)
heatmap_input<-column_to_rownames(heatmap_input, "dataset")
heatmap_input<-t(scale(t(heatmap_input)))
heatmap_input[is.na(heatmap_input)]<-0

order<-aggregate(collapsed$score, list(collapsed$celltype), median)
order<-rev(order[order(order[,2]),][,1])

heatmap_input<-heatmap_input[,order]
heatmap_input<-heatmap_input[dataset_info$dataset[order(dataset_info$precog_label)],]
heatmap_input<-t(heatmap_input)

h<-Heatmap(heatmap_input,
           cluster_rows = FALSE,
           cluster_columns = F,
           row_names_side = "left",
           col = circlize::colorRamp2(c(-1.5, 0, 1.5),c("darkblue","white","red")),
           heatmap_legend_param = list(title = "Scaled Prognostic Score", direction='horizontal'),
           column_names_rot = 45)

ComplexHeatmap::draw(h, heatmap_legend_side = "bottom")

gb = grid.grabExpr(ComplexHeatmap::draw(h, heatmap_legend_side = "bottom"))

b<-ggplot(collapsed, aes(x=celltype, y=score, fill=celltype)) +
  coord_flip() +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.size=0.5) +
  geom_jitter(position=position_jitter(0.1), size=1) +
  scale_fill_manual(values = get_palette("ucscgb",length(unique(collapsed$celltype)))) +
  xlab(NULL)+
  ylab("Average Prognostic Score") +
  #ggtitle("Average Celltype Prognostic Score Per Dataset") +
  theme_bw() +
  scale_x_discrete(limits=rev(order)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=10),
        axis.text.y = element_blank(),
        legend.position = "none") 


c<-ggarrange(b,nrow=2, heights=c(10,3.75))
ggarrange(gb,c,ncol=2, widths=c(3,2))
