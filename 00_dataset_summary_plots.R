# Armon Azizi
# Generate dataset summary plots for single cell data

library(ggplot2)
library(openxlsx)
library(ggpubr)

setwd("~/Bioinformatics/sc_prognosis/")

celltype_counts<-read.xlsx("outputs/celltype_counts.xlsx")
dataset_info<-read.xlsx("dataset_info/dataset_info.xlsx")

colors<-get_palette("Set1",11)[c(-4)]

# donut plot: number of cells per histology

plot_input<-data.frame(cancer_type=character(),
                       num_cells=numeric(),
                       stringsAsFactors = FALSE)
for(cancer_type in unique(dataset_info$precog_label)){
  datasets<-dataset_info$dataset[dataset_info$precog_label==cancer_type]
  num_cells<-sum(celltype_counts$count[celltype_counts$dataset%in%datasets])
  plot_input<-rbind(plot_input,data.frame(cancer_type=cancer_type,num_cells=num_cells))
}


plot_input$ymax = cumsum(plot_input$num_cells)
plot_input$ymin = c(0, head(plot_input$ymax, n=-1))

plot_input$label <- paste0(plot_input$num_cells)
plot_input$labelPosition <- (plot_input$ymax + plot_input$ymin) / 2

ggplot(plot_input, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=cancer_type)) +
  geom_rect(color="black") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  scale_fill_manual(values = colors) +
  geom_text(x=3.5, aes(y=labelPosition, label=label, fontface=2), size=3.5) +
  theme_void() +
  ggtitle("Number Of Cells Analyzed Per Cancer Type") +
  theme(plot.title = element_text(hjust = 0.5))


# donut plot: cells per platform

plot_input<-data.frame(platform=character(),
                       num_cells=numeric(),
                       stringsAsFactors = FALSE)
for(platform in unique(dataset_info$Platform)){
  datasets<-dataset_info$dataset[dataset_info$Platform==platform]
  num_cells<-sum(celltype_counts$count[celltype_counts$dataset%in%datasets])
  plot_input<-rbind(plot_input,data.frame(platform=platform,num_cells=num_cells))
}


plot_input$ymax = cumsum(plot_input$num_cells)
plot_input$ymin = c(0, head(plot_input$ymax, n=-1))

plot_input$label <- paste0(plot_input$num_cells)
plot_input$labelPosition <- (plot_input$ymax + plot_input$ymin) / 2

ggplot(plot_input, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=platform)) +
  geom_rect(color="black") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  scale_fill_manual(values = colors) +
  geom_text(x=3.5, aes(y=labelPosition, label=label, fontface=2), size=3.5) +
  theme_void() +
  ggtitle("Number Of Cells Analyzed Per Platform") +
  theme(plot.title = element_text(hjust = 0.5))


# stacked bar of number of cells per celltype in each dataset

plot_input<-celltype_counts
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

plot_input$celltype[plot_input$celltype%in%tcells]<-"T Cells"
plot_input$celltype[plot_input$celltype%in%bcells]<-"B Cells"
plot_input$celltype[plot_input$celltype%in%mono_mac]<-"Macrophage w/ (Mono/Myeloid)"
plot_input$celltype[plot_input$celltype%in%mast]<-"Mast Cells"
plot_input$celltype[plot_input$celltype%in%nk]<-"NK Cells"
plot_input$celltype[plot_input$celltype%in%other_immune]<-"Other Immune"
plot_input$celltype[plot_input$celltype%in%fibroblasts]<-"Fibroblasts"
plot_input$celltype[plot_input$celltype%in%dendritic]<-"Dendritic"
plot_input$celltype[plot_input$celltype%in%epithelial_malignant]<-"Epithelial/Malignant"
plot_input$celltype[plot_input$celltype%in%endothelial]<-"Endothelial"
plot_input$celltype[plot_input$celltype%in%other]<-"Other"

plot_input<-aggregate(count~dataset+celltype,data=plot_input,FUN=sum)

for(ds in unique(plot_input$dataset)){
  plot_input$count[plot_input$dataset==ds]<-plot_input$count[plot_input$dataset==ds]/sum(plot_input$count[plot_input$dataset==ds])
}

colors<-get_palette("ucscgb",length(unique(plot_input$celltype)))

ggplot(plot_input, aes(x=dataset,y=count)) +
  geom_bar(aes(fill=celltype),position="stack", stat="identity", color="black") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  ggtitle("Cell Type Distribution Per Dataset") +
  xlab("") +
  ylab("Fraction Of Cells") +
  theme(axis.text.x = element_text(angle = 45,hjust=1, size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust=0.5, size=15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15))
