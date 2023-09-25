# Armon Azizi
# aazizi@stanford.edu
#
# Single Cell PRECOG analysis functions

if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
if(!require(Seurat)){install.packages("Seurat")}
if(!require(reticulate)){install.packages("reticulate")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(viridis)){install.packages("viridis")}
if(!require(data.table)){install.packages("data.table")}
if(!require(tibble)){install.packages("tibble")}
if(!require(R.utils)){install.packages("R.utils")}
if(!require(fgsea)){BiocManager::install("fgsea")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(gplots)){install.packages("gplots")}
if(!require(ComplexHeatmap)){BiocManager::install("ComplexHeatmap")}
if(!require(ggpubr)){install.packages("ggpubr")}
if(!require(psych)){install.packages("psych")}
library(Seurat)
library(reticulate)
library(ggplot2)
library(viridis)
library(data.table)
library(tibble)
library(R.utils)
library(fgsea)
library(openxlsx)
library(gplots)
library(ComplexHeatmap)
library(ggpubr)
library(psych)
library(svd)

# Preprocess an SCrnaseq counts dataframe into a seurat object.
# Uses default settings; preps data for tSNE and uMAP.
#
# Data must be in dataframe format with colnames as cells and rownames as genes.
preprocess_seurat<-function(seurat_data, min_cells=5, min_features=200){
  seurat_data<- CreateSeuratObject(seurat_data, min.cells = min_cells, min.features = min_features)
  seurat_data<-NormalizeData(object=seurat_data)
  seurat_data<-FindVariableFeatures(object=seurat_data)
  seurat_data <- ScaleData(object=seurat_data)
  return(seurat_data)
}

# Run tSNE on preprocessed Seurat Object.
# Runs PCA amd then tSNE.
# Set findClusters to FALSE if cells have been previously annotated with set_cell_labels()
seurat_tsne<-function(processed_seurat, title, colors=NULL, findClusters=TRUE, recalculate=FALSE, file=NULL){
  processed_seurat <- RunPCA(object = processed_seurat, pc.genes = processed_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  
  if(findClusters){
    print("Finding Clusters...")
    # processed_seurat <- FindClusters(object = processed_seurat, reduction.type = "pca", dims.use = 1:10,
    #                                  resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc=recalculate)
    processed_seurat <- FindNeighbors(object = processed_seurat)
    processed_seurat <- FindClusters(object = processed_seurat)
  }
  
  print("Running tSNE...")
  processed_seurat <- RunTSNE(object = processed_seurat, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
  
  plot<-DimPlot(object = processed_seurat, reduction = "tsne", do.return=TRUE, cols=colors) +
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
  
  return(processed_seurat)
}

# Run uMAP on seurat object
# Requires seurat_tsne() to be run on object previously.
seurat_umap<-function(processed_seurat, title, colors=NULL, file=NULL, recalculate=TRUE){
  
  if(recalculate){
    processed_seurat <- RunUMAP(processed_seurat, dims = 1:10) 
  }
  
  plot<-DimPlot(processed_seurat, reduction.use = "umap", cols=colors) +
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
  
  return(processed_seurat)
}

# set labels for ecah cell for plotting.
# labels=dataframe with rownames as cell labels and one column of metadata
set_cell_labels<-function(processed_seurat, labels){
  processed_seurat <- AddMetaData(object = processed_seurat, metadata = labels)
  processed_seurat <- StashIdent(object = processed_seurat, save.name = "Cluster")
  processed_seurat <- SetIdent(object = processed_seurat, value = "celltype")
  return(processed_seurat)
}

get_precog_genes<-function(processed_seurat, precog_cancer_of_interest, num_genes, prognosis){
  subset_precog<-precog.data[rownames(processed_seurat),]
  subset_precog<-subset_precog[complete.cases(subset_precog),]
  
  if(prognosis=="favorable"){
    genes_of_interest<-rownames(subset_precog[order(subset_precog[,precog_cancer_of_interest]),][1:num_genes,])
  } else if(prognosis=="poor"){
    genes_of_interest<-rownames(subset_precog[rev(order(subset_precog[,precog_cancer_of_interest])),][1:num_genes,])
  }
  
  return(genes_of_interest)
}


### DEPRECATED ###
# score_cells_precog<-function(processed_seurat, precog_data, cancer_of_interest, z_cutoff=2, rank_space=FALSE){
#   precog_z_scores<-precog_data[,cancer_of_interest,drop=FALSE]
#   precog_z_scores<-precog_z_scores[abs(precog_z_scores[,cancer_of_interest])>z_cutoff,,drop=FALSE]
#   norm_expression_data<-as.data.frame(t(FetchData(processed_seurat, intersect(rownames(precog_z_scores), rownames(processed_seurat)))))
#   precog_z_scores<-precog_z_scores[rownames(norm_expression_data),,drop=FALSE]
#   #norm_expression_data[norm_expression_data==0]<-NA
#   correlations<-cor(precog_z_scores,norm_expression_data)
#   scores<-as.data.frame(t(correlations))
#   colnames(scores)<-c("precog_score")
#   
#   if(rank_space){
#     scores<-scores[order(scores$precog_score),,drop=FALSE]
#     scores$precog_score<-1:length(scores$precog_score)
#   }
#   
#   return(scores)
# }

score_cells_precog<-function(processed_seurat, precog_data, cancer_of_interest, z_cutoff=2, rank_space=FALSE){
  
  # Get precog data
  precog_z_scores<-precog_data[,cancer_of_interest,drop=FALSE]
  precog_z_scores<-precog_z_scores[abs(precog_z_scores[,cancer_of_interest])>z_cutoff,,drop=FALSE]
  
  # Get normalized scaled
  message("Extracting Expression Data...")
  norm_expression_data<-as.data.frame(t(FetchData(processed_seurat, intersect(rownames(precog_z_scores), rownames(processed_seurat)))))
  
  # Calculate z-scores for each cell
  # Dont take dropout or nonexpressed genes into consideration for z-score
  message("Calculating z-scores...")
  expression_zscores<-norm_expression_data
  expression_zscores[expression_zscores==0]<-NA
  expression_zscores<-as.data.frame(apply(expression_zscores, 2, function(x){scale(x)}))
  rownames(expression_zscores)<-rownames(norm_expression_data)
  precog_z_scores<-precog_z_scores[rownames(expression_zscores),,drop=FALSE]
  
  scores<-data.frame()
  
  message("Calculating Scores")
  count<-0
  
  # For each cell, calculate correlation between precog z scores and expression z scores
  for(cell in colnames(expression_zscores)){
    exp_vector<-expression_zscores[,cell,drop=FALSE]
    exp_vector<-exp_vector[!is.na(exp_vector[,cell]),,drop=FALSE]
    zscores<-precog_z_scores[rownames(exp_vector),,drop=FALSE]
    scores[cell,"precog_score"]<-cor(zscores[,1],exp_vector[,1])
    
    count<-count+1
    if(count%%1000==0){message(paste("Calculated Scores For:",as.character(count), "Cells"))}
  }
  
  if(rank_space){
    scores<-scores[order(scores$precog_score),,drop=FALSE]
    scores$precog_score<-1:length(scores$precog_score)
  }
  
  return(scores)
}


# smooth_scores<-function(processed_seurat, scores, rank_space=FALSE, subsample_size=1000){
#   smoothed<-cytoTRACE_modify(as.data.frame(t(FetchData(processed_seurat,vars = rownames(processed_seurat)))), enableFast = TRUE, subsamplesize=subsample_size, metric = scores, rank_space=rank_space)
#   cytotrace_scores<-as.data.frame(smoothed$CytoTRACE)
#   colnames(cytotrace_scores)<-c("cytotrace")
#   cytotrace_scores<-cytotrace_scores[colnames(processed_seurat),,drop=FALSE]
#   return(cytotrace_scores)
# }

score_cells_cytotrace<-function(processed_seurat,
                                subsample_size=1000,
                                rank_space=TRUE){
  
  cytotrace_results<-cytoTRACE(as.data.frame(t(FetchData(processed_seurat,vars = rownames(processed_seurat)))), enableFast = TRUE, subsamplesize=subsample_size, rank_space=rank_space)
  cytotrace_scores<-as.data.frame(cytotrace_results$CytoTRACE)
  colnames(cytotrace_scores)<-c("cytotrace")
  cytotrace_scores<-cytotrace_scores[colnames(processed_seurat),,drop=FALSE]
  return(cytotrace_scores)
  
}

# Plots the mean of expression of poor or favorable prognostic genes from PRECOG.
plot_combination_expression<-function(processed_seurat, 
                                      precog_cancer_of_interest, 
                                      num_genes=10, 
                                      prognosis="favorable", 
                                      reduction="tsne",
                                      title=NULL,
                                      return_average_score=FALSE,
                                      file=NULL){
  
  genes_of_interest<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis)
  
  print("Genes of interest:")
  print(genes_of_interest)
  
  if(return_average_score){
    return(mean(rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)))
  }
  
  # Get precog score and scale down outliers so that they don't wash out colorscale
  normalized_values<-log2(rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)+1)
  outliers<-boxplot(normalized_values, plot=FALSE)$out
  max_val<-max(normalized_values[!(normalized_values %in% outliers)])
  normalized_values[(normalized_values %in% outliers)]<-max_val
  
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = normalized_values, 
                                  col.name = "expScore")
  
  plot <- FeaturePlot(processed_seurat, "expScore", 
                      reduction=reduction)
  
  plot<-plot + 
    scale_color_viridis("PRECOG Score") + 
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
  
}

# Plots precog score
plot_precog_score_expression<-function(processed_seurat, 
                                       precog_cancer_of_interest, 
                                       z_cutoff=2,
                                       rank_space=FALSE,
                                       reduction="tsne",
                                       title=NULL,
                                       file=NULL){
  
  scores<-score_cells_precog(processed_seurat, 
                             cancer_of_interest = precog_cancer_of_interest,
                             precog_data = precog.data,
                             z_cutoff = z_cutoff,
                             rank_space = rank_space)
  
  scores<-scores[colnames(processed_seurat),]
  #scores[scores>0]<-scores[scores>0]/max(scores)
  #scores[scores<0]<-(-1*scores[scores<0]/min(scores))
  
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = scores, 
                                  col.name = "expScore")
  
  plot <- FeaturePlot(processed_seurat, "expScore", 
                      reduction=reduction)
  
  plot<-plot + 
    #scale_colour_gradient2("PRECOG Score", low="#440154FF", mid="grey50", high="#FDE725FF", midpoint=mean(c(max(scores),min(scores)))) +
    #scale_colour_gradient2("PRECOG Score", low="blue", mid="grey50", high="red", midpoint=mean(c(max(scores),min(scores)))) +
    #scale_colour_gradient2("PRECOG Score", low="navy", mid="grey50", high="red", midpoint=mean(c(max(scores),min(scores)))) +
    scale_color_viridis("PRECOG Score") +
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

# Plots the mean of expression of a gene set.
# applies precog score function to genes
plot_gene_set_expression<-function(processed_seurat, 
                                   gene_list=NULL,
                                   reduction="tsne",
                                   title=NULL,
                                   return_average_score=FALSE,
                                   file=NULL){
  
  genes_of_interest<-gene_list
  
  if(return_average_score){
    return(mean(rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)))
  }
  
  # Get precog score and scale down outliers so that they don't wash out colorscale
  normalized_values<-log2(rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)+1)
  outliers<-boxplot(normalized_values, plot=FALSE)$out
  max_val<-max(normalized_values[!(normalized_values %in% outliers)])
  normalized_values[(normalized_values %in% outliers)]<-max_val
  
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = normalized_values, 
                                  col.name = "expScore")
  
  plot <- FeaturePlot(processed_seurat, "expScore", 
                      reduction=reduction)
  
  plot<-plot + 
    scale_color_viridis("Expression Score") + 
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
  
}

plot_single_gene_expression<-function(processed_seurat, 
                                      gene=NULL, 
                                      reduction="tsne",
                                      title=NULL,
                                      file=NULL){
  
  plot <- FeaturePlot(processed_seurat, gene,
                      reduction=reduction)
  plot<-plot + 
    scale_color_viridis(paste(gene,"Expression")) +
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

plot_combination_violin<-function(processed_seurat, 
                                  precog_cancer_of_interest, 
                                  num_genes=10, 
                                  prognosis="favorable", 
                                  title=NULL,
                                  file=NULL){
  
  genes_of_interest<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis)
  
  print("Genes of interest:")
  print(genes_of_interest)
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = (rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)), 
                                  col.name = "expScore")
  
  plot<-VlnPlot(processed_seurat,c("expScore"), x.lab.rot = TRUE, do.return=TRUE)
  
  plot<-plot+ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

plot_gene_set_violin<-function(processed_seurat, 
                               genes, 
                               num_genes=10, 
                               prognosis="favorable", 
                               title=NULL,
                               file=NULL){
  
  genes_of_interest<-genes
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = (rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)), 
                                  col.name = "expScore")
  
  plot<-VlnPlot(processed_seurat,c("expScore"), x.lab.rot = TRUE, do.return=TRUE)
  
  plot<-plot+ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

plot_precog_violin<-function(processed_seurat, 
                             precog_cancer_of_interest, 
                             z_cutoff=2,
                             rank_space=FALSE,
                             title=NULL,
                             file=NULL){
  
  scores<-score_cells_precog(processed_seurat, 
                             cancer_of_interest = precog_cancer_of_interest,
                             precog_data = precog.data,
                             z_cutoff = z_cutoff,
                             rank_space = rank_space)
  
  scores<-scores[colnames(processed_seurat),]
  
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = scores, 
                                  col.name = "expScore")
  
  plot<-VlnPlot(processed_seurat,c("expScore"), x.lab.rot = TRUE, do.return=TRUE)
  
  plot<-plot+ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

plot_precog_histogram<-function(processed_seurat, 
                                precog_cancer_of_interest, 
                                z_cutoff=2,
                                rank_space=FALSE,
                                title=NULL,
                                file=NULL){
  
  scores<-score_cells_precog(processed_seurat, 
                             cancer_of_interest = precog_cancer_of_interest,
                             precog_data = precog.data,
                             z_cutoff = z_cutoff,
                             rank_space = rank_space)
  
  scores<-scores[colnames(processed_seurat),]
  
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = scores, 
                                  col.name = "expScore")
  
  plot<-ggplot(processed_seurat@meta.data, aes(x=expScore)) + 
    geom_histogram(bins=100, color="black",fill="darkgrey") +
    ylab("Frequency") +
    xlab("PRECOG Score") +
    ggtitle(title)
  
  # ggplot(processed_seurat@meta.data[processed_seurat@meta.data[,"expScore"]!=0,], aes(x=expScore)) + 
  #   geom_histogram(bins=100, color="black",fill="darkgrey") +
  #   ylab("Frequency") +
  #   xlab("PRECOG Score") +
  #   ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

plot_metadata_pca<-function(processed_seurat,
                            feature,
                            colors,
                            title=NULL){
  
  features<-processed_seurat@meta.data[,feature,drop=FALSE]
  pcs<-processed_seurat@reductions$pca@cell.embeddings[,c("PC_1","PC_2")]
  
  plot_input<-cbind(pcs,features)
  
  plot<-ggplot(plot_input, aes_string(x="PC_1",y="PC_2",color=feature)) +
    geom_point() +
    scale_color_manual(feature, values=colors) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          legend.text = element_text(color="black", size=12),
          legend.title = element_text(color="black", size=15),
          title=element_text(color="black", size=15)) +
    ggtitle(title)
  print(plot)
}


plot_gene_set_pca<-function(processed_seurat, 
                            gene_list=NULL,
                            reduction="tsne",
                            title=NULL,
                            return_average_score=FALSE,
                            file=NULL){
  
  genes_of_interest<-gene_list
  
  if(return_average_score){
    return(mean(rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)))
  }
  
  # Get precog score and scale down outliers so that they don't wash out colorscale
  normalized_values<-log2(rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)+1)
  outliers<-boxplot(normalized_values, plot=FALSE)$out
  max_val<-max(normalized_values[!(normalized_values %in% outliers)])
  normalized_values[(normalized_values %in% outliers)]<-max_val
  
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = normalized_values, 
                                  col.name = "expScore")
  
  features<-processed_seurat@meta.data[,feature,drop=FALSE]
  pcs<-processed_seurat@reductions$pca@cell.embeddings[,c("PC_1","PC_2")]
  
  plot_input<-cbind(pcs,features)
  
  plot<-ggplot(plot_input, aes_string(x="PC_1",y="PC_2",color=feature)) +
    geom_point() +
    scale_color_manual(feature, values=colors) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(color="black", size=12),
          axis.text.y = element_text(color="black", size=12),
          legend.text = element_text(color="black", size=12),
          legend.title = element_text(color="black", size=15),
          title=element_text(color="black", size=15)) +
    ggtitle(title)
  print(plot)
}




# Finds correlations between gene expressions and precog score
find_precog_gene_correlations<-function(processed_seurat,
                                        precog_slot="precog_score",
                                        z_cutoff=2,
                                        metric="pearson",
                                        rank_space=FALSE){
  
  scores<-processed_seurat@meta.data[,precog_slot,drop=FALSE]
  
  norm_expression_data<-as.data.frame(FetchData(processed_seurat, rownames(processed_seurat)))
  
  scores<-scores[rownames(norm_expression_data),,drop=FALSE]
  
  z_scores<-expression_zscores<-as.data.frame(t(apply(norm_expression_data, 1, function(x){scale(x)})))
  colnames(z_scores)<-colnames(norm_expression_data)
  
  scores[,precog_slot]<-scale(scores[,precog_slot])
  
  correlations<-psych::corr.test(norm_expression_data, scores, ci = FALSE, method=metric)
  
  result<-cbind(correlations$r, correlations$p)
  colnames(result)<-c("pearson_r", "p_val")
  
  return(result)
}

# Finds correlations between gene expressions and cytotrace score
find_cytotrace_gene_correlations<-function(processed_seurat,
                                           cytotrace_scores){
  
  scores<-cytotrace_scores
  
  norm_expression_data<-as.data.frame(FetchData(processed_seurat, rownames(processed_seurat)))
  
  scores<-scores[rownames(norm_expression_data),,drop=FALSE]
  
  correlations<-psych::corr.test(norm_expression_data, scores, ci = FALSE)
  
  result<-cbind(correlations$r, correlations$p)
  colnames(result)<-c("pearson_r", "p_val")
  
  return(result)
}

plot_combination_ridge<-function(processed_seurat, 
                                 precog_cancer_of_interest, 
                                 num_genes=10, 
                                 prognosis="favorable", 
                                 title=NULL,
                                 file=NULL){
  
  genes_of_interest<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis)
  
  print("Genes of interest:")
  print(genes_of_interest)
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = (rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)), 
                                  col.name = "expScore")
  
  plot<-RidgePlot(processed_seurat,c("expScore"), x.lab.rot = TRUE, do.return=TRUE)
  
  plot<-plot+ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}

precog_score_histogram<-function(processed_seurat, 
                                 precog_cancer_of_interest, 
                                 num_genes=10, 
                                 prognosis="favorable", 
                                 title=NULL,
                                 file=NULL){
  
  genes_of_interest<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis)
  
  print("Genes of interest:")
  print(genes_of_interest)
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = (rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)), 
                                  col.name = "expScore")
  
  plot<-ggplot(processed_seurat@meta.data, aes(x=expScore)) + 
    geom_histogram(bins=100, color="black",fill="darkgrey") +
    ylab("Frequency") +
    xlab("PRECOG Score") +
    ggtitle(title)
  
  # ggplot(processed_seurat@meta.data[processed_seurat@meta.data[,"expScore"]!=0,], aes(x=expScore)) + 
  #   geom_histogram(bins=100, color="black",fill="darkgrey") +
  #   ylab("Frequency") +
  #   xlab("PRECOG Score") +
  #   ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
  
}

add_precog_feature<-function(processed_seurat,
                             precog_cancer_of_interest, 
                             num_genes=10,
                             prognosis="both",
                             label="precog_score"){
  
  genes_of_interest<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis)
  
  # processed_seurat <- AddMetaData(object = processed_seurat, 
  #                                 metadata = (rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)), 
  #                                 col.name = label)
  
  score_data<-rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)
  score_data<-matrix(score_data, nrow=1,ncol=length(score_data))
  rownames(score_data)<-label
  
  processed_seurat@assays$RNA@data<-rbind(processed_seurat@assays$RNA@data, score_data)
  
  return(processed_seurat)
}


seurat_precog_heatmap<-function(processed_seurat,
                                precog_cancer_of_interest, 
                                num_genes=10,
                                gene_list=NULL,
                                gene_annot=NULL,
                                num_cells=NULL,
                                cell_order=NULL,
                                cluster_cells=TRUE,
                                prognosis="both",
                                show_gene_names=FALSE,
                                annotate=TRUE,
                                title=NULL){
  
  if(!is.null(gene_list)){
    genes<-gene_list
  }else{
    poor_genes<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis="poor")
    favorable_genes<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis="favorable")
    
    if(prognosis=="poor"){
      genes<-poor_genes
    } else if(prognosis=="favorable"){
      genes<-favorable_genes
    } else {
      genes<-c(poor_genes,favorable_genes) 
    }
  }
  
  annotations<-as.data.frame(processed_seurat@active.ident)
  
  norm_data<-processed_seurat@assays$RNA@data
  
  norm_data<-norm_data[genes,]
  if(!is.null(num_cells)){
    norm_data<-norm_data[,sample(ncol(norm_data), num_cells)]
  }
  norm_data<-log2(norm_data+1)
  
  if(!is.null(cell_order)){
    norm_data<-norm_data[,cell_order]
  }
  
  annotations<-annotations[colnames(norm_data),]
  
  
  colors<-colorRampPalette(c("darkblue","yellow"))
  cell_colors<-as.vector(get_palette(palette="ucscgb", k=length(levels(annotations))))
  names(cell_colors)<-levels(annotations)
  
  if(!is.null(gene_annot)){
    annot_pos<-match(gene_annot, rownames(norm_data))
    ga<-rowAnnotation(foo = anno_mark(at = annot_pos, labels = gene_annot))
    show_gene_names<-FALSE
  }else{
    ga<-NULL
  }
  
  if(annotate){
    Heatmap(as.matrix(norm_data), 
            col=colors(50),
            cluster_columns = cluster_cells,
            show_row_names=show_gene_names, 
            show_column_names=FALSE,
            top_annotation=HeatmapAnnotation(data.frame(CellType = annotations), col=list(CellType=cell_colors)),
            heatmap_legend_param = list(title = "Expression"),
            column_title=title,
            right_annotation=ga)
  }else{
    Heatmap(as.matrix(norm_data), 
            col=colors(50),
            cluster_columns = cluster_cells,
            show_row_names=show_gene_names, 
            show_column_names=FALSE,
            #top_annotation=HeatmapAnnotation(data.frame(CellType = annotations), col=list(CellType=cell_colors)),
            heatmap_legend_param = list(title = "Expression"),
            column_title=title,
            right_annotation=ga)
  }
  
}



precog_cells_diff_exp<-function(processed_seurat, 
                                precog_cancer_of_interest, 
                                num_genes=10, 
                                prognosis="favorable",
                                score_cutoff=NULL){
  
  genes_of_interest<-get_precog_genes(processed_seurat, precog_cancer_of_interest, num_genes, prognosis)
  
  print("Genes of interest:")
  print(genes_of_interest)
  
  processed_seurat <- AddMetaData(object = processed_seurat, 
                                  metadata = (rowSums(FetchData(processed_seurat, genes_of_interest))/length(genes_of_interest)), 
                                  col.name = "expScore")
  
  # create new column in the metadata for Neg/Pos status of gene
  processed_seurat@meta.data$precog_status <- 'Neg'
  processed_seurat@meta.data$precog_status[which(processed_seurat@meta.data[,"expScore"] >= score_cutoff)] <- 'Pos'
  
  # subset only the positive cells; change the ident first.
  processed_seurat <- SetIdent(object = processed_seurat, value = "precog_status")
  
  return(FindMarkers(object = processed_seurat, ident.1 = "Pos", ident.2 = "Neg", logfc.threshold=0))
}

cluster_diff_exp<-function(processed_seurat, cluster1=NULL,cluster2=NULL){
  return(FindMarkers(object = processed_seurat, ident.1 = cluster1, ident.2 = cluster2, logfc.threshold=0))
}


differential_gene_volcano<-function(findmarkers_output){
  findmarkers_output$p_val_adj<--log10(findmarkers_output$p_val_adj+.Machine$double.xmin)
  ggplot(findmarkers_output, aes(x=avg_logFC, y=p_val_adj))+
    geom_point() +
    theme_bw()
}


run_fgsea<-function(ranked_genes){
  
  ranked_genes$SCORE<--log10(ranked_genes$p_val+.Machine$double.xmin)*sign(ranked_genes$avg_logFC)
  ranked_genes<-ranked_genes[,"SCORE",drop=FALSE]
  ranked_genes<-rownames_to_column(ranked_genes)
  colnames(ranked_genes)<-c("Symbol","SCORE")
  ranked_genes <- setNames(ranked_genes$SCORE, ranked_genes$Symbol)
  
  # Load pathways
  pathways <- gmtPathways("/Users/majetilab/lungtmi_server/LTMI/gene_sets/msigdb_c2_c5_paper_gene_set.gmt")
  
  fgseaRes <- fgsea(pathways, ranked_genes, minSize=15, maxSize=500, nperm=1000)
  
  return(fgseaRes)
}

run_corr_fgsea<-function(ranked_genes){
  
  ranked_genes$SCORE<-ranked_genes$pearson_r
  ranked_genes<-ranked_genes[,"SCORE",drop=FALSE]
  ranked_genes<-rownames_to_column(ranked_genes)
  colnames(ranked_genes)<-c("Symbol","SCORE")
  ranked_genes <- setNames(ranked_genes$SCORE, ranked_genes$Symbol)
  
  # Load pathways
  pathways <- gmtPathways("/Users/majetilab/lungtmi_server/LTMI/gene_sets/msigdb_c2_c5_paper_gene_set.gmt")
  
  fgseaRes <- fgsea(pathways, ranked_genes, minSize=15, maxSize=500, nperm=1000)
  
  return(fgseaRes)
}

run_precog_fgsea<-function(precog_cancer_of_interest, prognosis="favorable"){
  
  ranked_genes<-precog.data[,precog_cancer_of_interest,drop=FALSE]
  ranked_genes<-rownames_to_column(ranked_genes)
  colnames(ranked_genes)<-c("Symbol","SCORE")
  if(prognosis=="favorable"){
    ranked_genes$SCORE<-ranked_genes$SCORE*-1
  }
  ranked_genes <- setNames(ranked_genes$SCORE, ranked_genes$Symbol)
  
  # Load pathways
  pathways <- gmtPathways("/Users/majetilab/lungtmi_server/LTMI/gene_sets/msigdb_c2_c5_paper_gene_set.gmt")
  
  fgseaRes <- fgsea(pathways, ranked_genes, minSize=15, maxSize=500, nperm=10000)
  
  return(fgseaRes)
}

find_purity<-function(processed_seurat, features=NULL){
  
  if(is.null(features)){
    features<-rownames(processed_seurat)
  }
  
  processed_seurat.subset <- GetAssayData(processed_seurat, slot="data")[features, ]
  processed_seurat.subset <- CreateSeuratObject(processed_seurat.subset)
  processed_seurat.subset<-NormalizeData(object=processed_seurat.subset)
  processed_seurat.subset<-FindVariableFeatures(object=processed_seurat.subset)
  processed_seurat.subset <- ScaleData(processed_seurat.subset)
  processed_seurat.subset <- RunPCA(object = processed_seurat.subset)
  processed_seurat.subset <- FindNeighbors(object = processed_seurat.subset)
  processed_seurat.subset <- FindClusters(object = processed_seurat.subset)
  
  cluster_cell_assignments<-as.data.frame(cbind(processed_seurat.subset@active.ident, processed_seurat@active.ident))
  colnames(cluster_cell_assignments)<-c("cluster", "cell_type")
  purity_df<-data.frame(cluster=unique(cluster_cell_assignments$cluster))
  
  for(i in rownames(purity_df)){
    cluster<-purity_df[i,"cluster"]
    purity_df[i,"assigned_cell_type"]<-names(sort(table(cluster_cell_assignments[cluster_cell_assignments$cluster==cluster,]$cell_type),decreasing=TRUE))[1]
    purity_df[i,"total_cells"]<-sum(cluster_cell_assignments$cluster==cluster)
    purity_df[i,"correct_cells"]<-sum(cluster_cell_assignments$cluster==cluster&cluster_cell_assignments$cell_type==purity_df[i,"assigned_cell_type"])
  }
  
  purity<-sum(purity_df$correct_cells)/sum(purity_df$total_cells)
  
  return(purity)
}

# Plots cytotrace
plot_cytotrace<-function(processed_seurat,
                         reduction="tsne",
                         subsample_size=1000,
                         title=NULL,
                         file=NULL){
  
  cytotrace_results<-cytoTRACE(as.data.frame(t(FetchData(processed_seurat,vars = rownames(processed_seurat)))), enableFast = TRUE, subsamplesize=subsample_size)
  cytotrace_scores<-as.data.frame(cytotrace_results$CytoTRACE)
  colnames(cytotrace_scores)<-c("cytotrace")
  cytotrace_scores<-cytotrace_scores[colnames(processed_seurat),]
  processed_seurat<-AddMetaData(object = processed_seurat, metadata = cytotrace_scores, col.name = "cytotrace")
  plot<-FeaturePlot(processed_seurat, "cytotrace", reduction=reduction) +
    scale_color_viridis("CytoTRACE Score") +
    ggtitle(title)
  
  if(is.null(file)){
    print(plot)
  } else {
    ggsave(file,plot=plot)
  }
}


# rename celltypes
rename_celltypes<-function(df, column){
  
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
  
  df[,column][df[,column]%in%tcells]<-"T Cells"
  df[,column][df[,column]%in%bcells]<-"B Cells"
  df[,column][df[,column]%in%mono_mac]<-"Macrophage w/ (Mono/Myeloid)"
  df[,column][df[,column]%in%mast]<-"Mast Cells"
  df[,column][df[,column]%in%nk]<-"NK Cells"
  df[,column][df[,column]%in%other_immune]<-"Other Immune"
  df[,column][df[,column]%in%fibroblasts]<-"Fibroblasts"
  df[,column][df[,column]%in%dendritic]<-"Dendritic"
  df[,column][df[,column]%in%epithelial_malignant]<-"Epithelial/Malignant"
  df[,column][df[,column]%in%endothelial]<-"Endothelial"
  df[,column][df[,column]%in%other]<-"Other"
  
  return(df)
}








preprocess_tirosh<-function(path){
  tirosh<-fread(path,sep="\t",header=TRUE,data.table=FALSE)
  metadata<-as.data.frame(t(tirosh[1:3,]), stringsAsFactors=FALSE)
  metadata<-metadata[2:length(rownames(metadata)),]
  colnames(metadata)<-c("tumor","malignant_status","cell_type")
  metadata$malignant_status[metadata$malignant_status==" 1"]<-"Non-Malignant"
  metadata$malignant_status[metadata$malignant_status==" 2"]<-"Malignant"
  metadata$malignant_status[metadata$malignant_status==" 0"]<-"Unresolved"
  metadata$cell_type[metadata$malignant_status=="Malignant"]<-"Malignant"
  metadata$cell_type[metadata$cell_type==" 0"]<-"Unknown"
  metadata$cell_type[metadata$cell_type==" 1"]<-"T-Cell"
  metadata$cell_type[metadata$cell_type==" 2"]<-"B-Cell"
  metadata$cell_type[metadata$cell_type==" 3"]<-"Macrophage"
  metadata$cell_type[metadata$cell_type==" 4"]<-"Endothelial"
  metadata$cell_type[metadata$cell_type==" 5"]<-"Cancer-Associated Fibroblast"
  metadata$cell_type[metadata$cell_type==" 6"]<-"NK"
  rownames(metadata)<-gsub("-","_",rownames(metadata))
  
  tirosh<-tirosh[-(1:3),]
  tirosh<-tirosh[tirosh$Cell!="MARCH1",]
  tirosh<-tirosh[tirosh$Cell!="MARCH2",]
  rownames(tirosh)<-NULL
  tirosh<-column_to_rownames(tirosh, "Cell")
  colnames(tirosh)<-gsub("-","_",colnames(tirosh))
  
  return(list("data" = tirosh, "metadata" = metadata))
}

preprocess_hnscc<-function(path){
  hnscc<-fread(path,sep="\t",header=TRUE,data.table=FALSE)
  metadata<-as.data.frame(t(hnscc[1:5,]), stringsAsFactors=FALSE)
  metadata<-metadata[2:length(rownames(metadata)),]
  colnames(metadata)<-c("enzyme","lymph_node","malignant","non-malignant","cell_type")
  metadata$cell_type[metadata$cell_type=="0"]<-"Malignant"
  metadata$lymph_node[metadata$lymph_node=="0"]<-"Non-Lymph Node"
  metadata$lymph_node[metadata$lymph_node=="1"]<-"Lymph Node"
  metadata$malignant[metadata$malignant=="0"]<-"Non Malignant"
  metadata$malignant[metadata$malignant=="1"]<-"Malignant"
  
  hnscc<-hnscc[-(1:5),]
  hnscc<-hnscc[hnscc$V1!="MARCH1",]
  hnscc<-hnscc[hnscc$V1!="MARCH2",]
  rownames(hnscc)<-NULL
  hnscc<-column_to_rownames(hnscc, "V1")
  rownames(hnscc)<-gsub("'","",rownames(hnscc))
  
  hnscc[] <- lapply(hnscc, as.numeric)
  
  return(list("data" = hnscc, "metadata" = metadata))
}

preprocess_li_colorectal<-function(path){
  colo<-as.data.frame(fread(path,sep=",",header=TRUE,data.table=FALSE))
  genes<-sapply(strsplit(colo$V1,split = "_"), `[`, 2)
  colo<-stats::aggregate(colo[,2:ncol(colo)], by=list(genes), mean)
  colo<-column_to_rownames(colo, var = "Group.1")
  
  metadata<-as.data.frame(t(as.data.frame(strsplit(colnames(colo),split = "__"))), stringsAsFactors = FALSE)
  rownames(metadata)<-colnames(colo)
  metadata<-metadata[,2,drop=FALSE]
  colnames(metadata)<-c("celltype")
  metadata$celltype[metadata$celltype=="NA"]<-"Unknown"
  
  return(list("data" = colo, "metadata" = metadata))
}


reorder_dataframe<-function(df,reference,sample_order, priority_features=NULL, return_assignments=FALSE){
  sample_order<-as.factor(sample_order)
  colnames(reference)<-c("CellType")
  reference$CellType<-reorder(reference$CellType, new.order=sample_order)
  so<-rownames(reference[order(reference$CellType),,drop=FALSE])
  
  
  collapsed_counts<-data.frame(matrix(nrow=nrow(df),ncol=0))
  
  for(celltype in sample_order){
    message(celltype)
    x<-colnames(collapsed_counts)
    collapsed_counts<-cbind(collapsed_counts,rowMeans(df[,rownames(reference)[reference$CellType==celltype],drop=FALSE]))
    colnames(collapsed_counts)<-c(x,celltype)
  }
  assignments<-colnames(collapsed_counts)[max.col(collapsed_counts)]
  
  if(return_assignments){return(data.frame(gene=rownames(collapsed_counts),CellType=assignments))}
  
  feature_order<-c()
  for(celltype in colnames(collapsed_counts)){
    regions<-rownames(collapsed_counts)[assignments==celltype]
    regions<-regions[rev(order(collapsed_counts[regions,celltype,drop=FALSE]))]
    
    if(!is.null(priority_features)){
      regions<-c(regions[regions%in%priority_features],regions[!regions%in%priority_features])
    }
    
    feature_order<-c(feature_order,regions)
  }
  return(df[feature_order,so])
}


smooth_scores<-function(counts, scores, num_sim_cells=10, ALPHA=1){
  
  # generate similarity matrix
  temp_mtx<-counts
  n_expr <- rowSums(temp_mtx > 0);
  temp_mtx_filt <- temp_mtx[n_expr >= 0.05 * ncol(temp_mtx),];
  vars <- apply(temp_mtx_filt, 1, var);
  means <- apply(temp_mtx_filt, 1, mean);
  disp <- vars / means;
  last_disp <- tail(sort(disp), 1000)[1];
  temp_mtx_filt <- temp_mtx_filt[disp >= last_disp,];
  sim_mtx<-HiClimR::fastCor(temp_mtx_filt)
  
  # Only use cells which are in top 10 similar cells per cell
  sim_mtx[which(sim_mtx < 0)] <- 0;
  sim_mtx<-t(apply(sim_mtx, 1, function(x) {x[x<rev(sort(x))[num_sim_cells]]<-0; return(x)}))
  
  temp <- sim_mtx
  sim_mtx <- sim_mtx / rowSums(sim_mtx);
  sim_mtx[which(rowSums(temp)==0),] <- 0
  
  
  smoothed_scores <- ALPHA * (sim_mtx %*% scores) + (1 - ALPHA) * scores;
}

# Do pca on sparce matrix using svd propack function
# modified from bogdan code
sc_propack <- function(x,n) {
  
  # filter genes
  use_genes <- which(colSums(x) > 1)
  m <- x[,use_genes]
  
  # scale data
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- sweep(m, 2, colMeans(m), '-')
  m <- sweep(m, 2, apply(m, 2, sd), '/')
  
  # compute pcs
  ppk<-propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  return(pca)
}

# do seurat pca
seurat_pca<-function(x,n){
  seurat_data<- CreateSeuratObject(x, min.cells = 5, min.features = 200)
  seurat_data<-NormalizeData(object=seurat_data, verbose = FALSE)
  seurat_data<-FindVariableFeatures(object=seurat_data, verbose = FALSE)
  seurat_data<-ScaleData(object=seurat_data, verbose = FALSE)
  seurat_data<-RunPCA(object = seurat_data, verbose=FALSE)
  return(seurat_data@reductions$pca@cell.embeddings[,1:n])
}


get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,unique(quantile(mean,seq(0.1,1,0.05)),Inf))))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}
