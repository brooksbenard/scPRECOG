# Armon Azizi (aazizi@stanford.edu)
#
#' Run Prognostic Scoring Algorithm
#'
#' The run_prognostic_scoring function assigns prognostic scores to single cells using a list of genes associated with outcomes.
#' Input is a single cell expression matrix, genes are rows, cells are columns.
#' Secondary input is a dataframe of genes and associated prognostic scores for each gene.
#' @param expression matrix of gene expression values where columns are cells and rows are genes
#' @param prognostic_data dataframe of gene scores. Rownames are genes, and the singular column is prognostic value
#' @param z_cutoff cutoff for which genes to include in score, based on z-score value (absolute z-score value)
#' @param normalize_gene_values if prognostic_data is not in z-score space, normalize it. 
#' @param cor_method default pearson, can be "pearson" or "spearman"
#' @return a dataframe of single cell labels and prognostic scores.
run_prognostic_scoring<-function(expression, prognostic_data, z_cutoff=2, normalize_gene_values=FALSE, cor_method="pearson"){
  
  # subset gene data to genes that are above z-score threshold.
  message("Finding Significant Genes...")
  colnames(prognostic_data)<-c("prognostic_value")
  prognostic_data<-prognostic_data[abs(prognostic_data[,"prognostic_value"])>z_cutoff,,drop=FALSE]
  
  # if needed, convert scores to z-space
  if(normalize_gene_values){
    prognostic_data$prognostic_value<-scale(prognostic_data$prognostic_value)
  }
  
  
  # if not log2 put expression into log2 space
  if(max(expression)>50){
    expression<-log2(expression+1)
  }
  
  
  # Subset expression data to relevant genes
  expression_data<-expression[intersect(rownames(expression), rownames(prognostic_data)),]
  
  # Calculate z-scores for each cell
  # Dont take dropout or nonexpressed genes into consideration for z-score
  # 1. convert 0s to NA so they are not taken into consideration when doing z-score transform
  # 2. Apply scale function across all columns to convert expression values to z score.
  message("Calculating z-scores...")
  expression_zscores<-expression_data
  expression_zscores[expression_zscores==0]<-NA
  expression_zscores<-as.data.frame(scale(expression_zscores))
  rownames(expression_zscores)<-rownames(expression_data)
  prognostic_data<-prognostic_data[rownames(expression_zscores),,drop=FALSE]
  
  
  
  scores<-data.frame()
  
  message("Calculating Scores")
  count<-0
  
  # For each cell, calculate pearson correlation between prognostic z scores and expression z scores
  # Exclude NA values from score (these were initally nonexpressed or dropout)
  for(cell in colnames(expression_zscores)){
    exp_vector<-expression_zscores[,cell,drop=FALSE]
    exp_vector<-exp_vector[!is.na(exp_vector[,cell]),,drop=FALSE]
    zscores<-prognostic_data[rownames(exp_vector),,drop=FALSE]
    scores[cell,"prognostic_score"]<-cor(zscores[,1],exp_vector[,1], method = cor_method)
    
    count<-count+1
    if(count%%1000==0){message(paste("Calculated Scores For:",as.character(count), "Cells"))}
  }
  
  return(scores)
}

# modified scoring function - fold change between adverse and favorable genes
run_prognostic_scoring_lfc<-function(expression, prognostic_data, z_cutoff=2, normalize_gene_values=FALSE, cor_method="pearson"){
  
  # subset gene data to genes that are above z-score threshold.
  message("Finding Significant Genes...")
  colnames(prognostic_data)<-c("prognostic_value")
  prognostic_data<-prognostic_data[abs(prognostic_data[,"prognostic_value"])>z_cutoff,,drop=FALSE]
  
  positive_genes<-rownames(prognostic_data)[prognostic_data[,"prognostic_value"]>z_cutoff]
  negative_genes<-rownames(prognostic_data)[prognostic_data[,"prognostic_value"]<(z_cutoff*-1)]
  
  positive_genes<-intersect(positive_genes,rownames(expression))
  negative_genes<-intersect(negative_genes,rownames(expression))
  
  # if log2 unlog
  if(max(expression)<30){
    expression<-(2^expression)-1
  }
  
  data.1 <- apply(X = expression[positive_genes,], MARGIN = 2, FUN = function(x) log2(x = mean(x) + 1)) 
  data.2 <- apply(X = expression[negative_genes,], MARGIN = 2, FUN = function(x) log2(x = mean(x) + 1))

  total.diff <- (data.1 - data.2)
  
  scores<-data.frame(prognostic_score=total.diff)
  rownames(scores)<-colnames(expression)
  
  return(scores)
}
