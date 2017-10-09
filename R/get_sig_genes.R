#' Output text files with lists of significant genes
#'
#' This function retrieves lists of gene identifiers based on the results of a differential expression
#' analysis. It provides options to retrieve based on FDR and logFC thresholds, positive and negative
#' logFC, and ranked lists.
#' @param topGenes a data frame, typically the output of a call to \code{topTable}. Must contain genes, log2 fold-change, and adjusted p-values.
#' @param method character, specifying the type of gene lists to return. Recognized values are "ranked_list", "combined", "directional", "up", "down". "ranked_list" returns a vector of all genes in ranked order by p-value (from smallest to largest). "combined" outputs a list of all significant genes meeting the threshold. "up" and "down" return lists of significant genes meeting the threshold that are up- or down-regulated, respectively. "directional" returns a list containing two vectors, one "up" and one "down", which are identical to those returned by the "up" and "down" methods. Partial matches are allowed.
#' @param adj_p_cut numeric, the cutoff for adjusted p-value. Genes with adjusted p-values greater than or equal to this value are not included in the result. Defaults to 0.01.
#' @param fc_cut numeric, the absolute value cutoff for log2 fold change. Genes with absolute value log2-FC less than or equal to this value are not included in the result. Defaults to log2(1.5). To include all genes, set to 0. To ignore logFC, set to NULL.
#' @param p_col name or number of the column in \code{topGenes} on which to sort. Generally the raw p-values, as adjusted p-values are often homogenized across a range of raw p-values. Defaults to "P.Value", which corresponds to the output from \code{topTable}. To include all genes, set to >1.
#' @param adj_p_col name or number of the column in \code{topGenes} containing the p-values to compare to \code{p_cut}. Defaults to "adj.P.Val", which corresponds to the output from \code{topTable}.
#' @param fc_col name or number of the column in \code{topGenes} containing the fold-change values to compare to \code{fc_cut}. Defaults to "logFC", which corresponds to the output from \code{topTable}.
#' @param threshold_col name or number of the column in \code{topGenes} containing the logical values indicating which genes meet thresholds. This is an alternate way to determine signficance of genes. If specified, \code{p_cut} and \code{fc_cut} are ignored.
#' @export
#' @details This function writes out lists of genes to text files. By default, it outputs a list ranked by p-value, lists of genes significant based on FDR and logFC thresholds (all, up, and down).
#' @usage \code{
#' get_sig_genes(
#'   topGenes,
#'   method=,
#'   adj_p_cut=0.01, fc_cut=log2(1.5),
#'   p_col="P.Value", adj_p_col="adj.P.Val", fc_col="logFC",
#'   threshold_col=NULL)}
get_sig_genes <-
  function(topGenes,
           method,
           adj_p_cut=0.01, fc_cut=log2(1.5), fc_adj_factor=1,
           p_col="P.Value", adj_p_col="adj.P.Val", fc_col="logFC",
           threshold_col=NULL) {
    if (!is.data.frame(topGenes)) stop("topGenes must be a data frame object")
    
    method <- match.arg(method, c("ranked_list", "combined", "directional", "up", "down"))
    
    topGenes <- topGenes[order(topGenes[,p_col]),] # order topGenes by p-value
    
    if (method == "ranked_list") { # return ranked vector
      genes.ranked <- rownames(topGenes)
      return(genes.ranked)
    } else if (method == "combined") { # return combined list of significant genes
      if (!is.null(threshold_col)) {
        genes.combined <-
          rownames(topGenes)[
            topGenes[,threshold_col]]
      } else if (is.null(fc_cut)) {
        genes.combined <-
          rownames(topGenes)[
            topGenes[,adj_p_col] < adj_p_cut]
      } else {
        genes.combined <-
          rownames(topGenes)[
            (topGenes[,adj_p_col] < adj_p_cut) & (abs(topGenes[,fc_col]) > fc_cut)]
      }
      return(genes.combined)
      
    } else if (method %in% c("directional","up","down")) { # return directional lists of significant genes
      if (is.null(fc_cut)) {
        stop(paste0("Cannot generate '", method, "' list of genes without a logFC threshold."))
      } else if (!is.null(threshold_col)) {
        genes.up <-
          rownames(topGenes)[
            topGenes[,threshold_col] &
              topGenes[,fc_col] > 0]
        genes.down <-
          rownames(topGenes)[
            topGenes[,threshold_col] &
              topGenes[,fc_col] < 0]
      } else {
        genes.up <-
          rownames(topGenes)[
            (topGenes[,adj_p_col] < adj_p_cut) &
              (abs(topGenes[,fc_col]) > fc_cut) &
              (topGenes[,fc_col] > 0)]
        genes.down <-
          rownames(topGenes)[
            (topGenes[,adj_p_col] < adj_p_cut) &
              (abs(topGenes[,fc_col]) > fc_cut) &
              (topGenes[,fc_col] < 0)]
      }
      return(
        switch(
          method,
          "directional" = list(up=genes.up, down=genes.down),
          "up" = genes.up,
          "down" = genes.down))
    }
  }
