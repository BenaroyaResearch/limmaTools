#' Extract counts of differentially expressed genes
#'
#' This function uses the results from a differential expression analysis to filter a counts
#' object to include only significant genes. Significance can be determined using p-value and log-fold-change
#' cuts specified in \code{p_cut} and \code{fc_cut}, respectively, or by including a threshold column in
#' \code{topGenes} and specifying that column using \code{threshold_col}.
#' @param counts a matrix or data frame containing the gene expression data, or an object from which counts can be exrracted (such as an EList or DGEList). Should have samples in columns and genes in rows. Rownames must contain gene names corresponding to gene names in \code{topGenes}.
#' @param topGenes a data frame, typically the output of a call to \code{topTable}. Must contain genes, log2 fold-change, and adjusted p-values.
#' @param p_cut numeric, the cutoff for adjusted p-value. Genes with adjusted p-values greater than or equal to this value are not included in the result. Defaults to 0.01. Ignored if \code{threshold_col} is specified.
#' @param fc_cut numeric, the absolute value cutoff for log2 fold change. Genes with absolute value log2-FC less than or equal to this value are not included in the result. Defaults to log2(1.5). Ignored if \code{threshold_col} is specified.
#' @param p_col name or number of the column in \code{topGenes} containing the p-values to compare to \code{p_cut}. Defaults to "adj.P.Val", which corresponds to the output from \code{topTable}.
#' @param fc_col name or number of the column in \code{topGenes} containing the fold-change values to compare to \code{fc_cut}. Defaults to "logFC", which corresponds to the output from \code{topTable}.
#' @param threshold_col name or number of the column in \code{topGenes} containing the logical values indicating which genes meet thresholds. This is an alternate way to determine signficance of genes. If specified, \code{p_cut} and \code{fc_cut} are ignored.
#' @import countSubsetNorm
#' @export
#' @return A matrix or data frame (matching the class of \code{counts}), with only rows with log-FC and adjusted p-values meeting specified thresholds.
#' @usage \code{
#' get_counts_sig_genes(
#'   counts, topGenes,
#'   p_cut=0.01, fc_cut=log2(1.5),
#'   p_col="adj.P.Val", fc_col="logFC",
#'   threshold_col=NULL)}
get_counts_sig_genes <-
  function(counts, topGenes,
           p_cut=0.01, fc_cut=log2(1.5),
           p_col="adj.P.Val", fc_col="logFC",
           threshold_col=NULL) {
    v_sig_fc <- if (!is.null(threshold_col)) {
      topGenes[topGenes[,threshold_col],]
    } else {
      topGenes[(topGenes[,p_col] < p_cut) & (abs(topGenes[,fc_col]) > fc_cut), ]
    }
    
    counts <- extract_counts(counts)
    counts_sig_fc <- counts[rownames(counts) %in% rownames(v_sig_fc),,drop=FALSE]
  }
