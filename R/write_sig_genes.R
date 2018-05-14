#' Output text files with lists of significant genes
#'
#' This function outputs text files with the results from a differential expression analysis. It provides
#' options to output based on FDR and logFC thresholds, positive and negative logFC, and ranked lists.
#' @param topGenes a data frame, typically the output of a call to \code{topTable}. Must contain genes, log2 fold-change, and adjusted p-values.
#' @param file_prefix name of the destination for files. Details of each output will be appended to this prefix.
#' @param method character, specifying the type of gene lists to output. "ranked_list" outputs a list of all genes in ranked order by p-value (from smallest to largest). "combined" outputs a list of all significant genes meeting the threshold. "directional" outputs lists of significant genes meeting the threshold that are up- and down-regulated. Partial matches are allowed.
#' @param adj_p_cut numeric, the cutoff for adjusted p-value. Genes with adjusted p-values greater than or equal to this value are not included in the result. Defaults to 0.01.
#' @param fc_cut numeric, the absolute value cutoff for log2 fold change. Genes with absolute value log2-FC less than or equal to this value are not included in the result. Defaults to log2(1.5). To include all genes, set to 0. To ignore logFC, set to NULL.
#' @param fc_adj_factor numeric, the adjustment factor used for log2-fold-change values with a numeric predictor. This is included so that the output file names can include the fold change prior to scaling. Defaults to 1, which is the appropriate value for categorical comparisons.
#' @param p_col name or number of the column in \code{topGenes} on which to sort. Generally the raw p-values, as adjusted p-values are often homogenized across a range of raw p-values. Defaults to "P.Value", which corresponds to the output from \code{topTable}. To include all genes, set to >1.
#' @param adj_p_col name or number of the column in \code{topGenes} containing the p-values to compare to \code{p_cut}. Defaults to "adj.P.Val", which corresponds to the output from \code{topTable}.
#' @param fc_col name or number of the column in \code{topGenes} containing the fold-change values to compare to \code{fc_cut}. Defaults to "logFC", which corresponds to the output from \code{topTable}.
#' @param threshold_col name or number of the column in \code{topGenes} containing the logical values indicating which genes meet thresholds. This is an alternate way to determine signficance of genes. If specified, \code{p_cut} and \code{fc_cut} are ignored.
#' @param input_type the input gene identifier class. Must match a variable type in \code{annotables} or \code{biomaRt}. Defaults to "symbol" with \code{annotables} or "hgnc_symbol" with \code{biomaRt}.
#' @param output_type the output gene identifier class. Must match a variable type in \code{annotables} or \code{biomaRt}. Defaults to
#' @param use_annotables logical, whether to use the annotables package to convert output_type, if necessary. Only used if output_type is specified. If annotables is not installed, the function defaults to using biomaRt.
#' @export
#' @details This function writes out lists of genes to text files. By default, it outputs a list ranked by p-value, lists of genes significant based on FDR and logFC thresholds (all, up, and down).
#' @usage \code{
#' write_sig_genes(
#'   topGenes, file_prefix,
#'   method=c("ranked_list", "combined", "directional"),
#'   adj_p_cut=0.01, fc_cut=log2(1.5), fc_adj_factor=1,
#'   p_col="P.Value", adj_p_col="adj.P.Val", fc_col="logFC",
#'   threshold_col=NULL,
#'   input_type, output_type=input_type, use_annotables=TRUE)}
write_sig_genes <-
  function(topGenes, file_prefix,
           method=c("ranked_list", "combined", "directional"),
           adj_p_cut=0.01, fc_cut=log2(1.5), fc_adj_factor=1,
           p_col="P.Value", adj_p_col="adj.P.Val", fc_col="logFC",
           threshold_col=NULL,
           input_type, output_type=input_type, use_annotables=TRUE) {
    if (!is.data.frame(topGenes)) stop("topGenes must be a data frame object")
    if (!missing(input_type))
      rownames(topGenes) <-
        convert_gene_names(rownames(topGenes), input_type, output_type, use_annotables=use_annotables)
             
    method <- match.arg(method, c("ranked_list", "combined", "directional"), several.ok=TRUE)
    
    for (method.tmp in method) { # iterate over selected methods
      genes_to_output <-
        get_sig_genes(topGenes, method=method.tmp, adj_p_cut=adj_p_cut, fc_cut=fc_cut,
                      p_col=p_col, adj_p_col=adj_p_col, fc_col=fc_col,
                      threshold_col=threshold_col)
      
      
      if (method.tmp == "ranked_list") { # output ranked list
        write.table(genes_to_output, file=paste0(file_prefix, ".all_genes_ranked_pval.txt"),
                    quote = FALSE, col.names=FALSE, row.names=FALSE)
      } else if (method.tmp %in% c("combined", "directional")) {
        
        # update filename
        if (!is.null(threshold_col)) { 
          threshold_text <- "_threshold"
        } else if (is.null(fc_cut)) {
          threshold_text <- paste0("_P", adj_p_cut)
        } else if ((fc_cut <= 0) & (adj_p_cut <= 1)) {
          threshold_text <- paste0("_P", adj_p_cut)
        } else if ((fc_cut > 0) & (adj_p_cut <= 1)) {
          threshold_text <- paste0("_FC", round(2^(fc_cut*fc_adj_factor), 3), "_and_P", adj_p_cut)
        } else if (fc_cut > 0) {
          threshold_text <- paste0("_FC", round(2^(fc_cut*fc_adj_factor), 3))
        } else {
          threshold_text <- ""
        }
        
        # output combined list of significant genes
        if (method.tmp == "combined") { 
          write.table(
            genes_to_output,
            file=paste0(file_prefix, ".genes", threshold_text, ".txt"),
            quote = FALSE, col.names=FALSE, row.names=FALSE)
        }
        
        # output directional lists of significant genes
        if (method.tmp == "directional") {
          write.table(
            genes_to_output[["up"]],
            file=paste0(file_prefix, ".genes", threshold_text, ".up.txt"),
            quote = FALSE, col.names=FALSE, row.names=FALSE)
          write.table(
            genes_to_output[["down"]],
            file=paste0(file_prefix, ".genes", threshold_text, ".down.txt"),
            quote = FALSE, col.names=FALSE, row.names=FALSE)
        }
      }
    }
  }
