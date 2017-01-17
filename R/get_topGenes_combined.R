#' Construction of combined topGenes object from a set of pairwise contrasts in limma models
#'
#' This function takes a set of pairwise limma model contrasts, and combines the results into a single
#' object.
#' @param topGenes a data frame, typically the output of a call to \code{topTable}. Must contain genes, log2 fold-change, and adjusted p-values.
#' @param topGenes.pairwise a list of data frames, each typically containing the output of a call to \code{topTable} for a single contrast. Each list element should be named with an identifier for the contrast, and must contain genes, log2 fold-change, and adjusted p-values. Can optionally include a "threshold" column, which should be logical indicating genes passing significance thresholds.
#' @param p_cut numeric, the cutoff for adjusted p-value. Genes with adjusted p-values greater than or equal to this value are not considered significant.
#' @param fc_cut numeric, the absolute value cutoff for log2 fold change. Genes with absolute value log2-FC less than or equal to this value are not considered significant.
#' @export
#' @return A data frame, with adjusted p-values and logical \code{threshold} values for all contrasts in \code{topGenes.pairwise}
#' @usage \code{
#' get_topGenes_combined(topGenes,
#'                       topGenes.pairwise,
#'                       p_cut=0.01, fc_cut=log2(1.5))}
get_topGenes_combined <-
  function(topGenes,
           topGenes.pairwise,
           p_cut=0.01, fc_cut=log2(1.5)) {
    topGenes.combined <-
      topGenes[topGenes$adj.P.Val < p_cut, c("F","P.Value", "adj.P.Val"), drop=FALSE]  # extract genes with significant F test, retaining average expression, F-stat (and its associated p-value and q-value)
    for (i in names(topGenes.pairwise)) {
      colnames(topGenes.pairwise[[i]]) <- 
        paste(colnames(topGenes.pairwise[[i]]), i, sep=".")
      topGenes.combined <-
        merge(topGenes.combined, topGenes.pairwise[[i]][,c(1,3,4), drop=FALSE],
              by="row.names", all.x=TRUE, sort=FALSE)
      rownames(topGenes.combined) <- topGenes.combined$Row.names
      topGenes.combined$Row.names <- NULL
    }
    adj.p.values <-
      matrix(
        data=p.adjust(as.matrix(topGenes.combined[,str_detect(names(topGenes.combined), "P.Value.")]),
                     method="BH"),
        ncol=length(topGenes.pairwise),
        dimnames=list(rownames(topGenes.combined),
                      paste("adj.P.Val", names(topGenes.pairwise), "combined", sep=".")))
    topGenes.combined <- cbind(topGenes.combined, adj.p.values)
    for (i in names(topGenes.pairwise)) {
      topGenes.combined[[paste0("threshold.", i)]] <-
        (abs(topGenes.combined[[paste0("logFC.", i)]]) > fc_cut) &
        (topGenes.combined[[paste0("adj.P.Val.", i, ".combined")]] < p_cut)
    }
    topGenes.combined$threshold.all <-
      apply(topGenes.combined[,grep("threshold", colnames(topGenes.combined), value=TRUE)],
          MARGIN=1, FUN=any)
  topGenes.combined
}
