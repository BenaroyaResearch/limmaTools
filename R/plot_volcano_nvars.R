#' Generate a set of volcano plots of genes from a differential expression (limma) analysis, for one or more contrasts
#'
#' Generate a set of volcano plots of genes from a differential expression (limma) analysis, for one or
#' more contrasts. These plots can be output to plotting windows, or to pdfs. The points can be colored
#' based on fold-change and p-value thresholds. Points can also be labeled with gene names, and the points
#' to be labeled can be set based on an ellipse oriented to the x- and y-axes.
#' @param topGenes.pairwise a list of data frames, each typically containing the output of a call to \code{topTable} for a single contrast. Each list element should be named with an identifier for the contrast, and must contain genes, log2 fold-change, and adjusted p-values. Can optionally include a "threshold" column, which should be boolean indicating genes passing significance thresholds.
#' @param my_cols a vector of colors for plotting points. First element provides the color for points not exceeding significance thresholds; second element provides the color for points exceeding significance thresholds.
#' @param file_prefix a character string. If provided, the function outputs pdfs of the plots, named "{file_prefix}.{list_element_name}.pdf".
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @param fc_cut numeric, the (absolute value) log2 fold-change threshold for determining significance of genes. This value is also plotted as vertical dotted lines.
#' @param p_cut numeric, the p-value threshold for determining significance of genes. This value is also plotted as a horizontal dotted line.
#' @param x_lim,y_lim either "auto", NULL, or numeric vectors. If "auto", x- and y-limits are determined from the data using \code{get_xy_lims}. If NULL, default plot limits are used. If provided as numeric vectors, the lower and upper limits of the plotting space along the x- and y-axes. Passed to \code{ggplot2::xlim}.
#' @param gene_labs logical, whether to include gene labels for genes with extreme logFC and p-value. If \code{TRUE}, genes with values outside the labeling ellipse will be labeled.
#' @param x_cut,y_cut numeric, the radii of the labeling ellipse along the x- and y-axes. Genes with values outside the ellipse are labeled with gene names. Default to 0, which results in all genes being labeled.
#' @param ... additional parameters passed to \code{pdf}.
#' @import ggplot2
#' @details A separate plot is generated for each element of topGenes.pairwise. 
#' @export
#' @usage \code{
#' plot_volcano_nvars(topGenes.pairwise, my_cols=c("darkcyan", "darkorange"),
#'                    file_prefix=NULL, plotdims=c(9,9),
#'                    fc_cut=log2(1.5), p_cut=0.01,
#'                    x_lim="auto", y_lim="auto",
#'                    gene_labs=FALSE, x_cut=0, y_cut=0,
#'                    ...)}
plot_volcano_nvars <- function(topGenes.pairwise, my_cols=c("darkcyan", "darkorange"),
                               file_prefix=NULL, plotdims=c(9,9),
                               fc_cut=log2(1.5), p_cut=0.01,
                               x_lim="auto", y_lim="auto",
                               gene_labs=FALSE, x_cut=0, y_cut=0,
                               ...) {
  if (identical(x_lim, "auto") | identical(y_lim, "auto")) {
    xy_lims <- get_xy_lims(topGenes.pairwise, pairwise=TRUE, min_x_abs=fc_cut, min_y2=-log10(p_cut))
    if (identical(x_lim, "auto")) x_lim <- xy_lims[["x"]]
    if (identical(y_lim, "auto")) y_lim <- xy_lims[["y"]]
  }
  
  for (i in names(topGenes.pairwise)) {
    topGenes.tmp <- topGenes.pairwise[[i]]
    colnames(topGenes.tmp) <- str_replace(colnames(topGenes.tmp), paste0(".", i), "") # strip out the de-ambiguation from the column names, if present
    topGenes.tmp$genes <- rownames(topGenes.tmp)
    if (is.null(topGenes.tmp$threshold)) topGenes.tmp$threshold <- (abs(topGenes.tmp$logFC) > fc_cut) & (topGenes.tmp$adj.P.Val < p_cut)
    
    volcano <- ggplot(data = topGenes.tmp,
                      aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
      geom_point(alpha=0.6, size=3) +
      theme(legend.position = "none") +
      xlab("log2 fold change") + ylab("-log10 Adj P") +
      geom_vline(xintercept = fc_cut, linetype="dotted",size=1.0) +
      geom_vline(xintercept = -fc_cut, linetype="dotted",size=1.0) +
      geom_hline(yintercept = -log10(p_cut), linetype="dotted",size=1.0) + 
      scale_colour_manual(values=my_cols)
    if (!is.null(x_lim)) {volcano <- volcano + xlim(x_lim)}
    if (!is.null(y_lim)) {volcano <- volcano + ylim(y_lim)}
    if (gene_labs) {
      volcano <- volcano +
        geom_text(data=topGenes.tmp[((topGenes.tmp$logFC^2)/(x_cut^2) + (log10(topGenes.tmp$adj.P.Val)^2)/(y_cut^2)) > 1,],
                  aes(label=genes),
                  color="black", size=3, vjust=1, hjust=0.5)}
    
    if (!is.null(file_prefix)) {
      pdf(file=paste(file_prefix, i, "pdf", sep="."), w=plotdims[1], h=plotdims[2], ...)
    } else {
      quartz(plotdims[1],plotdims[2])
      volcano <- volcano + ggtitle(i)
    }
    
    print(volcano)
    
    if (!is.null(file_prefix)) dev.off()
  }
}
