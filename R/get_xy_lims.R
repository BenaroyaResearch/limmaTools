#' Accessory function for determining x and y limits for volcano plots
#'
#' This function determines appropriate x and y limits for volcano plots or sets of volcano plots.
#' @param topGenes a data frame or list of data frame. Each data frame typically contains the output of a call to \code{topTable}. Must contain genes, log2 fold-change, and adjusted p-values.
#' @param pairwise logical, whether \code{topGenes} contains multiple comparisons. Defaults to FALSE. If TRUE, the function finds common limits for all plots based on the most extreme values across all comparisons.
#' @param x_symmetric logical, whether the x-axis should be made symmetric.
#' @param adjust_factor numeric, the factor for adjusting the limits of the x- and y-axes. The returned limits will be extended by this much to allow points at the extremes of the plot to fit. Defaults to 0.01.
#' @param y_floor_zero logical, whether to set the lower limit of y to 0. Defaults to TRUE.
#' @param x_col,y_col character or integer, the name or number of the columns containing the data to be plotted on the x- and y-axes. Default to "logFC" and "adj.P.Val", respectively.
#' @param neg_log10_ycol logical, whether the y-axis will be plotted as \code{(-log10(y))}. Defaults to TRUE, as the function is intended for use with p-values.
#' @param min_x_abs numeric, the minimum absolute value for the extent of the x-axis. Defaults to NULL, which sets no default. Typically used when plotting with thresholds, to ensure that threshold value is shown.
#' @param min_y2 numeric, the minimum value for upper limit the y-axis. Defaults to NULL, which sets no default. Typically used when plotting with thresholds, to ensure that threshold value is shown.
#' @export
#' @return A list, with vectors for x- and y-limits.
#' @details This function finds reasonable x- and y-limits for plots; it is specifically designed for use with volcano plots.
#' @usage \code{
#' get_xy_lims(topGenes,
#'             pairwise=FALSE,
#'             x_symmetric=TRUE,
#'             adjust_factor=0.01,
#'             y_floor_zero=TRUE,
#'             x_col="logFC", y_col="adj.P.Val",
#'             neg_log10_ycol=TRUE,
#'             min_x_abs=NULL,
#'             min_y2=NULL)}
get_xy_lims <- function(topGenes,
                        pairwise=FALSE,
                        x_symmetric=TRUE,
                        adjust_factor=0.01,
                        y_floor_zero=TRUE,
                        x_col="logFC", y_col="adj.P.Val",
                        neg_log10_ycol=TRUE,
                        min_x_abs=NULL,
                        min_y2=NULL) {
  if (pairwise) {
    x_lims <- range(unlist(lapply(topGenes, function(x) x[,grep(x_col, colnames(x), value=TRUE)])))
    y_lims <- range(unlist(lapply(topGenes, function(x) x[,grep(y_col, colnames(x), value=TRUE)])))
  } else {
    x_lims <- range(topGenes[,x_col])
    y_lims <- range(topGenes[,y_col])
  }
  
  if (neg_log10_ycol) y_lims <- rev(-log10(y_lims))
  if (!is.null(min_y2) & (y_lims[2] < min_y2)) y_lims[2] <- min_y2
  
  if (!is.null(min_x_abs)) {
    if (x_lims[1] > -min_x_abs) x_lims[1] <- min_x_abs
    if (x_lims[2] < min_x_abs) x_lims[2] <- min_x_abs
  }
  if (x_symmetric) x_lims <- c(-1,1) * max(abs(x_lims))
  x_adjustment <- (x_lims[2] - x_lims[1]) * adjust_factor
  x_lims <- x_lims + (c(-1,1) * x_adjustment)
  
  if (y_floor_zero & (y_lims[1] >= 0)) {
    y_lims[1] <- 0
    y_adjustment <- (y_lims[2] - y_lims[1]) * adjust_factor
    y_lims[2] <- y_lims[2] + y_adjustment
  } else {
    y_adjustment <- (y_lims[2] - y_lims[1]) * adjust_factor
    y_lims <- y_lims + (c(-1,1) * y_adjustment)
  }
  
  list(x=x_lims, y=y_lims)  
}
