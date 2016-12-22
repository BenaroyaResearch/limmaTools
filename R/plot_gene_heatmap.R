#' Plot a heatmap of gene expression counts
#'
#' This is a wrapper for \code{heatmap.2}, designed to make plotting of gene expression heatmaps more
#' automated. Samples can be order and/or colored based on variables in an optional \code{design} object.
#' It also provides a nicer default color palette than the heatmap.2 default.
#' 
#' @param counts a numeric matrix (or object that can be coerced to matrix) of gene expression values, with genes in rows and samples in columns.
#' @param design (optional) design object, used for ordering or coloring samples. If provided, must include a column (named in \code{libID_col}) with identifiers matching column names in \code{counts}, and columns matching \code{color_by_var} and/or \code{order_by_var}.
#' @param libID_col (optional) string, the name of the column in \code{design} containing sample identifiers matching the column names of \code{counts}.
#' @param order_by_var (optional) string, the name of the column in \code{design} by which to order the samples.
#' @param order_by_var_levels (optional) character vector, with the levels of \code{order_by_var} in the order to use for plotting. If not provided, elements are ordered based on first appearance in the design object. Ignored if \code{order_by_var} is numeric.
#' @param color_by_var (optional) string, the name of the column in \code{design} by which to color the sample labels. If \code{order_by_var} is specified and \code{color_by_var} is not, sample labels will be colored by \code{order_by_var}.
#' @param my_var_colors (optional) vector of colors for use in coloring column identifiers. Use varies depending on class of \code{color_by_var}. If \code{color_by_var} is numeric, values of \code{color_by_var} are broken up into even intervals with number specified by \code{color_by_var_levels}, and assigned colors from \code{my_var_colors}; if needed, additional colors are imputed using \code{colorRampPalette}. If not provided, the "YlGnBu" palette from \code{RColorBrewer} is used. If \code{color_by_var} is not numeric, \code{my_var_colors} should include at least as many colors as there are unique values in \code{color_by_var}; these colors are assigned to the values for plotting. If not provided, "Set1" from \code{RColorBrewer} is used.
#' @param color_by_var_levels (optional) character vector or numeric value, providing control over coloring of sample labels. Use varies depending on class of \code{color_by_var}. If \code{color_by_var} is NOT numeric, should include the unique values of \code{color_by_var}; this vector is matched to \code{my_var_colors} to allow control of color labels. If not provided, elements are colored based on first appearance in the design object. If \code{color_by_var} is numeric, should be an integer providing the number of intervals to split \code{color_by_var} into (defaults to 10).
#' @param norm.method name of the function to be used in normalizing the gene expression values in each row. Defaults to "range01", which normalizes each row to extend from 0 to 1. Can be any function that returns a numeric vector of the same length as its argument. Passed to \code{match.fun}. To use counts without normalization, use NULL or "identity".
#' @param scale alternative method for normalizing counts. Passed to \code{heatmap.2}. Defaults to "none", to allow control via norm.method. Can be "row" or "column", specifying centering and scaling over rows or columns. Use in combination with \code{norm.method} may yield unexpected results.
#' @param row_dendro,col_dendro boolean, whether to include the row and/or columns dendrogram(s)
#' @param filename a character string. If provided, the function outputs a pdf of the plot, named "{filename}.pdf". If not provided, the function prints to a plotting window.
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @param my_heatmap_cols a vector of color names, typically the result of a call to \code{colorRampPalette}.
#' @param add_legend boolean, whether to add a legend for the coloring of sample labels. Defaults to \code{FALSE}.
#' @param leg_x,leg_y x- and y-coordinates for the plot location of the legend for non-numeric variables. As the coordinates of heatmap.2 vary with the data set, the defaults may not provide a good location.
#' @param xl,xr,yb,yt x- and y-coordinates for the left, right, bottom, and top boundaries of the legend for numeric variables. As the coordinates of heatmap.2 vary with the data set, the defaults may not provide a good location.
#' @param key logical, whether to add a color key to the heatmap. Passed to \code{heatmap.2}. Defaults to \code{FALSE}.
#' @param ... (optional) additional arguments passed to \code{heatmap.2}.
#' @export
#' @usage \code{
#' plot_gene_heatmap(counts, design=NULL, libID_col="lib.id",
#'                   order_by_var=NULL, order_by_var_levels=NULL,
#'                   color_by_var=order_by_var, my_var_colors=NULL, color_by_var_levels=NULL,
#'                   norm.method="range01", scale="none", row_dendro=TRUE, col_dendro=FALSE,
#'                   filename=NULL, plotdims=c(9,9),
#'                   my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
#'                   add_legend=FALSE, leg_x=0.7, leg_y=1.1, xl=0.6, xr=0.7, yb=0.9, yt=1.1,
#'                   key=FALSE,
#'                   ...)}
plot_gene_heatmap <-
  function(counts, design=NULL, libID_col="lib.id",
           order_by_var=NULL, order_by_var_levels=NULL,
           color_by_var=order_by_var, my_var_colors=NULL, color_by_var_levels=NULL,
           norm.method="range01", scale="none", row_dendro=TRUE, col_dendro=FALSE,
           filename=NULL, plotdims=c(9,9),
           my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           add_legend=FALSE, leg_x=0.7, leg_y=1.1, xl=0.6, xr=0.7, yb=0.9, yt=1.1,
           key=FALSE,
           ...) {
    if (!is.null(order_by_var) & is.null(design)) stop("Cannot sort the libraries without annotation data.")
    
    # sort counts by specified variable, if applicable
    if (!is.null(order_by_var)) {
      if (!is.numeric(design[,order_by_var])) {
        if (is.null(order_by_var_levels))
          order_by_var_levels <- as.character(unique(design[,order_by_var]))
        design[,order_by_var] <- factor(design[,order_by_var], levels=order_by_var_levels) # set order
        }
      counts <- counts[,order(design[, order_by_var])] # put columns in new order
    }
    
    # put design in order of counts
    if (!is.null(design))
      design <- design[match(colnames(counts), design[,libID_col]),]
    
    if (is.null(norm.method))
      norm.method <- "identity"
    counts <- t(apply(counts, MARGIN=1, FUN=match.fun(norm.method, descend=FALSE))) # normalize counts
    
    if (!is.null(color_by_var)) {
      color_by_var <- design[,color_by_var]
      
      if (is.numeric(color_by_var)) {
        
        # I should update this to use my values2colors function
        
        if (is.null(color_by_var_levels)) color_by_var_levels <- 10
        if (is.null(my_var_colors)) {
          require(RColorBrewer)
          my_var_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(color_by_var_levels)
        }
        
        if (length(my_var_colors) != color_by_var_levels) my_var_colors <- colorRampPalette(my_var_colors)(color_by_var_levels)
        
        plot_colors <-
          my_var_colors[findInterval(color_by_var, 
                                     vec=seq(min(color_by_var), max(color_by_var),
                                             length.out=color_by_var_levels+1),
                                     rightmost.closed=TRUE)]
      } else {
        if (is.null(my_var_colors)) {
          require(RColorBrewer)
          my_var_colors <- brewer.pal(length(unique(color_by_var)), "Set1")
        }
        if (is.null(color_by_var_levels)) color_by_var_levels <- as.character(unique(color_by_var))
        color_by_var <- factor(color_by_var, levels=color_by_var_levels)
        
        my_var_colors <- my_var_colors[1:length(color_by_var_levels)]
        names(my_var_colors) <- color_by_var_levels
        plot_colors <- my_var_colors[color_by_var]
      }
    }
    
    # open plotting device
    if (!is.null(filename)) {
      pdf(file=filename, w=plotdims[1], h=plotdims[2])
    } else quartz(plotdims[1],plotdims[2])
    
    # generate heatmap
    if (exists("plot_colors")) {
      gplots::heatmap.2(counts, scale=scale, col=my_heatmap_cols,
                        Rowv=row_dendro, Colv=col_dendro,
                        trace="none", key=key, margins=c(8,12),
                        ColSideColors=plot_colors,
                        ...)
    } else
      gplots::heatmap.2(counts, scale=scale, col=my_heatmap_cols,
                        Rowv=row_dendro, Colv=col_dendro,
                        trace="none", key=key, margins=c(8,12),
                        ...)
    
    if (add_legend)
      if (is.numeric(color_by_var)) {
        legend_color_labels <- function(x) {
          x <- as.character(x)
          y <- ""
          for (i in 2:(length(x)-1)) {
            y[(i-1)*2] <- x[i]
            y[(i-0.5)*2] <- ""
          }
          return(y)
        }
        
        plotrix::color.legend(legend=
                                legend_color_labels(round(seq(min(color_by_var),
                                                              max(color_by_var),
                                                              length.out=color_by_var_levels+1),
                                                          round(-log10(diff(range(color_by_var)))+2))),
                              rect.col=my_var_colors,
                              gradient="y",
                              xl=0.6, xr=0.7, yb=0.9, yt=1.1,)
      } else {
        legend(legend=names(my_var_colors), fill=my_var_colors,
               bty="n", x=leg_x, y=leg_y, xpd=TRUE)
      }
    
    if (!is.null(filename)) dev.off() # close plotting device (if needed)
    
  }
