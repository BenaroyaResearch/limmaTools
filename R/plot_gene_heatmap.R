#' Plot a heatmap of gene expression counts
#'
#' This is a wrapper for \code{heatmap3}, designed to make plotting of gene expression heatmaps more
#' automated. Samples can be order and/or colored based on variables in an optional \code{design} object.
#' It also provides a nicer default color palette than the heatmap3 default.
#' 
#' @param counts a numeric matrix (or object that can be coerced to matrix) of gene expression values, with genes in rows and samples in columns, or an object from which counts can be extracted (such as an EList or DGEList).
#' @param design (optional) design object, used for ordering or coloring samples. If provided, must include a column (named in \code{libID_col}) with identifiers matching column names in \code{counts}, and columns matching \code{color_by_var} and/or \code{order_by_var}.
#' @param libID_col (optional) string, the name of the column in \code{design} containing sample identifiers matching the column names of \code{counts}.
#' @param order_by_vars (optional) character vector, the names of the columns in \code{design} by which to order the samples. If multiple columns are provided, ordering is done by the first variable, then the second, etc. by passing the columns as argument to \code{dplyr::arrange}
#' @param order_by_vars_levels (optional) list of character vectors, one for each element in \code{order_by_vars}, with the levels of \code{order_by_vars} in the order to use for plotting. If not provided, elements are ordered by factor levels (if a factor) or based on first appearance in the design object. Ignored if \code{order_by_vars} is numeric.
#' @param color_by_vars (optional) character vector, the names of the columns in \code{design} by which to color the sample labels. If \code{order_by_vars} is specified and \code{color_by_vars} is not, sample labels will be colored by \code{order_by_vars}.
#' @param my_vars_colors (optional) list of vectors of colors for use in coloring column identifiers, one for each element in \code{color_by_vars}. Use varies depending on class of \code{color_by_vars}. For each variable specified by \code{color_by_vars}: if the variable is numeric, values of are broken up into even intervals with number specified by \code{color_byvars_levels}, and assigned colors from the corresponding vector of colors in \code{my_vars_colors}; if needed, additional colors are imputed using \code{colorRampPalette}. If not provided, the "YlGnBu" palette from \code{RColorBrewer} is used. If the variable is not numeric, \code{my_vars_colors} should include at least as many colors as there are unique values for the variable; these colors are assigned to the values for plotting. If not provided, "Set1" from \code{RColorBrewer} is used.
#' @param color_by_vars_levels (optional) list of vectors, one for each element in \code{color_by_vars}, providing control over coloring of sample labels. Use varies depending on the class of the variable. If the variable is NOT numeric, should include the unique values of \code{color_by_vars}; this vector is matched to \code{my_vars_colors} to allow control of color labels. If an empty list is provided, elements are matched to colors in order by factor levels (if a factor) or based on first appearance in the design object. If the variable is numeric, should be an integer providing the number of intervals to split the variable into (defaults to 10).
#' @param norm.method name of the function to be used in normalizing the gene expression values in each row. Defaults to "range01", which normalizes each row to extend from 0 to 1. Can be any function that returns a numeric vector of the same length as its argument. Passed to \code{match.fun}. To use counts without normalization, use NULL or "identity".
#' @param scale alternative method for normalizing counts. Passed to \code{heatmap3}. Defaults to "none", to allow control via norm.method. Can be "row" or "column", specifying centering and scaling over rows or columns. Use in combination with \code{norm.method} may yield unexpected results.
#' @param row_dendro,col_dendro variables that specify the row and/or columns dendrogram(s). Set to NA to suppress dendrograms. if \code{order_by_vars} is specified, \code{col_dendro} is ignored.
#' @param filename a character string. If provided, the function outputs a pdf of the plot, named "{filename}.pdf". If not provided, the function prints to a plotting window.
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @param my_heatmap_cols a vector of color names, typically the result of a call to \code{colorRampPalette}.
#' @param add_legend logical, whether to add a legend for the coloring of sample labels. Defaults to \code{FALSE}.
#' @param leg_x,leg_y x- and y-coordinates for the plot location of the legend for non-numeric variables. As the coordinates of heatmap3 vary with the data set, the defaults may not provide a good location.
#' @param xl,xr,yb,yt x- and y-coordinates for the left, right, bottom, and top boundaries of the legend for numeric variables. As the coordinates of heatmap3 vary with the data set, the defaults may not provide a good location.
#' @param ... (optional) additional arguments passed to \code{heatmap3}.
#' @export
#' @usage \code{
#' plot_gene_heatmap(
#'      counts, design=NULL, libID_col="lib.id",
#'      order_by_vars=NULL, order_by_vars_levels=NULL,
#'      color_by_vars=order_by_vars, my_vars_colors=NULL, color_by_vars_levels=NULL,
#'      norm.method="range01", scale="none", row_dendro=NULL, col_dendro=NULL,
#'      filename=NULL, plotdims=c(9,9),
#'      my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
#'      add_legend=FALSE, leg_x=0.7, leg_y=1.1, xl=0.6, xr=0.7, yb=0.9, yt=1.1,
#'      ...)}
plot_gene_heatmap <-
  function(counts, design=NULL, libID_col="lib.id",
           order_by_vars=NULL, order_by_vars_levels=NULL,
           color_by_vars=order_by_vars, my_vars_colors=NULL, color_by_vars_levels=NULL,
           norm.method="range01", scale="none", row_dendro=NULL, col_dendro=NULL,
           filename=NULL, plotdims=c(9,9),
           my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           add_legend=FALSE, leg_x=0.7, leg_y=1.1, xl=0.6, xr=0.7, yb=0.9, yt=1.1,
           ...) {
    if (!is.null(order_by_vars) & is.null(design)) stop("Cannot sort the libraries without annotation data.")
    
    counts <- extract_counts(counts)
    
    if (any(dim(counts)<2)) stop("Counts object must include at least two genes and two libraries.")
    
    # sort counts by specified variable, if applicable
    if (!is.null(order_by_vars)) {
      for (i in order_by_vars) {
        if (!is.numeric(design[,order_by_vars[i]])) {
          if (is.null(order_by_vars_levels[[i]])) {
            if (is.factor(design[,order_by_vars[i]])) {
              order_by_vars_levels[[i]] <-
                levels(design[,order_by_vars[i]])
            } else {
              order_by_vars_levels[[i]] <-
                as.character(unique(design[,order_by_vars[i]]))
            }
          }
          design[,order_by_vars] <-
            factor(design[,order_by_vars], levels=order_by_vars_levels) # set order
        }
      }
      
      counts <- 
        counts[
          , match(
            design[do.call(order, unname(design[, order_by_vars, drop=FALSE])), libID_col],
            colnames(counts))] # put columns in new order
      col_dendro <- NA
    }
    
    # now mirror the updated order of counts in design
    if (!is.null(design))
      design <- design[match(colnames(counts), design[,libID_col]),]
    
    if (is.null(norm.method))
      norm.method <- "identity"
    counts <-
      t(apply(counts, MARGIN=1, FUN=match.fun(norm.method, descend=FALSE))) # normalize counts
    
    if (!is.null(color_by_vars)) {
      # make a data frame of colors for labeling columns
      plot_colors <- data.frame()
      
      # need to modify this to generate a data frame, where each column has the colors to pass to heatmap3 as ColSideColors, and is named with the variable name
      for (i in color_by_vars) {
        color_by_vars <- design[,color_by_vars]
        
        if (is.numeric(color_by_vars)) {
          
          # I should update this to use my values2colors function
          
          if (is.null(color_by_vars_levels)) color_by_vars_levels <- 10
          if (is.null(my_vars_colors)) {
            require(RColorBrewer)
            my_vars_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(color_by_vars_levels)
          }
          
          if (length(my_vars_colors) != color_by_vars_levels)
            my_vars_colors <- colorRampPalette(my_vars_colors)(color_by_vars_levels)
          
          plot_colors <-
            my_vars_colors[
              findInterval(
                color_by_vars, 
                vec=seq(min(color_by_vars), max(color_by_vars),
                        length.out=color_by_vars_levels+1),
                rightmost.closed=TRUE)]
        } else {
          if (is.null(my_vars_colors)) {
            require(RColorBrewer)
            my_vars_colors <- brewer.pal(length(unique(color_by_vars)), "Set1")
          }
          if (is.null(color_by_vars_levels)) {
            if (is.factor(color_by_vars)) {
              color_by_vars_levels <- levels(color_by_vars)
            } else {
              color_by_vars_levels <- as.character(unique(color_by_vars))
            }
            color_by_vars <- factor(color_by_vars, levels=color_by_vars_levels)
          }
          
          my_vars_colors <- my_vars_colors[1:length(color_by_vars_levels)]
          names(my_vars_colors) <- color_by_vars_levels
          plot_colors <- my_vars_colors[color_by_vars]
        }
      }
    }
    
    # open plotting device
    if (!is.null(filename)) {
      pdf(file=filename, w=plotdims[1], h=plotdims[2])
      on.exit(dev.off()) # close plotting device on exit
    } else quartz(plotdims[1],plotdims[2])
    
    # generate heatmap
    if (exists("plot_colors")) {
      heatmap3::heatmap3(
        counts, scale=scale, col=my_heatmap_cols,
        Rowv=row_dendro, Colv=col_dendro,
        margins=c(8,12),
        ColSideColors=plot_colors,
        ...)
    } else
      heatmap3::heatmap3(
        counts, scale=scale, col=my_heatmap_cols,
        Rowv=row_dendro, Colv=col_dendro,
        margins=c(8,12),
        ...)
    
    if (add_legend)
      if (is.numeric(color_by_vars)) {
        legend_color_labels <- function(x) {
          x <- as.character(x)
          y <- ""
          for (i in 2:(length(x)-1)) {
            y[(i-1)*2] <- x[i]
            y[(i-0.5)*2] <- ""
          }
          return(y)
        }
        
        plotrix::color.legend(
          legend=
            legend_color_labels(
              round(seq(min(color_by_vars),
                        max(color_by_vars),
                        length.out=color_by_vars_levels+1),
                    round(-log10(diff(range(color_by_vars)))+2))),
          rect.col=my_vars_colors,
          gradient="y",
          xl=0.6, xr=0.7, yb=0.9, yt=1.1,)
      } else {
        legend(legend=names(my_vars_colors), fill=my_vars_colors,
               bty="n", x=leg_x, y=leg_y, xpd=TRUE)
      }
    }
