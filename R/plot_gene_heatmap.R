#' Plot a heatmap of gene expression counts
#'
#' This is a wrapper for \code{heatmap3}, designed to make plotting of gene expression heatmaps more
#' automated. Samples can be order and/or colored based on variables in an optional \code{design} object.
#' It also provides a nicer default color palette than the heatmap3 default.
#' 
#' @param counts a numeric matrix (or object that can be coerced to matrix) of gene expression values, with genes in rows and samples in columns, or an object from which counts can be extracted (such as an EList or DGEList).
#' @param design (optional) design object, used for ordering or coloring samples. If provided, must include a column (named in \code{libID_col}) with identifiers matching column names in \code{counts}, and columns matching \code{color_by_var} and/or \code{order_by_var}.
#' @param libID_col (optional) string, the name of the column in \code{design} containing sample identifiers matching the column names of \code{counts}.
#' @param order_by_var (optional) character vector, the names of the columns in \code{design} by which to order the samples. If multiple columns are provided, ordering is done by the first variable, then the second, etc. by passing the columns as argument to \code{dplyr::arrange}
#' @param order_by_var_levels (optional) character vector, or list of character vectors, one for each element in \code{order_by_var}, with the levels of \code{order_by_var} in the order to use for plotting. If not provided, elements are ordered by factor levels (if a factor) or based on first appearance in the design object. Ignored if \code{order_by_var} is numeric.
#' @param color_by_var (optional) character vector, the names of the columns in \code{design} by which to color the sample labels. If \code{order_by_var} is specified and \code{color_by_var} is not, sample labels will be colored by \code{order_by_var}.
#' @param my_var_colors (optional) vector, or list of vectors, containing colors for use in coloring column identifiers, one for each element in \code{color_by_var}. Use varies depending on class of \code{color_by_var}. For each variable specified by \code{color_by_var}: if the variable is numeric, values of are broken up into even intervals with number specified by \code{color_byvars_levels}, and assigned colors from the corresponding vector of colors in \code{my_var_colors}; if needed, additional colors are imputed using \code{colorRampPalette}. If not provided, the "YlGnBu" palette from \code{RColorBrewer} is used. If the variable is not numeric, \code{my_var_colors} should include at least as many colors as there are unique values for the variable; these colors are assigned to the values for plotting. If not provided, "Set1" from \code{RColorBrewer} is used.
#' @param color_by_var_levels (optional) vector, or list of vectors, one for each element in \code{color_by_var}, providing control over coloring of sample labels. Use varies depending on the class of the variable. If the variable is NOT numeric, should include the unique values of \code{color_by_var}; this vector is matched to \code{my_var_colors} to allow control of color labels. If an empty list is provided, elements are matched to colors in order by factor levels (if a factor) or based on first appearance in the design object. If the variable is numeric, should be an integer providing the number of intervals to split the variable into (defaults to 10).
#' @param norm.method name of the function to be used in normalizing the gene expression values in each row. Defaults to "range01", which normalizes each row to extend from 0 to 1. Can be any function that returns a numeric vector of the same length as its argument. Passed to \code{match.fun}. To use counts without normalization, use NULL or "identity".
#' @param scale alternative method for normalizing counts. Passed to \code{heatmap3}. Defaults to "none", to allow control via norm.method. Can be "row" or "column", specifying centering and scaling over rows or columns. Use in combination with \code{norm.method} may yield unexpected results.
#' @param row_dendro,col_dendro variables that specify the row and/or columns dendrogram(s). Set to NA to suppress dendrograms. if \code{order_by_var} is specified, \code{col_dendro} is ignored.
#' @param filename a character string. If provided, the function outputs a pdf of the plot, named "{filename}.pdf". If not provided, the function prints to a plotting window.
#' @param plotdims a numeric vector, the width and height (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @param my_heatmap_cols a vector of color names, typically the result of a call to \code{colorRampPalette}.
#' @param add_legend logical, whether to add a legend for the coloring of sample labels. If samples are colored by more than one variable, legends are output as separate plots. Defaults to \code{FALSE}.
#' @param leg_x,leg_y x- and y-coordinates for the plot location of the legend for non-numeric variables. As the coordinates of heatmap3 vary with the data set, the defaults may not provide a good location. Used only if color_by_var has length 1.
#' @param xl,xr,yb,yt x- and y-coordinates for the left, right, bottom, and top boundaries of the legend for numeric variables. As the coordinates of heatmap3 vary with the data set, the defaults may not provide a good location. Used only if color_by_var has length 1.
#' @param ... (optional) additional arguments passed to \code{heatmap3}.
#' @export
#' @usage \code{
#' plot_gene_heatmap(
#'      counts, design=NULL, libID_col="lib.id",
#'      order_by_var=NULL, order_by_var_levels=NULL,
#'      color_by_var=order_by_var, my_var_colors=NULL, color_by_var_levels=NULL,
#'      norm.method="range01", scale="none", row_dendro=NULL, col_dendro=NULL,
#'      filename=NULL, plotdims=c(9,9),
#'      my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
#'      add_legend=FALSE, leg_x=0.7, leg_y=1.1, xl=0.6, xr=0.7, yb=0.9, yt=1.1,
#'      ...)}
plot_gene_heatmap <-
  function(counts, design=NULL, libID_col="lib.id",
           order_by_var=NULL, order_by_var_levels=NULL,
           color_by_var=order_by_var, my_var_colors=NULL, color_by_var_levels=NULL,
           norm.method="range01", scale="none", row_dendro=NULL, col_dendro=NULL,
           filename=NULL, plotdims=c(9,9),
           my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           add_legend=FALSE, leg_x=0.7, leg_y=1.1, xl=0.6, xr=0.7, yb=0.9, yt=1.1,
           ...) {
    if (!is.null(order_by_var) & is.null(design)) stop("Cannot sort the libraries without annotation data.")
    
    counts <- extract_counts(counts)
    
    if (any(dim(counts)<2)) stop("Counts object must include at least two genes and two libraries.")
    
    # sort counts by specified variable, if applicable
    if (!is.null(order_by_var)) {
      
      # format order_by_var_levels to work with one or more values of order_by_var
      if (is.null(order_by_var_levels)) {
        order_by_var_levels <- replicate(n=length(order_by_var), NULL)
      } else if (!is.list(order_by_var_levels)) {
        order_by_var_levels <- list(order_by_var_levels)
      }
      
      for (i in 1:length(order_by_var)) {
        if (!is.numeric(design[,order_by_var[i]])) {
          # determine order of factor levels
          if (is.null(order_by_var_levels[[i]])) {
            if (is.factor(design[,order_by_var[i]])) {
              order_by_var_levels[[i]] <-
                levels(design[,order_by_var[i]])
            } else {
              order_by_var_levels[[i]] <-
                as.character(unique(design[,order_by_var[i]]))
            }
          }
          # set order of factor levels
          design[,order_by_var[i]] <-
            factor(design[,order_by_var[i]], levels=order_by_var_levels[[i]]) 
        }
      }
      
      counts <- 
        counts[
          , match(
            design[do.call(order, unname(design[, order_by_var, drop=FALSE])), libID_col],
            colnames(counts))] # put columns in new order
      col_dendro <- NA
    }
    
    # now mirror the updated order of counts in design
    if (!is.null(design))
      design <- design[match(colnames(counts), design[,libID_col]),]
    
    # normalize data by rows
    if (is.null(norm.method))
      norm.method <- "identity"
    counts <-
      t(apply(counts, MARGIN=1, FUN=match.fun(norm.method, descend=FALSE))) # normalize counts
    
    # make a data frame of colors for labeling columns
    if (!is.null(color_by_var)) {
      
      if (is.null(color_by_var_levels)) {
        color_by_var_levels <- replicate(n=length(color_by_var), NULL)
      } else if (!is.list(color_by_var_levels)) {
        color_by_var_levels <- list(color_by_var_levels)
      }
      
      if (is.null(my_var_colors)) {
        my_var_colors <- replicate(n=length(color_by_var), NULL)
      } else if (!is.list(my_var_colors)) {
        my_var_colors <- list(my_var_colors)
      }
      
      # create data frame, with columns named by variable name,
      # containing vector of colors to pass to heatmap3 as ColSideColors
      plot_colors <- list()
      
      for (i in 1:length(color_by_var)) {
        
        if (is.numeric(design[,color_by_var[i]])) {
          
          # I should update this to use my values2colors function
          
          if (is.null(color_by_var_levels[[i]])) color_by_var_levels[[i]] <- 10
          if (is.null(my_var_colors[[i]])) {
            my_var_colors[[i]] <-
              viridis::viridis_pal()(color_by_var_levels[[i]])
          }
          
          if (length(my_var_colors[[i]]) != color_by_var_levels[[i]])
            my_var_colors[[i]] <-
              colorRampPalette(my_var_colors[[i]])(color_by_var_levels[[i]])
          
          plot_colors[[color_by_var[i]]] <-
            my_var_colors[[i]][
              findInterval(
                design[,color_by_var[i]], 
                vec=seq(min(design[,color_by_var[i]]), max(design[,color_by_var[i]]),
                        length.out=color_by_var_levels[[i]]+1),
                rightmost.closed=TRUE)]
        } else {
          if (is.null(my_var_colors[[i]])) {
            my_var_colors[[i]] <-
              RColorBrewer::brewer.pal(
                length(unique(design[,color_by_var[i]])), "Set1")
          }
          if (is.null(color_by_var_levels[[i]])) {
            if (is.factor(design[,color_by_var[i]])) {
              color_by_var_levels[[i]] <- levels(design[,color_by_var[i]])
            } else {
              color_by_var_levels[[i]] <- as.character(unique(design[,color_by_var[i]]))
            }
          }
          my_var_colors[[i]] <- my_var_colors[[i]][1:length(color_by_var_levels[[i]])]
          names(my_var_colors[[i]]) <- color_by_var_levels[[i]]
          
          plot_colors[[color_by_var[i]]] <-
            my_var_colors[[i]][
              design[,color_by_var[i]]]
        }
      }
      plot_colors <- do.call(cbind, plot_colors)
    }
    
    # open plotting device
    if (!is.null(filename)) {
      pdf(file=filename, w=plotdims[1], h=plotdims[2])
      on.exit(dev.off()) # close plotting device on exit
    } else quartz(w=plotdims[1], h=plotdims[2])
    
    # generate heatmap
    if (exists("plot_colors")) {
      heatmap3::heatmap3(
        counts, scale=scale, col=my_heatmap_cols,
        Rowv=row_dendro, Colv=col_dendro,
        margins=c(8,12),
        ColSideColors=plot_colors,
        ...)
    } else {
      heatmap3::heatmap3(
        counts, scale=scale, col=my_heatmap_cols,
        Rowv=row_dendro, Colv=col_dendro,
        margins=c(8,12),
        ...)
    }
    
    # generate legends in separate plotting devices
    # (because it's a pain to try to allow multiple legends within a single plot)
    if (add_legend)
      if (length(color_by_var) > 1) {
        for (i in 1:length(color_by_var)) {
          # open plotting device
          if (!is.null(filename)) {
            pdf(
              file=stringr::str_replace(
                filename, "\\.pdf$", paste0(".legend.", color_by_var[i], ".pdf")),
              w=2, h=3)
            on.exit(dev.off(), add=TRUE) # close plotting device on exit
          } else quartz(w=2, h=3)
          
          plot.new()
          
          if (is.numeric(design[,color_by_var[i]])) {
            legend_color_labels <- function(x) {
              x <- as.character(x)
              y <- ""
              for (j in 2:(length(x)-1)) {
                y[(j-1)*2] <- x[j]
                y[(j-0.5)*2] <- ""
              }
              return(y)
            }
            
            plotrix::color.legend(
              legend=
                legend_color_labels(
                  round(seq(min(design[,color_by_var[i]]),
                            max(design[,color_by_var[i]]),
                            length.out=color_by_var_levels[[i]]+1),
                        round(-log10(diff(range(design[,color_by_var[i]])))+2))),
              rect.col=my_var_colors[[i]],
              gradient="y",
              xl=0, xr=0.9, yb=0, yt=1.5)
          } else {
            legend(
              x="topleft",
              legend=names(my_var_colors[[i]]), fill=my_var_colors[[i]],
              bty="n", xpd=TRUE)
          }
        }
      } else if (length(color_by_var)==1) {
        if (is.numeric(design[,color_by_var[1]])) {
          legend_color_labels <- function(x) {
            x <- as.character(x)
            y <- ""
            for (j in 2:(length(x)-1)) {
              y[(j-1)*2] <- x[j]
              y[(j-0.5)*2] <- ""
            }
            return(y)
          }
          
          plotrix::color.legend(
            legend=
              legend_color_labels(
                round(seq(min(design[,color_by_var[1]]),
                          max(design[,color_by_var[1]]),
                          length.out=color_by_var_levels[[1]]+1),
                      round(-log10(diff(range(design[,color_by_var[1]])))+2))),
            rect.col=my_var_colors[[1]],
            gradient="y",
            xl=0.6, xr=0.7, yb=0.9, yt=1.1,)
        } else {
          legend(legend=names(my_var_colors[[1]]), fill=my_var_colors[[1]],
                 bty="n", x=leg_x, y=leg_y, xpd=TRUE)
        }
      }
  }
