#' Drop samples from EList object
#'
#' This function drops samples from an EList object, iterating over all components of the object. It uses the dimensions of the components to determine whether to drop rows, columns, etc., and thus will have difficulty if the number of samples is small.
#' @param EList an \code{EList} object containing expression values
#' @param keep_samples a numeric vector of the indices of samples to keep; alternately, a logical vector with length equal to number of libraries in EList, indicating which libraries to keep.
#' @param drop_samples a numeric vector of the indices of samples to drop; alternately, a logical vector with length equal to number of libraries in EList, indicating which libraries to drop. Ignored if \code{keep_samples} is provided.
#' @export
#' @details Either \code{keep_samples} or \code{drop_samples} must be specified.
#' @return an \code{EList} object with \code{drop_samples} removed from all components, or only \code{keep_samples} retained in all components.
#' @usage \code{trim_EList(EList, keep_samples, drop_samples=NULL)}
trim_EList <- function(EList, keep_samples, drop_samples=NULL) {
  if (!base::exists("keep_samples", where=environment(), inherits=FALSE) &
      is.null(drop_samples)) return(EList)
  
  nsamples <- ncol(EList[["E"]])
  
  if (!is.null(drop_samples)) {
    if (is.logical(drop_samples)) {
      if (length(drop_samples) != nsamples) stop("Vector of samples to drop does not match number of samples in EList object.")
      drop_samples <- which(drop_samples)
    }
    if (any(!(drop_samples %in% seq.int(nsamples))))
      stop("Vector of samples to drop does not match number of samples in EList object.")
    keep_samples <- (1:nsamples)[-drop_samples]
  } else {
    if (is.logical(keep_samples)) {
      if (length(keep_samples) != nsamples) stop("Vector of samples to retain does not match number of samples in EList object.")
      keep_samples <- which(keep_samples)
    }
    if (max(keep_samples) > nsamples)
      stop("Vector of samples to retain does not match number of samples in EList object.")
  }
  
  for (i in names(EList)) {
    if (i == "targets") {
      EList[[i]] <- EList[[i]][keep_samples,]
    } else if (i %in% c("E","weights")) {
      EList[[i]] <- EList[[i]][,keep_samples]
    } else if (i == "sample.weights") {
      EList[[i]] <- EList[[i]][keep_samples]
    } else if (i == "design") {
      EList[[i]] <- ordinal::drop.coef(EList[[i]][keep_samples,])
    } else {
      if (is.vector(EList[[i]])) {
        EList[[i]] <- EList[[i]][keep_samples]
      } else if ((dim(EList[[i]])[1] == nsamples) & (dim(EList[[i]])[2] != nsamples)) {
        EList[[i]] <- EList[[i]][keep_samples,]
      } else if ((dim(EList[[i]])[1] != nsamples) & (dim(EList[[i]])[2] == nsamples)) {
        EList[[i]] <- EList[[i]][,keep_samples]
      } else stop(paste("Cannot resolve dimensions of EList element", i))
    }
  }
  return(EList)
}
