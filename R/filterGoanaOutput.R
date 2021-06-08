#' Filter annotation enrichment output from goana
#'
#' This function filters the output of \code{goana} to include only terms that are either more specific or more strongly enriched. It uses the ancestor-dependent relationships of the Gene Ontology database to identify terms that are less specific, and culls them unless they are more strongly enriched in the analysis than the descendant term.
#' @param goanaResult a data frame containing the results of call to \code{goana}
#' @param idCol character or numeric, the column of \code{goanaResult} containing the ontology identifiers. Defaults to "GO_ID"
#' @param orderBy character or numeric, the column of \code{goanaResult} that contains information on the enrichment of the term. This is passed to \code{dplyr::arrange}. Defaults to "P.DE"
#'
#' @export
#' @return A data frame that is a subset of the input \code{goanaResult}, retaining only rows for terms that are more strongly enriched than any descendant terms
#' @usage \code{filterGoanaOutput(goanaResult, idCol = "GO_ID", orderBy = "P.DE")}
filterGoanaOutput <-
  function(goanaResult, idCol = "GO_ID", orderBy = "P.DE") {

    # load GO term ancestor-descendant relationships
    if (requireNamespace("GO.db", quietly = TRUE)) {
      goAncestorsAll <-
        c(as.list(GOBPANCESTOR), as.list(GOCCANCESTOR), as.list(GOMFANCESTOR))
    } else stop("This function requires the Bioconductor package 'GO.db'. Please install it to use this function.")

    # sort the input data frame by column of choice (default: significance)
    goanaResultFiltered <-
      dplyr::arrange(goanaResult, !!rlang::sym(orderBy))

    rowsToDropList <- list()
    # iterate over each term
    for (rowTmp in 1:nrow(goanaResultFiltered)) {
      # fetch ancestor terms for focal term
      goAncestorsTmp <- goAncestorsAll[[goanaResultFiltered[[idCol]][rowTmp]]]

      # identify less enriched ancestor terms in sorted data frame
      rowsToDropTmp <-
        which(goanaResultFiltered[[idCol]] %in% goAncestorsTmp)

      # avoid dropping more strongly enriched rows, using the sort order
      rowsToDropTmp <- setdiff(rowsToDropTmp, 1:rowTmp)

      # if any less enriched ancestor terms found, add them to the list to be dropped at end
      if (length(rowsToDropTmp) > 0)
        rowsToDropList[[goanaResultFiltered[[idCol]][rowTmp]]] <-
        rowsToDropTmp
    }

    # determine unique terms to drop, and remove them
    rowsToDrop <- unique(unlist(rowsToDropList))
    if (length(rowsToDrop) > 0)
      goanaResultFiltered <- goanaResultFiltered[-rowsToDrop,]

    # put filtered data frame back in original order by matching idCol to the input data frame
    goanaResultFiltered <-
      goanaResultFiltered[
        order(
          match(goanaResultFiltered[[idCol]], goanaResult[[idCol]])),]

    goanaResultFiltered
  }
