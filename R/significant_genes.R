#' Title
#'
#' @param .results
#' @param qval
#' @param l2fc
#'
#' @returns
#' @export
#'
#' @examples
significant_genes <- function(.results, qval=0.05, l2fc=0){
  .results <- unify_DE_format(.results)
  return(
    .results$Gene[
      which(
        .results$QVAL<qval,
        abs(.results$L2FC)>l2fc
      )
    ]
  )
}
