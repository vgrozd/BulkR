#' Title
#'
#' @param counts A counts matrix
#' @param factor Integer, a scaling factor
#' @param log Logical, log counts?
#'
#' @returns An integer of expected counts
#' @export
#'
#' @examples expected_counts(Seurat@assays$RNA@counts)
expected_counts <- function(counts, factor=1e6, log=FALSE, pseudocount=NULL){
  return(
    if(log){
      (rowSums(counts)+if(is.null(pseudocount))1 else pseudocount) * factor / sum(rowSums(counts)+if(is.null(pseudocount)) 1 else pseudocount)
    } else{
      rowSums(counts) * factor / sum(rowSums(counts))
    }

  )
}
