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
#' @examples expected_counts(Seurat@assays$RNA@counts, factor=1e3, log=TRUE, pseudocount=0.5)

expected_counts <- function(counts, factor=1e6, log=FALSE, pseudocount=NULL){
  return(
    if(log){
      (Matrix::rowSums(counts)+if(is.null(pseudocount))1 else pseudocount) * factor / sum(Matrix::rowSums(counts)+if(is.null(pseudocount)) 1 else pseudocount)
    } else{
      Matrix::rowSums(counts) * factor / sum(Matrix::rowSums(counts))
    }

  )
}
