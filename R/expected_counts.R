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

expected_counts <- function(counts, factor=1e6, log=FALSE, pseudocount=1){
  return(
    if(!log){
      Matrix::rowSums(counts, na.rm = TRUE) * factor / sum(Matrix::rowSums(counts, na.rm = TRUE), na.rm = TRUE)
    } else{
      (Matrix::rowSums(counts, na.rm = TRUE)+ pseudocount) * factor / sum(Matrix::rowSums(counts, na.rm = TRUE) + pseudocount, na.rm = TRUE)
    }

  )
}
