#' Title
#'
#' @param counts
#'
#' @returns An integer of expected counts
#' @export
#'
#' @examples expected_counts(Seurat@assays$RNA@counts)
expected_counts <- function(counts, factor=1e6){
  return(
    rowSums(counts) * factor / sum(rowSums(counts))
  )
}
