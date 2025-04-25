#' Title
#'
#' @param counts
#'
#' @returns An integer of expected counts
#' @export
#'
#' @examples expected_counts(Seurat@assays$RNA@counts)
expected_counts <- function(counts){
  expected <- rowSums(Counts)
  expected <- expected_RPM*1e6/sum(expected_RPM)
}
