#' Title
#'
#' @param counts Numeric vector, counts of a single cell/sample
#' @param expected_counts Numeric vector, expected counts per gene/feature in the same order as counts
#'
#' @returns dropout fraction per expected_counts bin
#' @export
#'
#' @examples dropout_fraq(Seurat@assays$RNA@counts[1,], expected_counts(Seurat@assays$RNA@counts))
dropout_fraq <- function(counts, expected_counts){

}
