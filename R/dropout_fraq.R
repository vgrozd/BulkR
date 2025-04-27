#' Title
#'
#' @param counts Numeric vector, counts of a single cell/sample
#' @param expected_counts Numeric vector, expected counts per gene/feature in the same order as counts
#' @param res Numeric, resolution (number of bins for expected count bins to calculate dropout fraction by), b/w 3 and 1e6
#' @param na.rm Logical, Ignore NAs in counts?
#' @returns dropout fraction per expected_counts bin
#' @export
#'
#' @examples dropout_fraq(Seurat@assays$RNA@counts[1,], expected_counts(Seurat@assays$RNA@counts))
dropout_fraq <- function(counts, expected_counts, res=1e3, na.rm = TRUE){

  stopifnot(is.numeric(res), res>3, res<1e6)
  ExpBin <- cut(
    expected_counts,
    breaks = seq(
      from = min(expected_counts),
      to = max(expected_counts),
      length.out = res+1
    ),
    include.lowest = TRUE,
    right = TRUE
  )

  return(
    data.frame(
      ExpBin = ExpBin,
      Counts = counts
    ) %>%
      group_by(ExpBin) %>%
      mutate(DropoutFraq = 1-mean(Counts>0, na.rm = na.rm)) %>%
      pull(DropoutFraq)
  )


}
