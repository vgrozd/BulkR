

#' Title
#'
#' @param DE_results_comb
#'
#' @returns
#' @export
#'
#' @examples
plot_L2FC_scatter <- function(DE_Results_comb, x=NULL, y=NULL){

  axes = grep("L2FC_", colnames(DE_Results_comb), value=TRUE)[1:2]
  p1 = ggplot(DE)
}
