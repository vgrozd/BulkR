#' Title
#'
#' @param a
#' @param b
#'
#' @returns
#' @export
#'
#' @examples
not_in <- function(in_a,b){
  return(in_a[which(! in_a %in% b)])
}
