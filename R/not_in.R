#' Title
#'
#' @param a
#' @param b
#'
#' @returns
#' @export
#'
#' @examples
not_in <- function(a_not_in, b){
  return(a[which(! a_not_in, b)])
}
