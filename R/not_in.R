#' Title
#'
#' @param a
#' @param b
#'
#' @returns
#' @export
#'
#' @examples
not_in <- function(a,b){
  return(a[which(! a %in% b)])
}
