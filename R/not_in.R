#' Title
#'
#' @param a
#' @param b
#'
#' @returns
#' @export
#'
#' @examples
not_in <- function(a,not_in_b){
  return(a[which(! a %in% not_in_b)])
}
