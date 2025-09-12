#' Title
#'
#' @param DE_results
#' @param format
#'
#' @returns
#' @export
#'
#' @examples
unify_DE_format <- function(DE_results, format=NULL){

  # Define format

  if(is.null(format)){

    # TODO: Auto-detect
    stop("Auto-detect not yet available, please specify 'format'")
  } else{

    if( ! format %in% c("DESeq2", "MAST")){
      stop("Format not familiar!")
    }
    if(format=="DESeq2"){
      # Check if results or dataframe
      # Convert
    }

    if(format=="MAST"){
      # Convert
    }

  }
}
