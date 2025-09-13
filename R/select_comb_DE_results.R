

#' Title
#'
#' @param DE_results_comb
#'
#' @returns
#' @export
#'
#' @examples
select_comb_DE_results <- function(DE_Results_comb, X=NULL, Y=NULL, type = "L2FC"){

  if(! type %in% c("L2FC", "PVAL", "QVAL")){
    stop(
      paste0(
        "\"",
        type,
        "\" is not a valid column type. \n,
          Please specify one of: \"L2FC\", \"PVAL\", \"QVAL\""
      )
    )
  }
  if(!is.null(X)){
    if(!is.character(X)){stop("X name not valid. Please specify X with a character string ")}
    if(paste0(type, "_", X) %in% colnames(DE_Results_comb)){
      X <- paste0(type, "_", X)
    } else{
      stop(
        paste0(
          X,
          " not found in data. Did you misspell it? ",
          "\nAvailable columns: ",
          apply(
            stringr::str_split(grep(paste0(type, "_"), colnames(DE_Results_comb), value=TRUE), "_", simplify=TRUE)[,-1],
            FUN = paste0,
            collapse="_",
            MARGIN = 1
          )
        )
      )
    }
  }

  if(!is.null(Y)){
    if(!is.character(Y)){stop("Y name not valid. Please specify Y with a character string ")}
    if(paste0(type, "_", Y) %in% colnames(DE_Results_comb)){
      Y <- paste0(type, "_", Y)
    } else{
      stop(
        paste0(
          Y,
          " not found in data. Did you misspell it? ",
          "\nAvailable columns: ",
          apply(
            stringr::str_split(grep(paste0(type,"_"), colnames(DE_Results_comb), value=TRUE), "_", simplify=TRUE)[,-1],
            FUN = paste0,
            collapse="_",
            MARGIN = 1
          )
        )
      )
    }
  }


  if(is.null(X)){
    if(length(grep(paste0(type, "_"), colnames(DE_Results_comb))) > 2){
      message(
        "More than three results available. Assigning first free to X"
      )
    }
    tryCatch({
      X = not_in(
        grep(paste0(type, "_"), colnames(DE_Results_comb), value=TRUE),
        Y
      )[1]
    },
    error = function(e){
      message("X: could not assign a column name ")
    })
  }

  if(is.null(Y)){
    if(length(grep(paste0(type, "_"), colnames(DE_Results_comb))) > 2){
      message(
        "More than three results available. Assigning next free to Y"
      )
    }
    tryCatch({
      Y = not_in(
        grep(paste0(type, "_"), colnames(DE_Results_comb), value=TRUE),
        X
      )[1]
    },
    error = function(e){
      message("Y: could not assign a column name ")
    })
  }

  return(c(X, Y))

}
