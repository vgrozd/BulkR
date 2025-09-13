

combine_results <- function(..., filter=TRUE, fdr=NULL, l2fc=NULL, type="list"){

  if(! type %in% c("list", "df", "tbl")){
    stop(
      "Unvalid return type. Please specify one of: \"list\", \"df\", \"tbl\""
    )
  }

  if(is.null(names(list(...)))){
    results_list = setNames(
      list(...),
      nm = sapply(substitute(list(...)), deparse)[-1]
    )
  }else{
    results_list = setNames(
      list(...),
      nm = names(list(...))
    )
  }

  if(filter){

    if(is.null(fdr)){
      message("Using default FDR value qval < 0.05")
      fdr = 0.05
    }
    if(length(fdr)==1){
      if(length(results_list) > 1){
        message(paste0("Using FDR qval<", fdr, " for all DE results"))
        fdr = rep(fdr, length(results_list))
      }
    } else{
      fdr = fdr[which(!is.na(fdr) & !is.null(fdr) & is.numeric(fdr) & (fdr>0) & (fdr<=1))]
      if(length(fdr)==length(results_list)){
        message("Using specified FDR values... ")
      } else {
        message(
          "Specified FDRs not same length as results. Reverting to default FDR value qval < 0.05... "
        )
        fdr = rep(0.05, length(results_list))
      }
    }

    if(is.null(l2fc)){
      message("Using default L2FC value <0, >0")
      l2fc = 0
    }
    if(length(l2fc)==1){
      if(length(results_list) > 1){
        message(paste0("Using L2FC value ", l2fc, " for all DE results"))
        l2fc = rep(l2fc, length(results_list))
      }
    } else{
      l2fc = l2fc[which(!is.na(l2fc) & !is.null(l2fc) & is.numeric(l2fc))]
      if(length(l2fc)==length(results_list)){
        message("Using specified L2FC values... ")
      } else {
        message(
          "Specified L2FCs not same length as results. Reverting to default L2FC value 0 ... "
        )
        l2fc = rep(0.05, length(results_list))
      }
    }

    fdr = fdr[!is.null(results_list)]
    l2fc = l2fc[!is.null(results_list)]

    results_list = results_list[!is.null(results_list)]
    message("---")
    tryCatch({
      nm = "?"
      res_l <- lapply(
        seq_along(results_list),
        function(i){
          nm <<- names(results_list)[i]
          significant_results(results_list[[i]], qval = fdr[i], l2fc = l2fc[i])
        }
      )
      res_l <- setNames(res_l, nm = names(results_list))
    }, error=function(e){
      stop(paste0(
        "Could not extract significant results from DE results \"", nm, "\"",
        "\nDid you mistype an argument '", nm, "'? "
      ))
      # TODO check here what happens if a result is null or doesn't have sign results
    }
    )

  } else{

    results_list = results_list[!is.null(results_list)]
    message("---")
    tryCatch({
      nm = "?"
      res_l <- lapply(
        seq_along(results_list),
        function(i){
          nm <<- names(results_list)[i]
          unify_DE_format(results_list[[i]])
        }
      )
      res_l <- setNames(res_l, nm = names(results_list))
    }, error=function(e){
      stop(paste0(
        "Could not extract unify format of DE results: \"", nm, "\""),
        ,
        "\nDid you mistype an argument '", nm, "'? "


      )
      # TODO check here what happens if a result is null or doesn't have sign results
    }
    )

  }


  if(type=="list") {
    return(res_l)
  }
  if(type %in% c("df", "tbl")){
    message(
      "Data frame and tibble return types not yet implemented"
    )
    return(invisible())
  }




}




