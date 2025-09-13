

combine_results <- function(..., filter=TRUE, fdr=NULL, l2fc=NULL, type="list", significant_genes="all"){

  if(! type %in% c("list", "df", "tbl")){
    stop(
      "Unvalid return type. Please specify one of: \"list\", \"df\", \"tbl\""
    )
  }

  if(! significant_genes %in% c("all", "common", "unique")){
    stop(
      "Unvalid significant genes type. Please specify 'significant_genes=' as one of: \"all\", \"common\", \"unique\""
    )
  }

  # if the first argument is a list and has length > 1, treat it as the datasets
  if (length(list(...)) == 1 && is.list(list(...)[[1]]) && !is.data.frame(list(...)[[1]])) {
    results_list <- list(...)[[1]]


    if(is.null(names(results_list))){
      results_list = setNames(
        results_list,
        nm = paste0(sapply(substitute(list(...)), deparse)[-1][1], "_", seq_along(results_list))
      )
    }


  } else {

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

    if(filter){

      message("---")
      tryCatch({
        nm = "?"
        full_res_l <- lapply(
          seq_along(results_list),
          function(i){
            nm <<- names(results_list)[i]
            unify_DE_format(results_list[[i]])
          }
        )
        full_res_l <- setNames(full_res_l, nm = names(results_list))
      }, error=function(e){
        stop(paste0(
          "Could not extract unify format of DE results: \"", nm, "\""),
          ,
          "\nDid you mistype an argument '", nm, "'? "


        )
        # TODO check here what happens if a result is null or doesn't have sign results
      }
      )

      if(significant_genes=="all"){
        Genes = unique(
          unlist(
            sapply(
              res_l,
              function(x){
                return(x$Gene)
              }
            )
          )
        )
        df = data.frame(
          Gene = Genes
        )
        sapply(
          names(res_l),
          function(r){
            df[[paste0("L2FC_", r)]] <<- full_res_l[[r]]$L2FC[match(df$Gene, full_res_l[[r]]$Gene)]
            df[[paste0("PVAL_", r)]] <<- full_res_l[[r]]$PVAL[match(df$Gene, full_res_l[[r]]$Gene)]
            df[[paste0("QVAL_", r)]] <<- full_res_l[[r]]$QVAL[match(df$Gene, full_res_l[[r]]$Gene)]
          }
        )
        message("---")
        message("Returning data from all datasets of all genes significant in at least one of the datasets... ")
        return(df)
      }

      if(significant_genes=="common"){
        Genes = Reduce(intersect, lapply(res_l, `[[`, "Gene"))
        df = data.frame(
          Gene = Genes
        )
        sapply(
          names(res_l),
          function(r){
            df[[paste0("L2FC_", r)]] <<- full_res_l[[r]]$L2FC[match(df$Gene, full_res_l[[r]]$Gene)]
            df[[paste0("PVAL_", r)]] <<- full_res_l[[r]]$PVAL[match(df$Gene, full_res_l[[r]]$Gene)]
            df[[paste0("QVAL_", r)]] <<- full_res_l[[r]]$QVAL[match(df$Gene, full_res_l[[r]]$Gene)]
          }
        )
        message("---")
        message("Returning data for all datasets of all genes significant in at least one of the datasets... ")
        return(df)
      }

      if(significant_genes=="unique"){
        Genes_oth = unique(
          unlist(
            sapply(
              res_l[-1],
              function(x){
                return(x$Gene)
              }
            )
          )
        )
        Genes_uni = setdiff(
          res_l[[1]]$Gene,
          Genes_oth
        )



        df = data.frame(
          Gene = Genes_uni
        )
        sapply(
          names(res_l),
          function(r){
            df[[paste0("L2FC_", r)]] <<- full_res_l[[r]]$L2FC[match(df$Gene, full_res_l[[r]]$Gene)]
            df[[paste0("PVAL_", r)]] <<- full_res_l[[r]]$PVAL[match(df$Gene, full_res_l[[r]]$Gene)]
            df[[paste0("QVAL_", r)]] <<- full_res_l[[r]]$QVAL[match(df$Gene, full_res_l[[r]]$Gene)]
          }
        )
        message("---")
        message("Returning data for all datasets of all genes significant in at least one of the datasets... ")
        return(df)
      }


    }else{
      Genes = Reduce(intersect, lapply(res_l, `[[`, "Gene"))
      df = data.frame(
        Gene = Genes
      )
      sapply(
        names(res_l),
        function(r){
          df[[paste0("L2FC_", r)]] <<- res_l[[r]]$L2FC[match(df$Gene, res_l[[r]]$Gene)]
          df[[paste0("PVAL_", r)]] <<- res_l[[r]]$PVAL[match(df$Gene, res_l[[r]]$Gene)]
          df[[paste0("QVAL_", r)]] <<- res_l[[r]]$QVAL[match(df$Gene, res_l[[r]]$Gene)]
        }
      )
      return(df)
    }


  }




}




