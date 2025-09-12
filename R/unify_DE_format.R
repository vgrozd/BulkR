#' Title
#'
#' @param DE_results
#' @param format
#'
#' @returns
#' @export
#'
#' @examples

unify_DE_format <- function(DE_results, format=NULL, inplace = FALSE){

  # Define format

  if(is.null(format)){

    # Auto-detect

    if(any(class(DE_results)=="DESeqResults")){
      format <- "DESeq2"
    }

    if(any(class(DE_results) %in% c("tbl_df", "tbl", "data.frame"))){
      format <- "DF"
    }

    if(format=="DF"){

      if(all(
        c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj") %in%
          colnames(DE_results)
      )){format <- "DESeq2"}

      if(all(
        c("gene", "p_value", "model_log2FC", "ci.hi", "ci.lo", "fdr") %in%
        colnames(DE_results)
      )){format <- "MAST"}
    }

    if(is.null(format)){stop("Format could not be detected, please specify 'format'")} else{message(paste0("Auto-detected DE results format: ", format))}
  } else{

    if( ! format %in% c("DESeq2", "MAST", "DF")){
      stop("Unknown format! ")
    }

  }

  # Check format

  if(
    format == "DESeq2" &
    any(! c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj") %in%
        colnames(DE_results))
  ){
    stop("Format set as DESeq2 but data is not compatible with this format")
  }

  if(
    format == "MAST" &
    any(! c("gene", "p_value", "model_log2FC", "ci.hi", "ci.lo", "fdr") %in%
        colnames(DE_results))
  ){
    stop("Format set as MAST but data is not compatible with this format")
  }

  # Unify format

  if(format=="DESeq2"){
    DE_results <- as.data.frame(DE_results)
    .df = data.frame(
        Gene = rownames(DE_results),
        L2FC = DE_results$log2FoldChange,
        PVAL = DE_results$pvalue,
        QVAL = DE_results$padj
      )
    if (inplace) {
      message(paste0("Modified ", deparse(substitute(DE_results)), " in place! "))
      assign(deparse(substitute(DE_results)), .df, envir = parent.frame())
      return(invisible())
    } else {
      return(.df)
    }
  }

  if(format=="MAST"){
    return(
      data.frame(
        Gene = DE_results$gene,
        L2FC = DE_results$model_log2FC,
        PVAL = DE_results$p_value,
        QVAL = DE_results$fdr
      )
    )
  }

  if(format=="DF"){
    message("Results in an uknown data.frame format. Trying to identify columns...")
    if(FALSE){return()}else{}
  }

  stop("DE results format could not be unified! ")

}

