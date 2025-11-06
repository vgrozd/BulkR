#' Serial SVA DESeq2 analysis
#'
#' @param DDS
#' @param n.sv
#' @param parallel
#'
#' @returns
#' @export
#'
#' @examples
serial_sva_DESeq2 <- function(DDS, n.sv = NULL, parallel = TRUE, qval = 0.05, l2fc = 0, nthr = NULL){

  if(!is.null(nthr)){
    BiocParallel::register(BiocParallel::MulticoreParam(nthr))
  }

  if(is.null(DESeq2::sizeFactors(DDS))){
    message(
      paste0(
        "Provided DDS object was not DESeq`ed. Running DESeq2 now..."
      )
    ); cat("\n")
    DDS <- DESeq2::DESeq(DDS, parallel = parallel, quiet = TRUE)
  } else {
    message(
      paste0(
        "DESeq object has been DESeq`ed already. Using supplied size facors... "
      )
    )
    message(
      paste0(
        "DESeq2 design: ",
        paste0(as.character(DESeq2::design(DDS)), collapse = " "),
        "... "
      )
    ); cat("\n")
  }

  dat  <- DESeq2::counts(DDS, normalized = TRUE)
  message(
    paste0(
      "Count matrix with ",
      dim(dat)[2],
      " samples and ",
      dim(dat)[1],
      " features. "
    )
  )
  idx  <- rowMeans(dat) > 1
  message(
    paste0(
      "Removing ",
      sum(rowMeans(dat) <= 1),
      " features due to low counts..."
    )
  ); cat("\n")
  dat  <- dat[idx, ]
  mod  <- model.matrix(DESeq2::design(DDS), SummarizedExperiment::colData(DDS))
  mod0 <- model.matrix(~ 1, SummarizedExperiment::colData(DDS))

  if(is.null(n.sv)){
    message(
      "Running SVA with automatic assessment of SVs..."
    ); cat("\n")
  } else{
    message(
      if(n.sv == 1) {"Running SVA with 1 SV..."} else{
        paste0(
          "Running SVA with 1-",
          n.sv,
          " SVs..."
        )
      }
    ); cat("\n")
  }

  suppressMessages(nSVs_Auto <- sva::svaseq(dat, mod, mod0, n.sv = NULL)$n.sv)
  svseq <- sva::svaseq(dat, mod, mod0, n.sv = n.sv)

  for(i in 1:ncol(svseq$sv)){
    DDS@colData[[
      paste0(
        "SV",
        i
      )
    ]] <- svseq$sv[,i]
  }

  rm(mod, mod0, dat, idx, i)


  DDS_with_serial_SVA <- list()

  Benchmark <- data.frame(
    Name = "No_SVA",
    nSVs = 0,
    Elapsed=system.time(
      suppressWarnings(
        DESeq2::DESeq(DDS)
      )
    )["elapsed"]
  )


  message(
    "Running DESeq with SVA with increasing numbers of SVs now... "
  ); cat("\n")

  for(i in 1:ncol(svseq$sv)){

    message(
      paste0(
        "DESeq2 with ", i, if(i==1){" SV... "} else{" SVs... "}
      )
    )

    DDS.tmp <- DDS
    design.tmp <- as.character(
      DESeq2::design(DDS.tmp)
    )
    design.tmp <- design.tmp[
      which(design.tmp != "~")
    ]
    DESeq2::design(DDS.tmp) <- as.formula(
      paste0(
        "~ ",
        paste0(paste0("SV", 1:i), collapse = " + "),
        " + ",
        design.tmp
      )
    )
    message(
      paste0(
        "DESeq2 Design: ",
        paste0(as.character(DESeq2::design(DDS.tmp)), collapse = " "),
        " ... "
      )
    ); cat("\n")

    Benchmark <- rbind(
      Benchmark,
      data.frame(
        Name = paste0("SVA_", i, "_SVs"),
        nSVs = i,
        Elapsed=system.time(
          DDS_with_serial_SVA[[i]] <- DESeq2::DESeq(DDS.tmp, parallel=parallel, quiet = TRUE)
        )["elapsed"]
      )
    )

    rm(DDS.tmp, design.tmp)
  }


  DDS_with_serial_SVA <- c(
    list(DDS),
    DDS_with_serial_SVA
  )

  names(DDS_with_serial_SVA) <- c(
    "No_SVA",
    paste0(
      "SVA_",
      1:n.sv,
      "_SVs"
    )
  )

  DESeq2_Results <- lapply(
    DDS_with_serial_SVA,
    function(x){
      DESeq2::results(x, alpha = qval, lfcThreshold = l2fc)
    }
  )

  Results <- data.frame(
    nSVs = 0:(length(DDS_with_serial_SVA)-1),
    nDEGs = sapply(
      DESeq2_Results,
      FUN = function(x){
        x |>
          BulkR::significant_genes() |>
          length()
      }
    ),
    DEGs = sapply(
      DESeq2_Results,
      FUN = function(x){
        x |>
          BulkR::significant_genes() |>
          paste0(collapse=",")
      }
    )
  )

  Serial_SVA <- list(
    DDS = DDS_with_serial_SVA,
    DESeq2_Results = DESeq2_Results,
    Results = Results,
    sva = svseq,
    SVs = svseq$sv,
    nSVs_Auto = nSVs_Auto,
    Benchmark = Benchmark
  )

  class(Serial_SVA) <- "Serial_SVA"

  return(
    Serial_SVA
  )

}
