#' Title
#'
#' @param counts A count vector/matrix to plot
#' @param expected_counts Vector of expected counts
#' @param dropout_fraq Array/vector of dropout fractions
#' @param per_cell Plot each single cell or average?
#' @param logistic_fit Fit a logistic curve?
#' @param fit.res Number of bins for the logistic curve
#' @param color_by_libsize Logical, plot line colors by lib size?
#' @param total_counts A vector of library sizes for each cell
#' @param force_sample Logical/numeric: Sample few cells to plot?
#' @param seed Numeric: Random sampling seed
#' @param ... Graphic parameters for the plot lines
#'
#' @returns Plot
#' @export
#'
#' @examples plot_dropout_vs_magnitude(expected_counts(Seurat@assays$RNA@counts),
#' @examples  dropout_fraq(
#'              Seurat@assays$RNA@counts,
#' @examples    expected_counts(Seurat@assays$RNA@counts)),
#' @examples    per_cell = FALSE,
#' @examples    log = TRUE
#' @examples  )


plot_dropout_vs_magnitude <- function(
    counts = NULL,
    expected_counts=NULL,
    dropout_fraq=NULL,
    per_cell = FALSE,
    logistic_fit = TRUE,
    fit.res = 300,
    color_by_libsize = TRUE,
    total_counts = NULL,
    force_sample = TRUE,
    seed = 142,
    ...
  ){

  if(
    !is.null(counts) & any(!is.null(c(expected_counts, dropout_fraq))) |
      is.null(counts) & any(is.null(c(expected_counts, dropout_fraq)))
  ){
    stop("Please provide on of both: Count matrix, or expected counts + dropout fractions! ")
  }

  if(!is.null(counts)){
    expected_counts = expected_counts(counts)
  }

  if(!per_cell & !is.null(dim(dropout_fraq))){
    dropout_fraq = dropout_fraq(counts, expected_counts)
    if(dim(dropout_fraq)[2]>1){
      dropout_fraq <- base::rowMeans(dropout_fraq, na.rm = TRUE)
    }
  }

  if(is.logical(force_sample)){
    if(isTRUE(force_sample) &
       any(
         c(
           dim(counts)[2] > 20,
           dim(dropout_fraq)[2] > 20
         )
       )
    ){
      tryCatch({
        set.seed(seed)
        counts <- counts[,sample(1:ncol(counts), size = 20, replace = FALSE)]
        message(paste0("More than 20 cells detected, reducing to 20 random cells... "))
      }, error = function(e){})
      tryCatch({
        set.seed(seed)
        dropout_fraq <- dropout_fraq[,sample(1:ncol(dropout_fraq), size = 20, replace = FALSE)]
        message(paste0("More than 20 cells detected, reducing to 20 random cells... "))
      }, error = function(e){})
    } else{}
  } else{
    if(is.numeric(force_sample) & force_sample > 0){
      if(any(
        c(
          dim(counts)[2] > force_sample,
          dim(dropout_fraq)[2] > force_sample
        )
      )){
        tryCatch({
          set.seed(seed)
          counts <- counts[,sample(1:ncol(counts), size = force_sample, replace = FALSE)]
          message(paste0("More than ", force_sample, " cells detected, reducing to 20 random cells... "))
        }, error = function(e){})
        tryCatch({
          set.seed(seed)
          dropout_fraq <- dropout_fraq[,sample(1:ncol(dropout_fraq), size = force_sample, replace = FALSE)]
          message(paste0("More than ", force_sample, " cells detected, reducing to 20 random cells... "))
        }, error = function(e){})
      } else{
        stop("More samples to draw specified than number of samples available")
      }
    } else{
      stop("'force_sample': Please provide either logical or a valid numerical value")
    }
  }

  if(!is.null(counts)){
    dropout_fraq = dropout_fraq(counts, expected_counts)
  }

  if(!per_cell){
    if(dim(dropout_fraq)[2]>1){
      dropout_fraq <- base::rowMeans(dropout_fraq, na.rm = TRUE)
    }
  }

  if(is.null(total_counts)){
    total_counts = rep(1e3, dim(as.matrix(dropout_fraq))[2])
  }

  if(!logistic_fit){

    return(


      ggplot2::ggplot(
        data.frame(
          expected_counts = rep(
            expected_counts,
            dim(as.matrix(dropout_fraq))[2]
          ),
          dropout_fraq = as.vector(dropout_fraq),
          Cell = rep(
            paste0("Cell_", 1:dim(as.matrix(dropout_fraq))[2]),
            each = length(expected_counts)
          )
        )
      ) +
        ggplot2::aes(expected_counts, dropout_fraq, fill = Cell) +
        (if(color_by_libsize){
          ggplot2::geom_smooth(aes(col=total_counts), se = FALSE, ...)
        } else{
          ggplot2::geom_smooth(col="blue", se = FALSE, ...)
        }) +
        (if(color_by_libsize){
          ggplot2::scale_color_viridis_b()
          } else {NULL}) +
        ggplot2::theme_classic()
    )

  } else {

    newdata = data.frame(
      expected_counts = seq(
        from = min(expected_counts, na.rm = TRUE),
        to = max(expected_counts, na.rm = TRUE),
        length.out = fit.res
      )
    )

    dat <- apply(
      as.matrix(dropout_fraq),
      FUN = function(x, expcs = expected_counts, nwdt = newdata){
        fit <- stats::glm(
          x ~ expected_counts,
          data = data.frame(
            x = x,
            expected_counts = expcs
          ),
          family = quasibinomial
        )

        return(
          predict(
            fit,
            nwdt,
            type = "response"
          )
        )

      },
      MARGIN = 2
    )


    return(

      ggplot2::ggplot(
        data.frame(
          expected_counts = rep(
            newdata$expected_counts,
            dim(as.matrix(dat))[2]
          ),
          dropout_fraq = as.vector(dat),
          Cell = rep(
            paste0("Cell_", 1:dim(as.matrix(dat))[2]),
            each = length(newdata$expected_counts)
          ),
          TotalCounts = rep(
            total_counts,
            each = length(newdata$expected_counts)
          )
        )
      ) +
        ggplot2::aes(expected_counts, dropout_fraq, fill=Cell) +
        (if(color_by_libsize){
          ggplot2::geom_line(ggplot2::aes(col=log10(TotalCounts)), size=0.5, alpha = ...)
        } else {
          ggplot2::geom_line(size=0.5, col = "blue", alpha = ...)
        }) +
        (if(color_by_libsize){
          ggplot2::scale_color_viridis_b()
        } else {
          NULL
        }) +
        ggplot2::theme_classic()
    )
  }
}



