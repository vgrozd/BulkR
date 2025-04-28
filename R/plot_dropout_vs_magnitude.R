#' Title
#'
#' @param expected_counts Vector of expected counts
#' @param dropout_fraq Array/vector of dropout fractions
#' @param per_cell Plot each single cell or average?
#' @param logistic_fit Fit a logistic curve?
#' @param fit.res Number of bins for the logistic curve
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
    expected_counts,
    dropout_fraq,
    per_cell = TRUE,
    logistic_fit = TRUE,
    fit.res = 300,
    total_counts = NULL,
    ...
  ){

  if(!per_cell & !is.null(dim(dropout_fraq))){
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
        ggplot2::aes(expected_counts, dropout_fraq, fill=Cell) +
        ggplot2::geom_smooth(..., se = TRUE) +
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
        ggplot2::geom_line(ggplot2::aes(col=log10(TotalCounts)), size=0.5, alpha = ...) +
        ggplot2::theme_classic()
    )

  }
}



