plot_dropout_vs_magnitude <- function(
    expected_counts,
    dropout_fraq,
    per_cell = TRUE,
    logistic_fit = TRUE,
    fit.res = 300,
    ...
  ){

  if(!per_cell & !is.null(dim(dropout_fraq))){
    if(dim(dropout_fraq)[2]>1){
      dropout_fraq <- base::rowMeans(dropout_fraq, na.rm = TRUE)
    }
  }

  if(!logistic_fit){

    return(


      ggplot2::ggplot(
        data.frame(
          expected_counts = rep(
            log10(expected_counts+0.5),
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

    apply(
      as.matrix(dropout_fraq),
      FUN = function(x){
        fit <- stats::glm(
          x ~ .expected_counts,
          data = data.frame(
            expected_counts = rep(
              .expected_counts,
              dim(as.matrix(x))[2]
            ),
            dropout_fraq = as.vector(x)
          ),
          family = quasibinomial
        )



        return(
          predict(fit, newdata = data.frame(expected_counts=seq(min(expected_counts), max(expected_counts), length.out = fit.res)), type = "response")
        )

      },
      .expected_counts=expected_counts, .fit.res = fit.res,
      MARGIN = 2
    )

    fit <- stats::glm(
      dropout_fraq ~ expected_counts,
      data = data.frame(
        expected_counts = rep(
          expected_counts,
          dim(as.matrix(dropout_fraq))[2]
        ),
        dropout_fraq = as.vector(dropout_fraq)
      ),
      family = quasibinomial
    )

   expected_counts = seq(min(expected_counts), max(expected_counts), length.out = fit.res)

   dropout_fraq <- do.call(
     cbind,
     list(apply(
       as.matrix(expected_counts),
       FUN = function(x){
         predict(fit, newdata = data.frame(x), type = "response")
       },
       MARGIN = 2
     )
   ))

    return(

      ggplot2::ggplot(
        data.frame(
          expected_counts = rep(
            newdata$expected_counts,
            dim(as.matrix(newdata$dropout_fraq))[2]
          ),
          dropout_fraq = as.vector(newdata$dropout_fraq),
          Cell = rep(
            paste0("Cell_", 1:dim(as.matrix(newdata$dropout_fraq))[2]),
            each = dim(as.matrix(newdata$dropout_fraq))[2]
          )
        )
      ) +
        ggplot2::aes(expected_counts, dropout_fraq, fill=Cell) +
        ggplot2::geom_line() +
        ggplot2::theme_classic()
    )

  }
}

