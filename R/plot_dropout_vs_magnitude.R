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
            log10(newdata$expected_counts+0.5),
            dim(as.matrix(dat))[2]
          ),
          dropout_fraq = as.vector(dat),
          Cell = rep(
            paste0("Cell_", 1:dim(as.matrix(dat))[2]),
            each = length(newdata$expected_counts)
          )
        )
      ) +
        ggplot2::aes(expected_counts, dropout_fraq, fill=Cell) +
        ggplot2::geom_line(..., se = FALSE) +
        ggplot2::theme_classic()
    )

  }
}

summary(colSums(Counts))
Counts3 <- Counts[,which(colSums(Counts) > 1000 & colSums(Counts) < 1100)]
plot_dropout_vs_magnitude(
  expected_counts(Counts3),
  dropout_fraq(Counts3, expected_counts(Counts3)), per_cell = FALSE, log = TRUE
)

