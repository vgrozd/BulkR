plot_dropout_vs_magnitude <- function(
    expected_counts,
    dropout_fraq,
    per_cell = TRUE,
    logistic_fit = TRUE,
    ...
  ){

  if(!logistic_fit){

    return(

      ggplot(
        data.frame(
          expected_counts = expected_counts,
          dropout_fraq = dropout_fraq
        )
      ) +
        aes(expected_counts, dropout_fraq) +
        geom_smooth(...)
    )

  }



  # Fit logistic model
  fit <- glm(DropoutFraq ~ ExpRPM_Log10, data = test, family = quasibinomial)
  newdata <- data.frame(ExpRPM_Log10 = seq(min(test$ExpRPM_Log10), max(test$ExpRPM_Log10), length.out = 300))
  newdata$DropoutPred <- predict(fit, newdata = newdata, type = "response")

  # Plot
  ggplot(test, aes(x = ExpRPM_Log10, y = DropoutFraq)) +
    geom_line(data = newdata, aes(x = ExpRPM_Log10, y = DropoutPred), color = "blue", size = 0.5) +
    theme_classic()
}
