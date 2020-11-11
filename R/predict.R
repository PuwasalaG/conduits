#' Predict method for conditional cross-correlation fits
#'
#' Compute predictions for conditional cross-correlations at lags $k=1,2,...$ from the fitted models
#' in \code{conditional_ccf} function.
#'
#' @param object an object returned from \code{conditional_ccf} function
#' @param new_data  a tibble for which the conditional cross-correlations need to be predicted.
#' This tibble should contain all the predictors passed to \code{conditional_ccf} along with the Timestamp.
#' @param ... further arguments passed to or from other methods
#'
#' @return returns a tibble with the predicted conditional cross-correlations
#' at lags $k=1,2,...$
#'
#' @author Puwasala Gamakumara
#'
#' @export
predict.conditional_ccf <- function(object, new_data, ...){

  k <- object$other$k

  ccf_gam_fit <- object$data_visualise$conditional_ccf$ccf_gam_fit

  data_predict_ccf <- matrix(0, ncol = length(k), nrow = nrow(new_data))

  for (i in k) {

    data_predict_ccf[,i] <- stats::predict.glm(ccf_gam_fit[[i]], newdata = new_data,
                                               type = "response")

  }

  colnames(data_predict_ccf) <- paste("k = ", k, sep = "")
  data_predict_ccf <- data_predict_ccf %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Timestamp = new_data$Timestamp) %>%
    dplyr::select(.data$Timestamp, tidyselect::everything())

  data_predict_ccf <- new_data %>%
    dplyr::left_join(data_predict_ccf, by = "Timestamp")

  return(data_predict_ccf)

}
