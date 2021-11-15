#' Estimating time delay between two sensors in a river system
#'
#' This function estimates the time that takes water to flow from an upstream location to a downstream
#' location conditional on the observed water-quality variables from the upstream sensor. That time lag is
#' defined as the lag that gives maximum cross-correlation conditional on upstream water-quality variables.
#'
#' @param object an object returned from \code{conditional_ccf} function
#' @param new_data data for which the time delay "dt" should be calculated
#' @param k_min minimum lag to be considered when calculating maximum cross-correlation.
#' default is $min(k)$ where k is the sequence of lags defined in \code{conditional_ccf}
#' @param k_max maximum lag to be considered when calculating maximum cross-correlation.
#' default is $max(k)$ where k is the sequence of lags defined in \code{conditional_ccf}
#'
#' @return an object of class "conditional_lag" with the following components
#'  \item{data_dt}{The original tibble passed to "new_data" appended with the estimated time lag "dt"
#'  and corresponding cross-correlation at maximum lag}
#'  \item{conditional_ccf_object}{The object passed to the \code{estimate_dt} function}
#'  \item{k_min}{minimum lag}
#'  \item{k_max}{maximum lag}
#'
#' @author Puwasala Gamakumara
#'
#' @export
estimate_dt <- function(object, new_data, k_min = NULL, k_max = NULL){


  k = object$other$k

  if(is.null(k_min)){
    k_min = min(k)
  }

  if(is.null(k_max)){
    k_max = max(k)
  }

  data_predict_ccf <- predict.conditional_ccf(object = object, new_data = new_data)

  DF_ccf_max <- data_predict_ccf %>%
    dplyr::select(.data$Timestamp, tidyselect::all_of(paste("k = ", k, sep = "")))
  colnames(DF_ccf_max) <- c("Timestamp", k)

  DF_ccf_max <- DF_ccf_max %>%
    tibble::as_tibble() %>%
    # tidyr::drop_na(.data$`1`) %>%
    tidyr::pivot_longer(cols = -.data$Timestamp, names_to = "dt", values_to = "ccf") %>%
    dplyr::filter(.data$dt %in% seq(k_min, k_max)) %>%
    dplyr::group_by(.data$Timestamp) %>%
    dplyr::slice(which.max(.data$ccf)) %>%
    dplyr::mutate(dt = as.numeric(.data$dt)) %>%
    dplyr::ungroup()

  new_data <- new_data %>%
    dplyr::left_join(DF_ccf_max, by = "Timestamp")

  # Add p value
  pval <-  data_predict_ccf %>%
    dplyr::select(.data$Timestamp, pvalue)
  new_data <- new_data %>%
    dplyr::left_join(pval, by = "Timestamp")

  return(structure(list(data_dt = new_data,
                        conditional_ccf_object = object,
                        k_min = k_min,
                        k_max = k_max),
                   class = "estimate_dt"))

}
