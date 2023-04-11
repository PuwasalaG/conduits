#' Normalize a series using conditional moments
#'
#' This function produces  a normalized series using conditional moments.
#' @param y The variable name
#' @param fit_mean 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} with information to append to observations.
#' @param fit_var 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_var}} with information to append to observations.
#' @param data a tsibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @return A vector of conditional normliased series
#' @importFrom dplyr ensym pull mutate
#' @importFrom mgcv predict.gam
#' @importFrom tsibble as_tsibble index
#' @examples
#' data <- NEON_PRIN_5min_cleaned |>
#'   dplyr::filter(site == "upstream") |>
#'   dplyr::select(Timestamp, turbidity, level, conductance, temperature) |>
#'   tsibble::as_tsibble(index = Timestamp)
#'
#' fit_mean <- data |>
#'   conditional_mean(turbidity ~ s(level, k = 8) +
#'     s(conductance, k = 8) + s(temperature, k = 8))
#'
#' fit_var <- data |>
#'   conditional_var(
#'     turbidity ~ s(level, k = 7) + s(conductance, k = 7) + s(temperature, k = 7),
#'     family = "Gamma",
#'     fit_mean = fit_mean
#'   )
#'
#' new_ts <- data |>
#'   dplyr::mutate(ystar = conduits::normalize(data, turbidity, fit_mean, fit_var))
#'
#' @export
#'
normalize <- function(data, y, fit_mean, fit_var) {
  y <- dplyr::ensym(y)
  cond_EY <- as.numeric(mgcv::predict.gam(fit_mean, newdata = data, type = "response"))
  cond_VY <- as.numeric(mgcv::predict.gam(fit_var,
    newdata = data, type = "response"
  ))
  y_star <- ((data |> dplyr::pull({{ y }})) - cond_EY) / sqrt(cond_VY)
  # y_norm <- data |>
  #  dplyr::mutate(y_star) |>
  #  tsibble::as_tsibble(index = tsibble::index(data))
  return(y_star)
}



#' Unnormalize a series using conditional moments
#'
#' This function produces  an unnormalized series using conditional moments.
#' @param ystar The normalized variable name
#' @param fit_mean 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} with information to append to observations.
#' @param fit_var 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_var}} with information to append to observations.
#' @param data a tsibble containing all the time series
#' which are uniquely identified by the corresponding
#' Timestamp.
#' @return A \code{\link[tsibble]{tsibble}} with the conditional normliased series
#' @importFrom dplyr ensym pull mutate
#' @importFrom mgcv predict.gam
#' @importFrom tsibble as_tsibble index
#' @examples
#' data <- NEON_PRIN_5min_cleaned |>
#'   dplyr::filter(site == "upstream") |>
#'   dplyr::select(Timestamp, turbidity, level, conductance, temperature) |>
#'   tsibble::as_tsibble(index = Timestamp)
#'
#' fit_mean <- data |>
#'   conditional_mean(turbidity ~ s(level, k = 8) +
#'     s(conductance, k = 8) + s(temperature, k = 8))
#'
#' fit_var <- data |>
#'   conditional_var(
#'     turbidity ~ s(level, k = 7) + s(conductance, k = 7) + s(temperature, k = 7),
#'     family = "Gamma",
#'     fit_mean = fit_mean
#'   )
#'
#' new_ts <- data |>
#'   dplyr::mutate(ystar = normalize(data, turbidity, fit_mean, fit_var))
#'
#' # For demonstrative purposes, declare three data points
#' # as missing values.
#' new_ts[3:5, 6] <- NA
#'
#' \dontrun{
#' library(fable)
#' library(dplyr)
#' impute_ts <- new_ts |>
#'   model(ARIMA(ystar)) |>
#'   interpolate(new_ts) |>
#'   rename(y_star_impt = ystar) |>
#'   full_join(new_ts, by = "Timestamp")
#' impute_ts <- impute_ts
#'   mutate(y = unnormalize(impute_ts, y_star_impt, fit_mean, fit_var))
#' }
#'
#' @export
#'
unnormalize <- function(data, ystar, fit_mean, fit_var) {
  ystar <- dplyr::ensym(ystar)
  cond_EY <- as.numeric(mgcv::predict.gam(fit_mean, newdata = data))
  cond_VY <- as.numeric(mgcv::predict.gam(fit_var,
    newdata = data,
    type = "response"
  ))
  y <- (data |> dplyr::pull({{ ystar }})) * sqrt(cond_VY) + cond_EY
  # y_trns <- data |>
  #  dplyr::mutate(y) |>
  #  tsibble::as_tsibble(index = tsibble::index(data))
  return(y)
}
