globalVariables(c(".fitted", ".se.fit", ".cond_m", "."))

#' Augment data with information from a conditional mean fit or conditional variance fit
#'
#' This function produces partial residuals for each predictor,
#' and the estimated conditional means, standard error and confidence limits.
#'
#' @param x 	Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} or  \code{\link[conduits]{conditional_var}}
#'  with information to append to observations.
#' @param level Confidence level. Default is set to 0.95.
#' @param ...	 Additional arguments, not currently used
#' @return A \code{\link[tibble]{tibble}} with information
#'  about data points.
#'
#' @seealso \code{\link[mgcv]{gam}}
#' @importFrom dplyr mutate bind_cols rename
#' @importFrom stats predict qnorm
#' @importFrom broom augment
#' @examples
#' data <- NEON_PRIN_5min_cleaned |>
#'   dplyr::filter(site == "upstream") |>
#'   dplyr::select(Timestamp, turbidity, level, conductance, temperature)
#'
#' fit_mean <- data |>
#'   conditional_mean(turbidity ~ s(level, k = 8) +
#'     s(conductance, k = 8) + s(temperature, k = 8))
#'
#' data_inf <- fit_mean |> augment()
#' @export augment.conditional_moment
#' @export
#'
augment.conditional_moment <- function(x, level = 0.95, ...) {
  # getting each component of the linear predictor from the fitted model
  fv <- stats::predict(x, type = "terms")
  # aug <- broom:::augment.gam(x)
  aug <- augment.gam(x)

  if (x$type == "conditional_var") {
    aug <- aug |>
      dplyr::mutate(.fitted = mgcv::predict.gam(x,
        newdata = aug,
        type = "response"
      ))
  }

  .partial.res <- fv + aug$.resid
  colnames(.partial.res) <- paste0(".presid_", colnames(.partial.res))

  zq <- abs(stats::qnorm((1 - level) / 2))
  data <- aug |>
    dplyr::mutate(
      .cond_m = .fitted,
      .LI = as.numeric(.fitted - zq * .se.fit),
      .UI = as.numeric(.fitted + zq * .se.fit),
    ) |>
    dplyr::bind_cols(.partial.res)

  if (x$type == "conditional_mean") {
    data <- data |> dplyr::rename(.cond_EX = .cond_m)
  }
  if (x$type == "conditional_var") {
    data <- data |> dplyr::rename(.cond_VAR = .cond_m)
  }

  return(data)
}

#' Augment data with information from a conditional cross-correlation fit
#'
#' This function produces estimated conditional cross-correlation between
#' $x_t$ and $y_t$ at lag $k$, i.e. $r_k = E(x_ty_{t+k}|z_t)$.
#'
#' @param x 	Model object of class "conditional_ccf" returned from
#'  \code{\link[conduits]{conditional_ccf}} with information to append to observations.
#' @param ...	 Additional arguments, not currently used.
#' @return A \code{\link[tibble]{tibble}} with information
#' about data points.
#' @export augment.conditional_ccf
#' @export
#' @importFrom purrr map_dfc
#' @examples
#'
#' old_ts <- NEON_PRIN_5min_cleaned |>
#'   dplyr::select(
#'     Timestamp, site, turbidity, level,
#'     conductance, temperature
#'   ) |>
#'   tidyr::pivot_wider(
#'     names_from = site,
#'     values_from = turbidity:temperature
#'   )
#'
#' fit_mean_y <- old_ts |>
#'   conditional_mean(turbidity_downstream ~
#'     s(level_upstream, k = 8) +
#'     s(conductance_upstream, k = 8) +
#'     s(temperature_upstream, k = 8))
#'
#' fit_var_y <- old_ts |>
#'   conditional_var(
#'     turbidity_downstream ~
#'       s(level_upstream, k = 7) +
#'       s(conductance_upstream, k = 7) +
#'       s(temperature_upstream, k = 7),
#'     family = "Gamma",
#'     fit_mean = fit_mean_y
#'   )
#'
#' fit_mean_x <- old_ts |>
#'   conditional_mean(turbidity_upstream ~
#'     s(level_upstream, k = 8) +
#'     s(conductance_upstream, k = 8) +
#'     s(temperature_upstream, k = 8))
#'
#' fit_var_x <- old_ts |>
#'   conditional_var(
#'     turbidity_upstream ~
#'       s(level_upstream, k = 7) +
#'       s(conductance_upstream, k = 7) +
#'       s(temperature_upstream, k = 7),
#'     family = "Gamma",
#'     fit_mean = fit_mean_x
#'   )
#'
#' fit_c_ccf <- old_ts |>
#'   tidyr::drop_na() |>
#'   conditional_ccf(
#'     I(turbidity_upstream * turbidity_downstream) ~ splines::ns(
#'       level_upstream,
#'       df = 5
#'     ) +
#'       splines::ns(conductance_upstream, df = 5),
#'     lag_max = 10,
#'     fit_mean_x = fit_mean_x, fit_var_x = fit_var_x,
#'     fit_mean_y = fit_mean_y, fit_var_y = fit_var_y,
#'     df_correlation = c(5, 5)
#'   )
#'
#' data_inf <- fit_c_ccf |> augment()
augment.conditional_ccf <- function(x, ...) {
  lag_max <- length(x) - 2
  data_NEW <- x$data

  predict_ccf_gam <- function(k) {
    cond_ccf <- stats::predict.glm(x[[k]],
      newdata = data_NEW,
      type = "response"
    )
    return(cond_ccf)
  }

  cond_ccf_est <- purrr::map_dfc(seq(lag_max), predict_ccf_gam) |>
    stats::setNames(paste("c", seq(lag_max), sep = ""))

  return(dplyr::bind_cols(data_NEW, cond_ccf_est))
}

#' Augment data with information from a conditional auto-correlation fit
#'
#' This function produces estimated conditional autocorrelation between
#' $x_t$ and $y_t$ at lag $k$, i.e. $r_k = E(x_ty_{t+k}|z_t)$.
#'
#' @param x 	Model object of class "conditional_acf" returned from
#'  \code{\link[conduits]{conditional_acf}} with information to append to observations.
#' @param ...	 Additional arguments, not currently used.
#' @return A \code{\link[tibble]{tibble}} with information
#' about data points.
#' @export augment.conditional_acf
#' @export
#' @importFrom purrr map_dfc
#' @importFrom scales rescale
#' @examples
#' old_ts <- NEON_PRIN_5min_cleaned |>
#'   dplyr::select(
#'     Timestamp, site, turbidity, level,
#'     conductance, temperature
#'   ) |>
#'   tidyr::pivot_wider(
#'     names_from = site,
#'     values_from = turbidity:temperature
#'   )
#'
#' fit_mean <- old_ts |>
#'   conditional_mean(turbidity_downstream ~
#'     s(level_upstream, k = 8) +
#'     s(conductance_upstream, k = 8) +
#'     s(temperature_upstream, k = 8))
#'
#' fit_var <- old_ts |>
#'   conditional_var(
#'     turbidity_downstream ~
#'       s(level_upstream, k = 7) +
#'       s(conductance_upstream, k = 7) +
#'       s(temperature_upstream, k = 7),
#'     family = "Gamma",
#'     fit_mean = fit_mean
#'   )
#' fit_c_acf <- old_ts |>
#'   tidyr::drop_na() |>
#'   conditional_acf(
#'     turbidity_upstream ~ splines::ns(level_upstream, df = 5) +
#'       splines::ns(conductance_upstream, df = 5),
#'     lag_max = 10, fit_mean = fit_mean, fit_var = fit_var,
#'     df_correlation = c(5, 5)
#'   )
#'
#' data_inf <- fit_c_acf |> augment()
augment.conditional_acf <- function(x, ...) {
  lag_max <- length(x) - 1
  data_NEW <- x$data

  predict_acf_gam <- function(k) {
    cond_acf <- stats::predict.glm(x[[k]], newdata = data_NEW, type = "response")
    cond_acf <- scales::rescale(cond_acf, to = c(-1, 1))
    return(cond_acf)
  }

  cond_acf_est <- purrr::map_dfc(seq(lag_max), predict_acf_gam) |>
    stats::setNames(paste("r", seq(lag_max), sep = ""))

  return(dplyr::bind_cols(data_NEW, cond_acf_est))
}

augment.gam <- utils::getFromNamespace("augment.gam", "broom")
