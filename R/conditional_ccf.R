#' Computing conditional cross-correlations at given lags
#'
#' This function computes cross correlation between $x_t$ and $y_{t+k}$ at $k = 1,2,...$
#' conditional on a set of time series $z_t$
#'
#' @param data a tibble containing all the time series including
#' ystar*xstar which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula A GAM formula. The response variable should be in the format of
#' I(x*y) ~ . See \code{\link[mgcv]{formula.gam}}.
#' @param lag_max Maximum lag at which to calculate the conditional ccf
#' @param fit_mean_x Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} for series x
#' @param fit_var_x Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_var}} for series x
#' @param fit_mean_y Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}} for series y
#' @param fit_var_y Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_var}} for series y
#' @param df_correlation a vector specifying the degrees of freedom to be considered for each numerical
#' predictor when fitting additive models for conditional cross-correlations. Each component of the
#' vector should corresponds to the degrees of freedom each predictor.
#' @return The function returns a list of objects of class
#' "glm" as described in \code{\link[stats]{glm}}. the length og the list is equal to lag_max
#'
#' @details{ Suppose $x_t$ and $y_t$ are conditionally normalised with respect
#' to $z_t$ using \code{conditional_mean} and \code{conditional_var}. Then
#' we can estimate the conditional cross-correlation between $x_t$ and $y_t$ at lag $k$, i.e. $r_k = E(x_ty_{t+k}|z_t)$
#' via generalised additive models (GAM). \code{conditional_ccf} uses natural splines implemented
#' in \code{splines} package to estimate the conditional cross-correlations between two
#' time series given a set of time series predictors. Users first need  to
#' normalise $x_t$ and $y_t$ at lag $k$ using \code{conditional_mean} and \code{conditional_var}}
#'
#' @seealso \code{\link[stats]{glm}}
#' @importFrom stats update glm as.formula
#' @importFrom purrr map
#'
#' @export
#'
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
#'       splines::ns(temperature_upstream, df = 5),
#'     lag_max = 10,
#'     fit_mean_x = fit_mean_x, fit_var_x = fit_var_x,
#'     fit_mean_y = fit_mean_y, fit_var_y = fit_var_y,
#'     df_correlation = c(5, 5)
#'   )
conditional_ccf <- function(data, formula, lag_max = 10, fit_mean_x,
                            fit_var_x, fit_mean_y, fit_var_y,
                            df_correlation) {
  vars <- all.vars(formula)
  x_name <- vars[1]
  y_name <- vars[2]

  # Calculate x_t*y_{t+k}
  tmp <- data |>
    dplyr::mutate(xstar = normalize(
      data, {{ x_name }}, fit_mean_x, fit_var_x
    ))
  tmp2 <- purrr::map_dfc(
      seq(lag_max), calc_xyk_star, tmp, {{ y_name }}, fit_mean_y,
      fit_var_y
    ) |>
    stats::setNames(paste("xystar_t", seq(lag_max), sep = ""))
  new_ts <- dplyr::bind_cols(data, tmp2)

  xynames <- colnames(new_ts)[grepl("xystar", names(new_ts))]
  corrl <- corrlink()

  # Fit GAM model for   x_t*y_{t+k}
  fit_ccf_gam <- function(k) {
    fk <- paste(xynames[k], "~.")
    formula_k <- stats::update(formula, stats::as.formula(fk))
    ccf_gam_fit_k <- stats::glm(
      formula = formula_k,
      data = new_ts,
      family = stats::gaussian(link = corrl),
      start = rep(0, (sum(df_correlation) + 1)),
      control = stats::glm.control(maxit = 400)
    )
    return(ccf_gam_fit_k)
  }

  ccf_gam_fit <- purrr::map(seq(lag_max), fit_ccf_gam)

  ccf_gam_fit$data <- new_ts
  ccf_gam_fit$lag_max <- lag_max
  class(ccf_gam_fit) <- c("conditional_ccf")
  return(ccf_gam_fit)
}

corrlink <- function() {
  ## link
  linkfun <- function(mu) {
    log((1 + mu) / (1 - mu))
  }
  ## inverse link
  linkinv <- function(eta) {
    (exp(eta) - 1) / (exp(eta) + 1)
  }
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) {
    2 * exp(eta) / (exp(eta) + 1)^2
  }
  valideta <- function(eta) TRUE
  link <- "corrlink"
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = valideta,
      name = link
    ),
    class = "link-glm"
  )
}

## Calculate x_t*y_{t+k}
calc_xyk_star <- function(k, data, y, fit_mean_y, fit_var_y) {
  old_ts_lead <- data |>
    dplyr::mutate_at({{ y }}, dplyr::lead, n = k)
  data$xstar * normalize(old_ts_lead,
    y = {{ y }}, fit_mean = fit_mean_y, fit_var = fit_var_y)
}
