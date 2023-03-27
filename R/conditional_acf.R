#' Computing conditional autocorrelations at given lags
#'
#' This function computes autocorrelation between $x_t$ and $y_{t+k}$ at $k = 1,2,...$
#' conditional on a set of time series $z_t$
#'
#' @param data a tibble containing all the time series including
#' $ystar*ystar_{t-k}$ which are uniquely identified by the corresponding
#' Timestamp.
#' @param formula A GAM formula. See \code{\link[mgcv]{formula.gam}}.
#' @param lag_max Maximum lag at which to calculate the conditional acf
#' @param fit_mean Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_mean}}
#' @param fit_var Model object of class "conditional_moment" returned from
#'  \code{\link[conduits]{conditional_var}}
#' @param df_correlation a vector specifying the degrees of freedom to be considered for each numerical
#' predictor when fitting additive models for conditional auto-correlations. Each component of the
#' vector should corresponds to each predictor specified in "z_numeric".
#' @return The function returns a list of objects of class
#' "glm" as described in \code{\link[stats]{glm}}.
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
#' @import dplyr
#' @importFrom tidyr drop_na
#'
#' @examples
#'
#' old_ts <- NEON_PRIN_5min_cleaned %>%
#'   dplyr::select(
#'     Timestamp, site, turbidity, level,
#'     conductance, temperature
#'   ) %>%
#'   tidyr::pivot_wider(
#'     names_from = site,
#'     values_from = turbidity:temperature
#'   )
#'
#' fit_mean <- old_ts %>%
#'   conditional_mean(turbidity_downstream ~
#'   s(level_upstream, k = 8) +
#'     s(conductance_upstream, k = 8) +
#'     s(temperature_upstream, k = 8))
#'
#' fit_var <- old_ts %>%
#'   conditional_var(
#'     turbidity_downstream ~
#'     s(level_upstream, k = 7) +
#'       s(conductance_upstream, k = 7) +
#'       s(temperature_upstream, k = 7),
#'     family = "Gamma",
#'     fit_mean
#'   )
#'
#'
#' fit_c_acf <- old_ts %>%
#'   tidyr::drop_na() %>%
#'   conditional_acf(
#'     turbidity_upstream ~ splines::ns(level_upstream, df = 5) +
#'       splines::ns(conductance_upstream, df = 5),
#'     lag_max = 10, fit_mean, fit_var,
#'     df_correlation = c(5, 5)
#'   )
#' @export
#'
conditional_acf <- function(data, formula, lag_max, fit_mean, fit_var, df_correlation) {
  vars <- all.vars(formula)
  y_name <- vars[1]

  # Calculate y_t*y_{t+k}
  new_ts <- data %>%
    dplyr::mutate(ystar = normalize(
      ., {{ y_name }}, fit_mean, fit_var
    )) %>%
    purrr::map_dfc(
      1:lag_max, calc_yyk_star, ., {{ y_name }},
      fit_mean, fit_var
    ) %>%
    stats::setNames(paste("ystarystar", 1:lag_max, sep = "")) %>%
    dplyr::bind_cols(data, .)

  yynames <- colnames(new_ts)[grepl("ystarystar", names(new_ts))]
  corrl <- corrlink()

  fit_acf_gam <- function(k) {
    fk <- paste(yynames[k], "~.")
    formula_k <- stats::update(formula, stats::as.formula(fk))
    new_ts_ystar <- new_ts %>%
      dplyr::select(yynames[k], vars[-1]) %>%
      tidyr::drop_na()
    acf_gam_fit_k <- stats::glm(
      formula = formula_k,
      data = new_ts_ystar,
      # family = stats::gaussian(link = corrl),
      family = stats::gaussian(), # Check the link function
      start = rep(0, (sum(df_correlation) + 1)),
      control = stats::glm.control(maxit = 400)
    )
    return(acf_gam_fit_k)
  }

  acf_gam_fit <- purrr::map(1:lag_max, fit_acf_gam)
  acf_gam_fit$data <- new_ts
  class(acf_gam_fit) <- c("conditional_acf")
  return(acf_gam_fit)
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
  structure(list(
    linkfun = linkfun,
    linkinv = linkinv,
    mu.eta = mu.eta,
    valideta = valideta,
    name = link
  ),
  class = "link-glm"
  )
}

calc_yyk_star <- function(k, data, y, fit_mean, fit_var) {
  old_ts_lead <- data %>%
    dplyr::mutate_at({{ y }}, dplyr::lag, n = k) %>%
    normalize(., {{ y }}, fit_mean, fit_var) * data$ystar
}
