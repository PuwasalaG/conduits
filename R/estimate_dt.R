#' Estimating time delay between two sensors in a river system
#'
#' This function estimates the time that takes water to flow from an upstream location to a downstream
#' location conditional on the observed water-quality variables from the upstream sensor. That time lag is
#' defined as the lag that gives maximum cross-correlation conditional on upstream water-quality variables.
#'
#' @param x 	Model object of class "conditional_ccf" returned from
#'  \code{\link[conduits]{conditional_ccf}}
#' @return A \code{\link[tibble]{tibble}} with estimated time lag "dt"
#'  and corresponding maximum cross-correlation
#' @importFrom dplyr mutate
#' @importFrom stats predict.glm pnorm
#' @importFrom purrr map_dfc
#'
#' @author Puwasala Gamakumara & Priyanga Dilini Talagala
#'
#' @examples
#'
#' old_ts <- NEON_PRIN_5min_cleaned |>
#'   dplyr::select(
#'     Timestamp, site, turbidity, level, temperature
#'   ) |>
#'   tidyr::pivot_wider(
#'     names_from = site,
#'     values_from = turbidity:temperature
#'   )
#'
#' fit_mean_y <- old_ts |>
#'   conditional_mean(turbidity_downstream ~
#'     s(level_upstream, k = 5) +
#'     s(temperature_upstream, k = 5))
#'
#' fit_var_y <- old_ts |>
#'   conditional_var(
#'     turbidity_downstream ~
#'       s(level_upstream, k = 4) +
#'       s(temperature_upstream, k = 4),
#'     family = "Gamma",
#'     fit_mean = fit_mean_y
#'   )
#'
#' fit_mean_x <- old_ts |>
#'   conditional_mean(turbidity_upstream ~
#'     s(level_upstream, k = 5) +
#'     s(temperature_upstream, k = 5))
#'
#' fit_var_x <- old_ts |>
#'   conditional_var(
#'     turbidity_upstream ~
#'       s(level_upstream, k = 4) +
#'       s(temperature_upstream, k = 4),
#'     family = "Gamma",
#'     fit_mean = fit_mean_x
#'   )
#'
#' fit_c_ccf <- old_ts |>
#'   tidyr::drop_na() |>
#'   conditional_ccf(
#'     I(turbidity_upstream * turbidity_downstream) ~
#'       splines::ns(level_upstream, df = 3) +
#'       splines::ns(temperature_upstream, df = 3),
#'     lag_max = 10,
#'     fit_mean_x = fit_mean_x, fit_var_x = fit_var_x,
#'     fit_mean_y = fit_mean_y, fit_var_y = fit_var_y,
#'     df_correlation = c(3, 3)
#'   )
#'
#' new_data <- fit_c_ccf |> estimate_dt()
#' @export
#'
estimate_dt <- function(x) {
  data_inf <- x |> augment()
  k <- x$lag_max
  cnames <- paste("c", seq(k), sep = "")
  data_sub <- data_inf
  data_sub$dt <- max.col(data_sub[cnames], ties.method = "first")
  data_sub$max_ccf <- apply(data_sub[cnames], 1, max)

  # p value calculation
  predict_ccf_gam_stdz <- function(t) {
    cond_ccf <- stats::predict.glm(x[[t]],
      newdata = x$data,
      type = "response",
      se.fit = TRUE
    )
    cond_ccf_std <- cond_ccf$fit / cond_ccf$se.fit
    return(cond_ccf_std)
  }

  c_ccf_est_std <- purrr::map_dfc(seq(k), predict_ccf_gam_stdz)
  row_max <- apply(c_ccf_est_std, 1, max)
  data_sub$pval <- 1 - (stats::pnorm(row_max))^k

  return(data_sub)
}

#' Computing bootstrapped confidence intervals for dt
#'
#' This function computes the  bootstrapped confidence intervals for dt.
#' It resample the residuals from the various models used in the
#' conditional cross-correlation calculation to generate new data.
#' As the  residuals are serially correlated,  a sieve bootstrap
#' approach to capture the autocorrelation  structure in the data.
#'
#' @param x 	Model object of class "conditional_ccf" returned from
#'  \code{\link[conduits]{conditional_ccf}}
#' @param m number of replications for boostrap confidence intervals
#' @param new_data the dataset with the some predictors that are
#'  set to the median value (if required). Default is set to NULL.
#' @return A \code{\link[tibble]{tibble}} with estimated time lag "dt"
#'
#' @importFrom broom augment_columns
#' @importFrom splines ns
#' @importFrom purrr map
#' @importFrom forecast auto.arima
#' @importFrom stats quantile
#'
#' @author Priyanga Dilini Talagala & Puwasala Gamakumara
#'
#' @examples
#' \dontrun{
#' old_ts <- NEON_PRIN_5min_cleaned |>
#'   dplyr::select(
#'     Timestamp, site, turbidity, level, temperature
#'   ) |>
#'   tidyr::pivot_wider(
#'     names_from = site,
#'     values_from = turbidity:temperature
#'   )
#' fit_mean_y <- old_ts |>
#'   conditional_mean(turbidity_downstream ~
#'     s(level_upstream, k = 5) +
#'     s(temperature_upstream, k = 5)
#'   )
#' fit_var_y <- old_ts |>
#'   conditional_var(
#'     turbidity_downstream ~
#'       s(level_upstream, k = 4) +
#'       s(temperature_upstream, k = 4),
#'     family = "Gamma",
#'     fit_mean = fit_mean_y
#'   )
#' fit_mean_x <- old_ts |>
#'   conditional_mean(turbidity_upstream ~
#'     s(level_upstream, k = 5) +
#'     s(temperature_upstream, k = 5)
#'   )
#' fit_var_x <- old_ts |>
#'   conditional_var(
#'     turbidity_upstream ~
#'       s(level_upstream, k = 4) +
#'       s(temperature_upstream, k = 4),
#'     family = "Gamma",
#'     fit_mean = fit_mean_x
#'   )
#' fit_c_ccf <- old_ts |>
#'   tidyr::drop_na() |>
#'   conditional_ccf(
#'     I(turbidity_upstream * turbidity_downstream) ~
#'       splines::ns(level_upstream, df = 3) +
#'       splines::ns(temperature_upstream, df = 3),
#'     lag_max = 10,
#'     fit_mean_x = fit_mean_x, fit_var_x = fit_var_x,
#'     fit_mean_y = fit_mean_y, fit_var_y = fit_var_y,
#'     df_correlation = c(3, 3)
#'   )
#' df_dt <- fit_c_ccf |> calc_dt_CI(100)
#'
#' # Calculate  dt vs an  upstream covariate while holding the
#' # remaining upstream covariates at their medians
#' new_data <- fit_c_ccf$data
#' new_data <- new_data |>
#'   dplyr::mutate(temperature_upstream = median(temperature_upstream))
#' df_dt2 <- fit_c_ccf |> calc_dt_CI(100, new_data)
#' }
#'
#'
#' @export
#'
calc_dt_CI <- function(x, m, new_data = NULL) {
  lag_max <- x$lag_max
  corrl <- corrlink()

  calc_dt_CI_m <- function(i) {
    # Add fitted values, residuals, and other common
    # outputs to an augment call for each k
    aug_ccf_gam <- function(k) {
      fit <- x[[k]]
      aug <- broom::augment_columns(fit, fit$data)
      return(aug)
    }
    df_ccf_aug <- purrr::map(seq(lag_max), aug_ccf_gam)

    # Extract residuals of the fitted GAMs models for each k
    res_ccf_gam <- function(k) {
      aug <- df_ccf_aug[[k]]
      res_cond_ccf <- aug$.resid
      return(res_cond_ccf)
    }
    df_ccf_resid <- purrr::map(seq(lag_max), res_ccf_gam)

    # Fit kth order autoregressive model for the
    # residual series for each k
    arma_fit <- function(k) {
      res <- df_ccf_resid[[k]]
      fit <- forecast::auto.arima(res,
        approximation = FALSE,
        stepwise = FALSE,
        max.q = 0, d = 0,
        max.order = 20
      )
      res <- stats::residuals(fit) |> as.vector()
      fitted <- stats::fitted(fit) |> as.vector()
      df <- cbind(res, fitted)
      return(df)
    }
    df_ar_aug <- purrr::map(seq(lag_max), arma_fit)

    # Generate bootstrap samples from the residuals
    # of the fitted AR models
    bootstrap_resid <- function(k) {
      ar_aug <- df_ar_aug[[k]] |> as_tibble()
      res <- ar_aug$res
      bs_resid <- sample(res,
        size = length(res),
        replace = TRUE
      )

      return(bs_resid)
    }
    df_bs_resid_k <- purrr::map(seq(lag_max), bootstrap_resid)

    # Create a bootstrap residual series of the AR models
    # fitted  for each k
    bs_ar_res <- function(k) {
      res <- df_ar_aug[[k]][, "fitted"] + df_bs_resid_k[[k]]
      return(res)
    }
    df_bs_ar_k <- purrr::map(seq(lag_max), bs_ar_res)

    # Create a bootstrap response series for each k
    bs_yx_k <- function(k) {
      yx_bs <- df_ccf_aug[[k]]$.fitted + df_bs_ar_k[[k]]
      return(yx_bs)
    }
    df_bs_yx_k <- purrr::map(seq(lag_max), bs_yx_k)

    cbind_bs <- function(k) {
      xy_k <- df_bs_yx_k[[k]]
      df <- cbind(df_ccf_aug[[k]], xy_k)
      return(df)
    }
    df_bs <- purrr::map(seq(lag_max), cbind_bs)

    # Fit a GAM to the bootstrapped data for each k
    fit_GAM_bs <- function(k) {
      df <- df_bs[[k]]
      formula <- x[[k]]$formula
      fk <- paste("xy_k", "~.")
      f <- stats::update(formula, stats::as.formula(fk))
      ccf_gam_fit <- stats::glm(
        formula = f,
        data = df,
        family = stats::gaussian(link = corrl),
        start = rep(0, (x[[k]]$rank)),
        control = stats::glm.control(maxit = 400)
      )
      return(ccf_gam_fit)
    }
    bs_gam_fit <- purrr::map(seq(lag_max), fit_GAM_bs)

    predict_ccf_gam <- function(k) {
      if (is.null(new_data)) {
        new_data <- bs_gam_fit[[k]]$data
      }

      cond_ccf <- stats::predict.glm(bs_gam_fit[[k]],
        newdata = new_data,
        type = "response"
      )
      df <- cbind(new_data["Timestamp"], cond_ccf)
      return(df)
    }
    cond_ccf_est <- purrr::map(seq(lag_max), predict_ccf_gam)

    # Compute time delay, dt, for the boostrap samples
    cnames <- paste("c", seq(lag_max), sep = "")
    df <- purrr::reduce(cond_ccf_est, dplyr::left_join, by = "Timestamp")
    #df <- plyr::join_all(cond_ccf_est, by = "Timestamp", type = "left")
    colnames(df) <- c("Timestamp", cnames)
    data_sub <- df |>
      dplyr::mutate(
        dt =
          max.col(.[cnames],
            ties.method = "first"
          )
      ) |>
      select(Timestamp, dt)

    cat("Iteration : ", i, "\n")
    return(data_sub)
  }

  dt_bs <- purrr::map(seq(m), calc_dt_CI_m)
  df_dt_bs <- purrr::reduce(dt_bs, dplyr::left_join, by = "Timestamp")
  #df_dt_bs <- plyr::join_all(dt_bs, by = "Timestamp", type = "left")
  colnames(df_dt_bs) <- c("Timestamp", paste("dt", seq(m), sep = ""))

  df_bootstrap_CI <- df_dt_bs |>
    dplyr::as_tibble() |>
    tidyr::pivot_longer(cols = starts_with("dt")) |>
    dplyr::group_by(Timestamp) |>
    dplyr::summarise(
      dt_95_LI = stats::quantile(value, probs = c(0.025), na.rm = TRUE),
      dt_95_UI = stats::quantile(value, probs = c(0.975), na.rm = TRUE),
      dt_80_LI = stats::quantile(value, probs = c(0.1), na.rm = TRUE),
      dt_80_UI = stats::quantile(value, probs = c(0.9), na.rm = TRUE)
    )

  # data_with_dt<- x |> estimate_dt()

  # names <- colnames(data_with_dt)[!grepl("xystar", names(data_with_dt)) &
  #                                !grepl("c", names(data_with_dt) ) ]

  df <- df_dt_bs |>
    #  dplyr::select(names) |>
    dplyr::left_join(df_bootstrap_CI, by = "Timestamp")

  return(df)
}
