#' Visualising the models fitted for conditional moments
#'
#' A function for visualising the fitted generalised additive models for
#' conditional means and variances from \code{conditional_moments} function
#'
#' @param object an object returned from \code{conditional_moments} function
#' @param type the type of moment want to visualise. This can take one of
#' "mean" or "variance"
#' @param show_points Should partial residuals be shown around the fitted line?
#' @param ... ignored
#'
#' @return returns plots visualising fitted gam model for conditional
#' mean or variance vs each predictor.
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_ribbon ylab geom_rug
#' @importFrom tidyselect eval_select
#' @importFrom dplyr enquo
#'
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[visreg]{visreg}}
#'
#' @author Puwasala Gamakumara
#'
#' @export
autoplot.conditional_moments <- function(
  object, type = c("mean", "variance"), show_points = TRUE,
  ...) {
  type <- match.arg(type)
  x <- object$x
  z_numeric <- object$z_numeric
  z_factors <- object$z_factors

  names_x <- names(eval_select(enquo(x), object$data_conditional_moments))
  names_z_numeric <- names(eval_select(enquo(z_numeric), object$data_conditional_moments))
  names_z_factors <- names(eval_select(enquo(z_factors), object$data_conditional_moments))

  p <- length(names_z_numeric) + length(names_z_factors)
  ncol <- max(min(p, 3), 2)

  if (type == "mean") {
    data <- object$mean_gam$data
    fit <- object$mean_gam$fit
  } else if (type == "variance") {
    data <- object$var_gam$data
    fit <- object$var_gam$fit
  }

  data_vis <- data %>%
    dplyr::select(.data$Timestamp, {{ x }}, .data$X_Ex2, {{ z_numeric }}, {{ z_factors }}) %>%
    tidyr::drop_na()
  summary_fit <- fit %>% summary()
  fit_edf <- summary_fit$edf %>% round(digits = 3)
  # getting each component of the linear predictor from the fitted model
  fv <- mgcv::predict.gam(fit, type = "terms")

  if (is.null(z_factors)) {
    colnames(fv) <- paste("s_", names_z_numeric, sep = "")
  } else {
    colnames(fv) <- c(
      paste("s_", names_z_factors, sep = ""),
      paste("s_", names_z_numeric, sep = "")
    )
  }

  # residuals
  resid <- mgcv::residuals.gam(fit, type = "working")
  if(type == "var")
    resid <- exp(resid)

  # combining linear predictor and residuals
  resid <- cbind(fv, resid) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Timestamp = data_vis$Timestamp)
  ## -- visualising fitted model
  plot_fit <- list()

  for (i in seq_along(names_z_numeric)) {
    ## -- preparing data for visualisation for numerical variables
    data_vis_z <- data_vis %>%
      dplyr::mutate_at(
        names_z_numeric[-i],
        ~ median(., na.rm = T)
      )
    if (!is.null(z_factors)) {
      data_vis_z <- data_vis_z %>%
        dplyr::mutate_at(names_z_factors, ~ names(which.max(table(.))))
    }

    # predicted values
    predict_fit <- mgcv::predict.gam(fit,
        type = "link", newdata = data_vis_z, se.fit = TRUE
      ) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        LI = as.numeric(.data$fit - 1.96*.data$se.fit),
        UI = as.numeric(.data$fit + 1.96*.data$se.fit),
        Timestamp = data_vis_z$Timestamp
      )

    s_name_zi <- paste("s_", names_z_numeric[i], sep = "")
    data_vis_z <- list(data_vis_z, resid, predict_fit) %>%
      purrr::reduce(dplyr::left_join, by = "Timestamp") %>%
      dplyr::mutate(partial_resid = resid + !!rlang::sym(s_name_zi))

    name_zi <- names_z_numeric[i]
    edf_zi <- fit_edf[i]
    range_zi <- data_vis_z %>% dplyr::pull(!!rlang::sym(name_zi))

    plot_fit[[i]] <- data_vis_z %>%
      ggplot()
    if(show_points) {
      plot_fit[[i]] <- plot_fit[[i]] +
        geom_point(
          aes(x = !!rlang::sym(name_zi), y = .data$partial_resid),
          color = "#666666", size = 0.5
        )
    }
    plot_fit[[i]] <- plot_fit[[i]] +
      geom_line(
        color = "#0099FF",
        aes(x = !!rlang::sym(name_zi), y = (.data$fit)),
        size = 1
      ) +
      geom_ribbon(
        aes(x = !!rlang::sym(name_zi), y = (.data$fit),
            ymin = .data$LI, ymax = .data$UI),
        alpha = 0.5,
        fill = "grey70"
      ) +
      ylab(paste("f(", name_zi, ",", edf_zi, ")", sep = "")) +
      scale_x_continuous(breaks = scales::pretty_breaks()(range_zi))
  }

  invisible(plot_fit)
}
