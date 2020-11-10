#' Visualising the models fitted for conditional means
#'
#' A function for visualising the fitted generalised additive models for
#' conditional means and variances from \code{conditional_moments} function
#'
#' @param object an object returned from \code{conditional_moments} function
#' @param type the type of moment want to visualise. This can take one of
#' "mean" or "variance"
#' @param ... further arguments passed to or from other methods
#'
#' @return returns plots visualising fitted gam model for conditional
#' mean or variance vs each predictor.
#'
#' @examples
#'
#' @export
autoplot.conditional_moments <- function(object,
                                         type = c("mean", "variance"),
                                         ...){
  type <- match.arg(type)
  x <- object$x
  z_numeric <- object$z_numeric
  z_factors <- object$z_factors

  mean_fit <- object$mean_gam$fit
  var_fit <- object$var_gam$fit

  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), object$data_conditional_moments))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), object$data_conditional_moments))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), object$data_conditional_moments))

  p <- length(names_z_numeric) + length(names_z_factors)

  if(p <= 2){
    n_col <- 2
  }

  if(p > 2){
    n_col <- 3
  }



  if(type == "mean"){

    data_vis <- object$mean_gam$data %>%
      dplyr::select(.data$Timestamp, {{x}}, {{z_numeric}}, {{z_factors}}) %>%
      tidyr::drop_na()

    # getting each component of the linear predictor from the fitted model
    fv_mean <- mgcv::predict.gam(mean_fit, type = "terms")

    if(is.null(z_factors)){
      colnames(fv_mean) <- paste("s_", names_z_numeric, sep = "")
    }else{
      colnames(fv_mean) <- c(paste("s_", names_z_factors, sep = ""),
                             paste("s_", names_z_numeric, sep = ""))
    }

    # residuals
    resid_mean <- mgcv::residuals.gam(mean_fit, type = "working")

    # combining linear predictor and residuals
    resid_mean <- cbind(fv_mean, resid_mean) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Timestamp = data_vis$Timestamp)

    ##-- visualising fitted mean model
    plot_means <- list()
    data <- list()

    for (i in 1:length(names_z_numeric)) {
      ##-- preparing data for visualisation for numerical variables
      data_vis_z <- data_vis %>%
        dplyr::mutate_at(names_z_numeric[-i],
                  ~median(.,na.rm = T))
      if(!is.null(z_factors)){
        data_vis_z <- data_vis_z %>%
          dplyr::mutate_at(names_z_factors, ~names(which.max(table(.))))
      }

      # predicted values
      predict_fit_mean <- mgcv::predict.gam(mean_fit, newdata = data_vis_z,
                                      se.fit = T) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(LI = as.numeric(.data$fit - .data$se.fit),
               UI = as.numeric(.data$fit + .data$se.fit),
               Timestamp = data_vis_z$Timestamp)

      s_name_zi <- paste("s_", names_z_numeric[i],
                         sep = "")
      data_vis_z <- list(data_vis_z, resid_mean, predict_fit_mean) %>%
        purrr::reduce(dplyr::left_join, by = "Timestamp") %>%
        dplyr::mutate(partial_resid = resid_mean + !!rlang::sym(s_name_zi))

      name_zi <- names_z_numeric[i]
      plot_means[[i]] <- data_vis_z %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = !!rlang::sym(name_zi),
                                         y = .data$partial_resid),
                            color = "#666666", size = 0.5) +
        ggplot2::geom_line(color = "#0099FF",
                           ggplot2::aes(x = !!rlang::sym(name_zi),
                                        y = .data$fit),
                           size = 1) +
        ggplot2::geom_ribbon(ggplot2::aes(x = !!rlang::sym(name_zi),
                                          y = .data$fit,
                                          ymin = .data$LI,
                                          ymax = .data$UI),
                    alpha = 0.5,
                    fill = "grey70") +
        ggplot2::ylab(paste("f(",name_zi,")", sep = ""))

      data[[i]] <- data_vis_z %>%
        dplyr::select(.data$Timestamp)

    }

    gridExtra::grid.arrange(grobs = plot_means,
                                           ncol = n_col,
                                           top = paste("conditional mean of",
                                                       names_x, sep = " "))

  }


  if(type == "variance"){

    data_vis <- object$var_gam$data %>%
      dplyr::select(.data$Timestamp, {{x}}, .data$X_Ex2, {{z_numeric}}, {{z_factors}}) %>%
      tidyr::drop_na()
    # getting each component of the linear predictor from the fitted model
    fv_var <- mgcv::predict.gam(var_fit, type = "terms")

    if(is.null(z_factors)){
      colnames(fv_var) <- paste("s_", names_z_numeric, sep = "")
    }else{
      colnames(fv_var) <- c(paste("s_", names_z_factors, sep = ""),
                            paste("s_", names_z_numeric, sep = ""))
    }

    # residuals
    resid_var <- mgcv::residuals.gam(var_fit, type = "working")

    # combining linear predictor and residuals
    resid_var <- cbind(fv_var, resid_var) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(Timestamp = data_vis$Timestamp)

    ##-- visualising fitted variance model
    plot_var <- list()

    for (i in 1:length(names_z_numeric)) {
      ##-- preparing data for visualisation for numerical variables
      data_vis_z <- data_vis %>%
        dplyr::mutate_at(names_z_numeric[-i],
                  ~median(.,na.rm = T))
      if(!is.null(z_factors)){
        data_vis_z <- data_vis_z %>%
          dplyr::mutate_at(names_z_factors, ~names(which.max(table(.))))
      }

      # predicted values
      predict_fit_var <- mgcv::predict.gam(var_fit, type = "link",
                                           newdata = data_vis_z,
                                           se.fit = T) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(LI = as.numeric(.data$fit - .data$se.fit),
                      UI = as.numeric(.data$fit + .data$se.fit),
                      Timestamp = data_vis_z$Timestamp)

      s_name_zi <- paste("s_", names_z_numeric[i],
                         sep = "")
      data_vis_z <- list(data_vis_z, resid_var, predict_fit_var) %>%
        purrr::reduce(dplyr::left_join, by = "Timestamp") %>%
        dplyr::mutate(partial_resid = exp(resid_var) + !!rlang::sym(s_name_zi))

      name_zi <- names_z_numeric[i]
      plot_var[[i]] <- data_vis_z %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = !!rlang::sym(name_zi),
                                         y = .data$partial_resid),
                            color = "#666666", size = 0.5) +
        ggplot2::geom_line(color = "#0099FF",
                           ggplot2::aes(x = !!rlang::sym(name_zi),
                                        y = (.data$fit)),
                           size = 1) +
        ggplot2::geom_ribbon(ggplot2::aes(x = !!rlang::sym(name_zi),
                                          y = (.data$fit),
                                          ymin = .data$LI,
                                          ymax = .data$UI),
                             alpha = 0.5,
                             fill = "grey70") +
        ggplot2::ylab(paste("f(",name_zi,")", sep = ""))

    }

    gridExtra::grid.arrange(grobs = plot_var,
                            ncol = n_col,
                                           top = paste("conditional variance of",
                                                       names_x, sep = " "))

  }

}
