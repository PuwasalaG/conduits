#' Visualising the models fitted for conditional moments
#'
#' A function for visualising the fitted generalised additive models for
#' conditional means and variances from \code{conditional_moments} function
#'
#' @param object an object returned from \code{conditional_moments} function
#' @param type the type of moment want to visualise. This can take one of
#' "mean" or "variance"
#' @param title if TRUE returning the plot with a title
#' @param ... ignored
#'
#' @return returns plots visualising fitted gam model for conditional
#' mean or variance vs each predictor.
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_ribbon ylab geom_rug
#'
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[visreg]{visreg}}
#'
#' @author Puwasala Gamakumara
#'
#' @export
autoplot.conditional_moments <- function(object,
                                         type = c("mean", "variance"),
                                         title = TRUE,
                                         ...){
  type <- match.arg(type)
  x <- object$x
  z_numeric <- object$z_numeric
  z_factors <- object$z_factors

  mean_fit <- object$mean_gam$fit
  var_fit <- object$var_gam$fit

  names_x <- names(tidyselect::eval_select(dplyr::enquo(x), object$data_conditional_moments))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric),
                                                   object$data_conditional_moments))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors),
                                                   object$data_conditional_moments))

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
    summary_mean <- object$mean_gam$fit %>% summary
    mean_edf <- summary_mean$edf %>% round(digits = 3)

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
      edf_zi <- mean_edf[i]
      range_zi <- data_vis_z %>% pull(!!rlang::sym(name_zi))

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
        ggplot2::ylab(paste("f(",name_zi, ",", edf_zi, ")", sep = "")) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()(range_zi))

      data[[i]] <- data_vis_z %>%
        dplyr::select(.data$Timestamp)

    }

    if(title==TRUE){
      gridExtra::grid.arrange(grobs = plot_means,
                                           ncol = n_col,
                                           top = paste("conditional mean of",
                                                       names_x, sep = " "))
    }else{
      gridExtra::grid.arrange(grobs = plot_means,
                              ncol = n_col)
    }
  }


  if(type == "variance"){

    data_vis <- object$var_gam$data %>%
      dplyr::select(.data$Timestamp, {{x}}, .data$X_Ex2, {{z_numeric}}, {{z_factors}}) %>%
      tidyr::drop_na()
    summary_var <- object$var_gam$fit %>% summary
    var_edf <- summary_var$edf %>% round(digits = 3)
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
      edf_zi <- var_edf[i]
      range_zi <- data_vis_z %>% pull(!!rlang::sym(name_zi))

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
        ggplot2::ylab(paste("f(",name_zi,",",edf_zi, ")", sep = "")) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()(range_zi))

    }

    if(title==TRUE){
      gridExtra::grid.arrange(grobs = plot_var,
                              ncol = n_col,
                              top = paste("conditional variance of",
                                          names_x, sep = " "))
    }else{
      gridExtra::grid.arrange(grobs = plot_var,
                              ncol = n_col)
    }
  }

}


#' Visualising the models fitted for conditional cross-correlations
#'
#' A function for visualising the fitted generalised additive models for
#' conditional means, variances and cross-correlations from
#' \code{conditional_ccf} function
#'
#' @param object an object returned from \code{conditional_ccf} function
#' @param type the type of moment or cross-correlation to visualise.
#' This can take one of "mean", "variance" or "cross-correlation"
#' @param ... further arguments passed to or from other methods
#'
#' @return returns plots visualising fitted gam models for conditional
#' mean or variance vs each predictor.
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_ribbon ylab geom_rug
#'
#' @seealso \code{\link[visreg]{visreg}}
#'
#' @author Puwasala Gamakumara
#'
#' @export
autoplot.conditional_ccf <- function(object,
                                     type = c("mean", "variance",
                                              "cross-correlation"),
                                     ...){

  type <- match.arg(type)
  # x <- object$other$x
  # y <- object$other$y

  z_numeric <- object$other$z_numeric
  z_factors <- object$other$z_factors

  k <- object$other$k
  cond_moment_x <- object$data_visualise$conditional_moments$x
  cond_moment_y <- object$data_visualise$conditional_moments$y

  # names_x <- names(eval_select(enquo(x), object$data_ccf))
  names_z_numeric <- names(tidyselect::eval_select(dplyr::enquo(z_numeric), object$data_ccf))
  names_z_factors <- names(tidyselect::eval_select(dplyr::enquo(z_factors), object$data_ccf))

  p <- length(names_z_numeric) + length(names_z_factors)

  if(p <= 2){
    n_col <- 2
  }

  if(p > 2){
    n_col <- 3
  }

  #Conditional mean of x

  oask <- grDevices::devAskNewPage(TRUE)
  on.exit(grDevices::devAskNewPage(oask))


  # visualising conditional means
  if(type == "mean"){
    autoplot.conditional_moments(object = cond_moment_x,
                                 type = "mean")
    autoplot.conditional_moments(object = cond_moment_y,
                                 type = "mean")
  }

  # conditional variance of x

  # visualising conditional variance
  if(type == "variance"){
    autoplot.conditional_moments(object = cond_moment_x,
                                 type = "variance")
    autoplot.conditional_moments(object = cond_moment_y,
                                 type = "variance")
  }

  # visualising conditional cross-correlations
  if(type == "cross-correlation"){

    for (j in seq_along(k)) {

      ccf_gam_fit_k <- object$data_visualise$conditional_ccf$ccf_gam_fit[[j]]

      data_vis <- ccf_gam_fit_k$data %>%
        tidyr::drop_na()

      # getting each component of the linear predictor from the fitted model
      fv_ccf <- stats::predict.glm(ccf_gam_fit_k, type = "terms")

      if(is.null(z_factors)){
        colnames(fv_ccf) <- paste("s_", names_z_numeric, sep = "")
      }else{
        colnames(fv_ccf) <- c(paste("s_", names_z_factors, sep = ""),
                              paste("s_", names_z_numeric, sep = ""))
      }

      # residuals
      resid_ccf <- stats::residuals.glm(ccf_gam_fit_k, type = "working")

      # combining linear predictor and residuals
      resid_ccf <- cbind(fv_ccf, resid_ccf) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Timestamp = data_vis$Timestamp)

      plot_link <- list()


      plot_response <- list()

      ## link
      linkfun <- function(mu) {log((1+mu)/(1-mu))}
      ## inverse link
      linkinv <- function(eta) {(exp(eta) - 1)/(exp(eta) + 1)}

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
        predict_fit_ccf_link <- stats::predict.glm(ccf_gam_fit_k, type = "link",
                                            newdata = data_vis_z,
                                            se.fit = T) %>%
          tibble::as_tibble() %>%
          dplyr::rename("fit_link" = .data$fit,
                 "se.fit_link" = .data$se.fit) %>%
          dplyr::mutate(LI_link = as.numeric((.data$fit_link - .data$se.fit_link)),
                 UI_link = as.numeric((.data$fit_link + .data$se.fit_link)),
                 Timestamp = data_vis_z$Timestamp) %>%
          dplyr::select(-.data$residual.scale)

        predict_fit_ccf_res <- stats::predict.glm(ccf_gam_fit_k,
                                           type = "response",
                                           newdata = data_vis_z,
                                           se.fit = T) %>%
          tidyr::as_tibble() %>%
          dplyr::rename("fit_res" = .data$fit,
                 "se.fit_res" = .data$se.fit) %>%
          dplyr::mutate(LI_res = as.numeric((.data$fit_res - .data$se.fit_res)),
                 UI_res = as.numeric((.data$fit_res + .data$se.fit_res)),
                 Timestamp = data_vis_z$Timestamp) %>%
          dplyr::select(-.data$residual.scale)

        s_name_zi <- paste("s_", names_z_numeric[i],
                           sep = "")
        data_vis_z <- list(data_vis_z, resid_ccf, predict_fit_ccf_link,
                           predict_fit_ccf_res) %>%
          purrr::reduce(dplyr::left_join, by = "Timestamp")

        name_zi <- names_z_numeric[i]
        range_zi <- data_vis_z %>% pull(!!rlang::sym(name_zi))

        plot_link[[i]] <- data_vis_z %>%
          ggplot() +
          geom_point(aes(x = !!rlang::sym(name_zi), y = .data$XY),
                     color = "#666666", size = 0.5) +
          geom_line(color = "#0099FF", aes(x = !!rlang::sym(name_zi),
                                           y = .data$fit_link),
                    size = 1) +
          geom_ribbon(aes(x = !!rlang::sym(name_zi), y = .data$fit_link,
                          ymin = .data$LI_link, ymax = .data$UI_link),
                      alpha = 0.5,
                      fill = "grey70") +
          scale_x_continuous(breaks = scales::pretty_breaks()(range_zi)) +
          ylab(paste("f(",name_zi,")", sep = ""))

        plot_response[[i]] <- ggplot() +
          geom_line(data = data_vis_z,
                    color = "#0099FF", aes(x = !!rlang::sym(name_zi),
                                           y = .data$fit_res),
                    size = 1) +
          geom_ribbon(data = data_vis_z,
                      aes(x = !!rlang::sym(name_zi), y = .data$fit_res,
                          ymin = .data$LI_res, ymax = .data$UI_res),
                      alpha = 0.5,
                      fill = "grey70") +
          geom_rug(data = data_vis_z %>%
                     dplyr::filter(.data$resid_ccf >= 0),
                   aes(x = !!rlang::sym(name_zi), y = linkinv(.data$resid_ccf)),
                   sides = "t", color = "#666666") +
          geom_rug(data = data_vis_z %>%
                     dplyr::filter(.data$resid_ccf < 0),
                   aes(x = !!rlang::sym(name_zi), y = linkinv(.data$resid_ccf)),
                   sides = "b", color = "#666666") +
          scale_x_continuous(breaks = scales::pretty_breaks()(range_zi)) +
          ylab(bquote(c[.(j)]))

      }

      gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plot_link, ncol = n_col,
                               top = "predictor scale"),
                              gridExtra::arrangeGrob(grobs = plot_response, ncol = n_col,
                               top = "response scale"),
                              top = paste("k=", j, sep = ""))
    }


  }

}
