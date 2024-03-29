% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_dt.R
\name{estimate_dt}
\alias{estimate_dt}
\title{Estimating time delay between two sensors in a river system}
\usage{
estimate_dt(x)
}
\arguments{
\item{x}{Model object of class "conditional_ccf" returned from
\code{\link[conduits]{conditional_ccf}}}
}
\value{
A \code{\link[tibble]{tibble}} with estimated time lag "dt"
and corresponding maximum cross-correlation
}
\description{
This function estimates the time that takes water to flow from an upstream location to a downstream
location conditional on the observed water-quality variables from the upstream sensor. That time lag is
defined as the lag that gives maximum cross-correlation conditional on upstream water-quality variables.
}
\examples{

old_ts <- NEON_PRIN_5min_cleaned |>
  dplyr::select(
    Timestamp, site, turbidity, level, temperature
  ) |>
  tidyr::pivot_wider(
    names_from = site,
    values_from = turbidity:temperature
  )

fit_mean_y <- old_ts |>
  conditional_mean(turbidity_downstream ~
    s(level_upstream, k = 5) +
    s(temperature_upstream, k = 5))

fit_var_y <- old_ts |>
  conditional_var(
    turbidity_downstream ~
      s(level_upstream, k = 4) +
      s(temperature_upstream, k = 4),
    family = "Gamma",
    fit_mean = fit_mean_y
  )

fit_mean_x <- old_ts |>
  conditional_mean(turbidity_upstream ~
    s(level_upstream, k = 5) +
    s(temperature_upstream, k = 5))

fit_var_x <- old_ts |>
  conditional_var(
    turbidity_upstream ~
      s(level_upstream, k = 4) +
      s(temperature_upstream, k = 4),
    family = "Gamma",
    fit_mean = fit_mean_x
  )

fit_c_ccf <- old_ts |>
  tidyr::drop_na() |>
  conditional_ccf(
    I(turbidity_upstream * turbidity_downstream) ~
      splines::ns(level_upstream, df = 3) +
      splines::ns(temperature_upstream, df = 3),
    lag_max = 10,
    fit_mean_x = fit_mean_x, fit_var_x = fit_var_x,
    fit_mean_y = fit_mean_y, fit_var_y = fit_var_y,
    df_correlation = c(3, 3)
  )

new_data <- fit_c_ccf |> estimate_dt()
}
\author{
Puwasala Gamakumara & Priyanga Dilini Talagala
}
