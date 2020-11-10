
<!-- README.md is generated from README.Rmd. Please edit that file -->

# conduits (CONDitional UI for Time Series normalization)

<!-- badges: start -->

<!-- badges: end -->

Package `conduits` provides an user interface for conditionally
normalising a time series. This also facilitates functions to produce
conditional cross-correlations between two normalised time series at
different lags while providing some graphycal tools for visualisation.

`conduits` can also be used to estimate the time delay between two
sensor locations in river systems.

## Installation

You can install `conduits` from github with:

``` r
# install.packages("devtools")
devtools::install_github("PuwasalaG/conduits")
```

## Example

This is a basic example which shows you how to use functions in
`conduits`:

Still working on the
examples

<!-- ```{r example} -->

<!-- library(conduits) -->

<!-- library(lubridate) -->

<!-- head(NEON_PRIN_5min) -->

<!-- ## Normalising upstream turbidity conditional on upstream conductance, temperature and level -->

<!-- NEON_PRIN_5min <- NEON_PRIN_5min %>% -->

<!--   mutate(Timestamp = ymd_hms(roundedTime)) %>%  -->

<!--   filter(Timestamp >= ymd("2019-10-01")  -->

<!--          & Timestamp < ymd("2020-01-01")) %>%  -->

<!--   select(Timestamp, site, turbidity, surfacewaterElevMean, -->

<!--          surfWaterTempMean, ) -->

<!-- ``` -->