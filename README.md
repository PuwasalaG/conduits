
<!-- README.md is generated from README.Rmd. Please edit that file -->

# conduits (CONDitional UI for Time Series normalisation)

<img src="man/figures/logo.png" align="right" height="190"/>
<!-- badges: start --> <!-- badges: end -->

Package `conduits` provides an user interface for conditionally
normalising a time series. This also facilitates functions to produce
conditional cross-correlations between two normalised time series at
different lags while providing some graphical tools for visualisation.

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

`conduits` contains a set of water-quality variables measured by in-situ
sensors from the Pringle Creek located in Wise County, Texas. This is
one of the aquatic NEON field sites hosted by the US Forest Service.

This data contains water-quality variables such as, turbidity, specific
conductance, dissolved oxygen, pH and fDOM along with surface elevation
and surface temperature from two sites located about 200m apart. Data
are available from
![2019-07-01](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2019-07-01 "2019-07-01")
to
![2019-12-31](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2019-12-31 "2019-12-31")
at every 5 mintues.

In this example we choose turbidity from upstream and downstream sites
to calculate the cross-correlation while conditioning on level,
temperature and conductance from the upstream location.

Let us first prepare data as follows

``` r
library(conduits)
library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.6     ✓ dplyr   1.0.8
#> ✓ tidyr   1.2.0     ✓ stringr 1.4.0
#> ✓ readr   2.1.2     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

# Data

`conduits` contains `NEON_PRIN_5min_cleaned`, a set of water-quality
variables measured by in-situ sensors from Pringle Creek located in Wise
County, Texas. This is one of the aquatic NEON field sites hosted by the
US Forest Service.

This data contains water-quality variables such as turbidity, specific
conductance, dissolved oxygen, pH and fDOM along with surface elevation
and surface temperature from two sites located about 200m apart. Data
are available from
![2019-07-01](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2019-07-01 "2019-07-01")
to
![2019-12-31](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2019-12-31 "2019-12-31")
at every 5 minutes.

In this example, we choose turbidity from upstream and downstream sites
to calculate the cross-correlation while conditioning on the water
level, temperature and conductance from the upstream location.

Let us first prepare data as follows

``` r
data <- NEON_PRIN_5min_cleaned %>%
  dplyr::filter(site == "upstream") %>%
  dplyr::select(Timestamp, turbidity, level,
                conductance, temperature) %>%
  tsibble::as_tsibble(index = Timestamp)
head(data)
#> # A tsibble: 6 x 5 [5m] <UTC>
#>   Timestamp           turbidity level conductance temperature
#>   <dttm>                  <dbl> <dbl>       <dbl>       <dbl>
#> 1 2019-07-01 00:00:00     1.55   251.        821.        29.6
#> 2 2019-07-01 00:05:00     1.1    251.        821.        29.6
#> 3 2019-07-01 00:10:00     0.855  251.        821.        29.5
#> 4 2019-07-01 00:15:00     0.87   251.        821.        29.5
#> 5 2019-07-01 00:20:00     1.01   251.        821.        29.4
#> 6 2019-07-01 00:25:00     0.828  251.        821.        29.4
```

### Conditional normalisation

The following code shows how to use the `conditional_moments` function
to normalise turbidity from upstream sites.

``` r
# Estimating conditional mean of the turbidity from the upstream site 

fit_mean <- data %>%
  conditional_mean(
    turbidity ~ s(level, k = 8) + s(conductance, k = 8) + s(temperature, k = 8))

summary(fit_mean)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> turbidity ~ s(level, k = 8) + s(conductance, k = 8) + s(temperature, 
#>     k = 8)
#> 
#> Parametric coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  4.25221    0.02438   174.4   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                  edf Ref.df       F p-value    
#> s(level)       6.945  6.998 3059.76  <2e-16 ***
#> s(conductance) 6.974  7.000 1596.37  <2e-16 ***
#> s(temperature) 6.979  7.000   65.78  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.687   Deviance explained = 68.7%
#> GCV = 29.604  Scale est. = 29.591    n = 49796
class(fit_mean)
#> [1] "conditional_moment" "gam"                "glm"               
#> [4] "lm"

# Estimating conditional variance of the turbidity from the upstream site 

fit_var <- data %>%
  conditional_var(
    turbidity ~ s(level, k = 7) + s(conductance, k = 7) + s(temperature, k = 7),
    family = "Gamma",
    fit_mean
  )

class(fit_var)
#> [1] "conditional_moment" "gam"                "glm"               
#> [4] "lm"
summary(fit_var)
#> 
#> Family: Gamma 
#> Link function: log 
#> 
#> Formula:
#> Y_Ey2 ~ s(level, k = 7) + s(conductance, k = 7) + s(temperature, 
#>     k = 7)
#> 
#> Parametric coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  1.86626    0.06525    28.6   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                  edf Ref.df     F p-value    
#> s(level)       5.977  6.000 59.28  <2e-16 ***
#> s(conductance) 5.855  5.990 18.93  <2e-16 ***
#> s(temperature) 5.960  5.999 33.12  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  -0.0112   Deviance explained =   38%
#> GCV = 4.9694  Scale est. = 212.02    n = 49796

# Normalize the series using conditional moments
new_ts <- data %>%
  dplyr::mutate(
    ystar = normalize(., turbidity, fit_mean, fit_var))
```
