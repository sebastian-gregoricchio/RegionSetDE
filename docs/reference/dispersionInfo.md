# dispersionInfo

Returns the summary of the dispersion estimated, or supplied, when a
model was fitted: the common value, whether it was held fixed, whether
the design had no replicates, and where the null rows came from.

## Usage

``` r
dispersionInfo(fit)

# S4 method for class 'RegionSetDE.fit'
dispersionInfo(fit)
```

## Arguments

- fit:

  `RegionSetDE.fit` object returned by
  [`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md).

## Value

A named list. Which elements it holds depends on the engine and on how
the dispersion was obtained.

## See also

[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
[`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md),
[`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
names(dispersionInfo(fit))
#> [1] "common"        "trended"       "fixed"         "no.replicates"
#> [5] "prior.df"      "source"        "holdout.index" "holdout.type" 
dispersionInfo(fit)$common
#> [1] 0.1394794
```
