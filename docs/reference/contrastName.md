# contrastName

Returns the contrast that produced a `RegionSetDE.results` object.

## Usage

``` r
contrastName(results)

# S4 method for class 'RegionSetDE.results'
contrastName(results)

# S4 method for class 'RegionSetDE.setResults'
contrastName(results)

# S4 method for class 'RegionSetDE.resultsList'
contrastName(results)

# S4 method for class 'RegionSetDE.setResultsList'
contrastName(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A string.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

contrastName(results)
#> [1] "condition: SHR vs BN"
```
