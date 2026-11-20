# contrastName

Returns the contrast that produced a result.

## Usage

``` r
contrastName(results)

# S4 method for class 'RegionSetDE.results'
contrastName(results)

# S4 method for class 'RegionSetDE.setResults'
contrastName(results)

# S4 method for class 'RegionSetDE.geneResults'
contrastName(results)

# S4 method for class 'RegionSetDE.resultsList'
contrastName(results)

# S4 method for class 'RegionSetDE.setResultsList'
contrastName(results)
```

## Arguments

- results:

  `RegionSetDE.results`, `RegionSetDE.setResults` or
  `RegionSetDE.geneResults` object, or one of the two list classes
  holding several contrasts.

## Value

A string, or one string per contrast for a list of results.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

contrastName(results)
#> [1] "condition: SHR vs BN"
```
