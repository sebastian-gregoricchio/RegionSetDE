# RegionSetDE.setResultsList class

S4 class holding the set level outcome of several contrasts run on the
same fit. Every function that takes a `RegionSetDE.setResults` object
also takes this one, together with a `contrast` argument naming which of
them to use.

## Usage

``` r
# S4 method for class 'RegionSetDE.setResultsList'
show(object)
```

## Arguments

- object:

  `RegionSetDE.setResultsList` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `results`:

  Named list of `RegionSetDE.setResults` objects.

- `contrasts`:

  Character vector with the names of the contrasts, in the order they
  were run.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# A named list of contrasts gives one set level result per contrast
setResultsList <- testRegionSets(fit,
    contrast = list(shr = c("condition", "SHR", "BN"), bn = c("condition", "BN", "SHR")),
    verbose = FALSE
)
setResultsList
#> An object of class 'RegionSetDE.setResultsList'
#>   contrasts       : 2 
#> 
#>  name             contrast n.sets n.significant
#>   shr condition: SHR vs BN      4             0
#>    bn condition: BN vs SHR      4             0
names(setResultsList)
#> [1] "shr" "bn" 
```
