# RegionSetDE.resultsList class

S4 class holding the outcome of several contrasts run on the same fit.
Since all of them come from one model, the design, the offsets and the
dispersion are shared and the contrasts can be compared with each other
directly. Every function that takes a `RegionSetDE.results` object also
takes this one, together with a `contrast` argument naming which of them
to use.

## Usage

``` r
# S4 method for class 'RegionSetDE.resultsList'
show(object)
```

## Arguments

- object:

  `RegionSetDE.resultsList` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `results`:

  Named list of `RegionSetDE.results` objects.

- `contrasts`:

  Character vector with the names of the contrasts, in the order they
  were run.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# A named list of contrasts gives one result per contrast
resultsList <- testRegions(fit,
    contrast = list(shr = c("condition", "SHR", "BN"), bn = c("condition", "BN", "SHR")),
    verbose = FALSE
)
resultsList
#> An object of class 'RegionSetDE.resultsList'
#>   contrasts       : 2 
#> 
#>  name             contrast n.regions up down
#>   shr condition: SHR vs BN      1895  3    9
#>    bn condition: BN vs SHR      1895  9    3
resultsList[["shr"]]
#> An object of class 'RegionSetDE.results'
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   regions         : 1895 
#>   counts carried  : 4 samples
#>   thresholds      : FDR < 0.05 | |log2FC| > 0 
#>   changing regions:
#>     promoterNonCpG: 1 up, 2 down
#>     intergenic: 2 up, 3 down
#>     geneBody: 0 up, 3 down
#>     promoterCpG: 0 up, 1 down
```
