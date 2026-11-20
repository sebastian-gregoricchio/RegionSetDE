# Accessors of RegionSetDE.resultsList

Extract one contrast, all of them, or the names of the contrasts, from a
`RegionSetDE.resultsList` or a `RegionSetDE.setResultsList` object.

## Usage

``` r
# S4 method for class 'RegionSetDE.setResultsList,ANY,ANY'
x[[i, j, ...]]

# S4 method for class 'RegionSetDE.setResultsList'
x$name

# S4 method for class 'RegionSetDE.setResultsList'
names(x)

# S4 method for class 'RegionSetDE.setResultsList'
length(x)

# S4 method for class 'RegionSetDE.resultsList,ANY,ANY'
x[[i, j, ...]]

# S4 method for class 'RegionSetDE.resultsList'
x$name

# S4 method for class 'RegionSetDE.resultsList'
names(x)

# S4 method for class 'RegionSetDE.resultsList'
length(x)

# S4 method for class 'RegionSetDE.resultsList'
as.list(x, ...)

# S4 method for class 'RegionSetDE.setResultsList'
as.list(x, ...)
```

## Arguments

- x:

  `RegionSetDE.resultsList` or `RegionSetDE.setResultsList` object.

- i:

  String with the name of a contrast, or its position.

- j:

  Not used, present because the `[[` generic carries it.

- ...:

  Not used.

- name:

  String with the name of a contrast.

## Value

A `RegionSetDE.results` object, or a `RegionSetDE.setResults` object for
the list of set level results. `names` returns a character vector,
`length` the number of contrasts, and `as.list` a named list with every
result, one per contrast.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

resultsList <- testRegions(fit,
    contrast = list(shr = c("condition", "SHR", "BN"), bn = c("condition", "BN", "SHR")),
    verbose = FALSE
)
names(resultsList)
#> [1] "shr" "bn" 
length(resultsList)
#> [1] 2
resultsList$bn
#> An object of class 'RegionSetDE.results'
#>   contrast        : condition: BN vs SHR 
#>   engine          : edgeR 
#>   regions         : 1895 
#>   counts carried  : 4 samples
#>   thresholds      : FDR < 0.05 | |log2FC| > 0 
#>   changing regions:
#>     promoterNonCpG: 2 up, 1 down
#>     intergenic: 3 up, 2 down
#>     geneBody: 3 up, 0 down
#>     promoterCpG: 1 up, 0 down
```
