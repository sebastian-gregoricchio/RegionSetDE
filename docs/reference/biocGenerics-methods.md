# Methods for the generics of BiocGenerics

The objects of the package answer the generics other Bioconductor
packages define for the same kind of data, so that code written for
`DESeq2` or `edgeR` objects reads them too. `counts` returns the raw
count matrix, `design` the design matrix of a fit, `plotMA` the MA plot
and `plotPCA` the principal component plot. They call the functions of
the package of the same purpose, which are documented in full on their
own pages.

## Usage

``` r
# S4 method for class 'RegionSetDE.counts'
counts(object, ...)

# S4 method for class 'RegionSetDE.fit'
counts(object, ...)

# S4 method for class 'RegionSetDE.results'
counts(object, ...)

# S4 method for class 'RegionSetDE.setResults'
counts(object, ...)

# S4 method for class 'RegionSetDE.resultsList'
counts(object, ...)

# S4 method for class 'RegionSetDE.setResultsList'
counts(object, ...)

# S4 method for class 'RegionSetDE.fit'
design(object, ...)

# S4 method for class 'RegionSetDE.results'
plotMA(object, ...)

# S4 method for class 'RegionSetDE.resultsList'
plotMA(object, ...)

# S4 method for class 'RegionSetDE.counts'
plotMA(object, ...)

# S4 method for class 'RegionSetDE.counts'
plotPCA(object, ...)

# S4 method for class 'RegionSetDE.fit'
plotPCA(object, ...)

# S4 method for class 'RegionSetDE.results'
plotPCA(object, ...)

# S4 method for class 'RegionSetDE.resultsList'
plotPCA(object, ...)

# S4 method for class 'RegionSetDE.setResults'
plotPCA(object, ...)
```

## Arguments

- object:

  For `counts`, a `RegionSetDE.counts`, `RegionSetDE.fit`,
  `RegionSetDE.results` or `RegionSetDE.setResults` object, or one of
  the two list classes holding several contrasts. For `design`, a
  `RegionSetDE.fit` object. For `plotMA`, a `RegionSetDE.results` or
  `RegionSetDE.resultsList` object, drawn by
  [`plotResultsMA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotResultsMA.md),
  or a normalised `RegionSetDE.counts` object, drawn by
  [`plotSetMA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md).
  For `plotPCA`, any object accepted by
  [`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md).

- ...:

  For `plotMA` and `plotPCA`, the arguments of
  [`plotResultsMA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotResultsMA.md),
  [`plotSetMA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md)
  or
  [`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md),
  such as `set` or `contrast`. Not used by `counts` and `design`.

## Value

`counts` returns a numeric matrix with one row per region, or per tile,
and one column per sample, holding the raw counts. For a result they are
the counts it carries, which needs the test to have been run with
`carryCounts = TRUE`;
[`countTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countTable.md)
returns the normalised values or a table with the coordinates. `design`
returns the matrix of the fixed effects, one row per sample, as
[`designMatrix`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitAccessors.md)
does. `plotMA` and `plotPCA` return a `ggplot` object.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

dim(counts(fit))
#> [1] 1895    4
head(design(fit))
#>                                 (Intercept) conditionSHR
#> lv-H3K4me3-BN-female-bio1-tech1           1            0
#> lv-H3K4me3-BN-male-bio2-tech1             1            0
#> lv-H3K4me3-SHR-male-bio2-tech1            1            1
#> lv-H3K4me3-SHR-male-bio3-tech1            1            1

results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
plotMA(results)

plotPCA(fit)

```
