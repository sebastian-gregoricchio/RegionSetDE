# RegionSetDE.geneResults class

S4 class holding the per-gene table built by
[`aggregateToGenes`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/aggregateToGenes.md)
from a per-region result. One row per gene, with the regions that
reached it counted and their directions reported rather than resolved,
and the filters applied upstream still readable in the inherited
provenance slots.

## Usage

``` r
# S4 method for class 'RegionSetDE.geneResults'
show(object)
```

## Arguments

- object:

  `RegionSetDE.geneResults` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `results`:

  `DataFrame` with one row per gene. Read it with
  [`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md).

- `contrast`:

  String with the contrast the regions were tested on.

- `method`:

  String with the aggregation used, `"simes"` or `"potential"`.

- `links`:

  `DataFrame` with the region to gene links the aggregation was built
  on, one row per link.

- `weighting`:

  String saying what weighted the regions of a gene, one among
  `"link.weight"`, `"distance"` and `"none"`.

- `engine`:

  String with the engine the regions were fitted with.

- `thresholds`:

  List with the cut-offs the directions were read at.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# Links of your own, four regions per gene
regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fitCounts(fit)))[, 1:3]
geneNumber <- ceiling(nrow(regionTable) / 4)
geneIds <- rep(paste0("gene", seq_len(geneNumber)), each = 4)[seq_len(nrow(regionTable))]

annotatedFit <- annotateRegions(fit, linkTable = data.frame(regionTable, gene.id = geneIds), verbose = FALSE)
results <- testRegions(annotatedFit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

geneResults <- aggregateToGenes(results, method = "simes", verbose = FALSE)
geneResults
#> An object of class 'RegionSetDE.geneResults'
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   aggregation     : simes 
#>   weighting       : none 
#>   genes           : 474 
#>   links           : 1895 
#>   direction       : up 3, down 9, mixed 0, unchanged 462 
#> 
#>  gene.id gene.symbol n.regions n.up n.down direction  p.value      FDR
#>   gene54        <NA>         4    0      1      down 4.59e-08 2.18e-05
#>  gene166        <NA>         4    0      1      down 1.10e-05 2.60e-03
#>    gene4        <NA>         4    0      1      down 5.48e-05 6.58e-03
#>  gene309        <NA>         4    0      1      down 5.55e-05 6.58e-03
#>  gene295        <NA>         4    0      1      down 1.34e-04 1.27e-02
#>  gene162        <NA>         4    1      0        up 1.68e-04 1.33e-02
```
