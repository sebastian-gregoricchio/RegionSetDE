# Gene result accessors

Return how a `RegionSetDE.geneResults` object was built: the aggregation
of the regions of a gene and what weighted them. The table itself comes
from
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
the links from
[`regionLinks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md).

## Usage

``` r
geneMethod(results)

# S4 method for class 'RegionSetDE.geneResults'
geneMethod(results)

geneWeighting(results)

# S4 method for class 'RegionSetDE.geneResults'
geneWeighting(results)
```

## Arguments

- results:

  `RegionSetDE.geneResults` object returned by
  [`aggregateToGenes`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/aggregateToGenes.md).

## Value

`geneMethod` returns the aggregations that were run, `"simes"`,
`"potential"` or both. `geneWeighting` returns what weighted the regions
of a gene, one among `"link.weight"`, `"distance"` and `"none"`.

## See also

[`aggregateToGenes`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/aggregateToGenes.md),
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
[`regionLinks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# Links of your own, four regions per gene
regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fitCounts(fit)))[, 1:3]
geneNumber <- ceiling(nrow(regionTable) / 4)

annotatedFit <- annotateRegions(fit,
    linkTable = data.frame(regionTable,
        gene.id = rep(paste0("gene", seq_len(geneNumber)),
            each = 4
        )[seq_len(nrow(regionTable))]
    ),
    verbose = FALSE
)

results <- testRegions(annotatedFit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
geneResults <- aggregateToGenes(results, method = c("simes", "potential"), verbose = FALSE)

geneMethod(geneResults)
#> [1] "simes"     "potential"
geneWeighting(geneResults)
#> [1] "none"
```
