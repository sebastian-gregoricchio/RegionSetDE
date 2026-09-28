# updateThresholds

Labels the regions of a result again with other cut-offs, without
running the test a second time. The `diff.status` column is filled anew
from the `FDR` and `log2FC` given, and the thresholds stored in the
object are replaced, so that the tables, the plots and the export all
follow the new ones.

## Usage

``` r
updateThresholds(
  results,
  FDR = NULL,
  log2FC = NULL,
  contrast = NULL,
  verbose = TRUE
)
```

## Arguments

- results:

  `RegionSetDE.results` or `RegionSetDE.resultsList` object returned by
  [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
  or `RegionSetDE.setResults` or `RegionSetDE.setResultsList` object
  returned by
  [`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  or
  [`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md).

- FDR:

  Numeric value with the new adjusted p-value cut-off. Default: `NULL`,
  the one stored in the object.

- log2FC:

  Numeric value with the new absolute log2 fold change cut-off. Not
  available for the results of the region sets, which carry no label per
  region. Default: `NULL`, the one stored in the object.

- contrast:

  Character vector with the names or the positions of the contrasts to
  update, when `results` holds several of them. Default: `NULL`, all of
  them.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

An object of the same class as `results`, with the new labels in
`diff.status` and the new cut-offs in the `thresholds` slot and in the
parameters recorded for the test. The previous cut-offs are kept in
`parameters$updateThresholds`.

## Details

The cut-offs only decide the labels. The p-values, the adjusted p-values
and the fold changes come from the model and do not move, which is why
nothing has to be fitted or tested again. What does need a new test is
`lfcThreshold` of
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
which moves the fold change inside the test itself and changes the
p-values; it stays as it was and is reported next to the new cut-offs.

Everything reading the labels afterwards sees the new ones:
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md),
[`plotVolcano`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md)
and
[`plotResultsMA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotResultsMA.md),
which draw the dashed lines at the stored cut-offs,
[`peakOccupancyTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/peakOccupancyTable.md),
[`plotProfile`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotProfile.md)
and
[`exportResults`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/exportResults.md),
which writes the cut-offs down with the rest of the parameters. On tiled
results the region labels are updated, while `n.tiles.up` and
`n.tiles.down` keep counting the tiles changing within the region at the
fixed rate csaw uses, since they come from the combination of the tiles
and not from the labels.

For the results of the region sets only `FDR` applies: it decides which
sets
[`plotSetEffect`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md)
and the printed summary call significant.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

table(resultsTable(results)$diff.status)
#> 
#> down null   up 
#>    9 1883    3 

# Stricter cut-offs, with no new test
strictResults <- updateThresholds(results, FDR = 0.01, log2FC = 1)
#> Relabelled with FDR < 0.01 and |log2FC| > 1: 0 up and 4 down out of 1895 regions.
table(resultsTable(strictResults)$diff.status)
#> 
#> down null   up 
#>    4 1891    0 

# The cut-offs travel with the object, the plots draw them
plotVolcano(strictResults)

```
