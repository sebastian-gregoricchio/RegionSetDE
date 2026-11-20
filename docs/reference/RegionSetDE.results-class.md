# RegionSetDE.results class

S4 class holding the outcome of a per-region contrast. The table in the
`results` slot has one row per region and follows the coordinates in the
`regions` slot, so the two can be bound together without any matching
step. When the counts were tiled, the tile level statistics are kept in
the `tiles` slot next to the combined ones, since a region that changes
only over part of its length is easier to read from the tiles than from
the summary.

## Usage

``` r
# S4 method for class 'RegionSetDE.results'
show(object)
```

## Arguments

- object:

  `RegionSetDE.results` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `results`:

  `DataFrame` with one row per region. Read it with
  [`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md).

- `tiles`:

  `DataFrame` with one row per tile, empty when the counts were not
  tiled. Read it with
  [`tileTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/tileTable.md).

- `regions`:

  `GRanges` with the coordinates of the rows of `results`.

- `contrast`:

  String describing the contrast that was tested.

- `contrast.vector`:

  Numeric vector with the coefficients of the contrast over the columns
  of the design.

- `contrast.groups`:

  List with the `column` of the `colData` the contrast separates and the
  two `groups` it compares, the first one being the level the fold
  change is positive for, and `n.samples`, the number of samples in each
  of the two. Empty when the contrast is not a difference between two
  levels of one variable.
  [`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
  reads it to draw the model brackets, and
  [`contrastInfo`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastInfo.md)
  reports it.

- `engine`:

  String with the engine that produced the statistics.

- `counting.level`:

  String indicating whether the model was fitted on regions or on tiles.

- `combination`:

  List with the method used to combine the tiles and whether it was
  applied.

- `thresholds`:

  List with the FDR and log2 fold change cut-offs used to fill the
  `diff.status` column.

- `counts`:

  `RegionSetDE.counts` object the contrast was computed on, so that the
  values behind a result can be drawn without carrying a second object
  around. Empty when `testRegions` was called with
  `carryCounts = FALSE`.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
results
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
head(resultsTable(results), 3)
#>       region.set    region.id tile.id seqnames start   end width     log2FC
#> 1 promoterNonCpG region_00012      NA    chr12 26988 27987  1000 -0.4744875
#> 2 promoterNonCpG region_00017      NA    chr12 39449 40448  1000 -1.6911510
#> 3 promoterNonCpG region_00019      NA    chr12 44116 45115  1000  0.1328365
#>   average.signal average.signal.BN average.signal.SHR       stat
#> 1       3.098695          3.154609           3.071096 0.16629041
#> 2       3.141973          3.512708           2.774385 2.33758830
#> 3       3.273341          3.149796           3.390540 0.01002874
#>   stat.distribution df1      df2   p.value       FDR diff.status     regionId
#> 1                 f   1 18.94543 0.6880013 0.9200864        null region_00012
#> 2                 f   1 18.73556 0.1429932 0.7997195        null region_00017
#> 3                 f   1 18.45274 0.9213104 0.9807286        null region_00019
```
