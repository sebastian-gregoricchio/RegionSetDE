# peakOccupancyTable

Crosses the consensus a region came from with what happened to it in the
test: how many of the regions called in both groups changed, how many of
the ones called in a single group did, and in which direction. It works
on the results of regions built by
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
which carry the occupancy of every group.

## Usage

``` r
peakOccupancyTable(
  results,
  contrast = NULL,
  groups = NULL,
  by = "group",
  set = NULL,
  FDR = NULL,
  log2FC = NULL
)
```

## Arguments

- results:

  `RegionSetDE.results` object, or a `RegionSetDE.resultsList` together
  with `contrast`.

- contrast:

  String with the name of a contrast, or its position, when `results`
  holds several of them. Default: `NULL`.

- groups:

  Character vector with the consensus groups the regions are classified
  by. Default: `NULL`, the two groups the contrast compares when they
  are among them, every group otherwise.

- by:

  String indicating what the regions are grouped by, either `"group"`,
  the combination of consensus groups covering the region, or
  `"samples"`, the number of samples that carried a peak on it. Default:
  `"group"`.

- set:

  Character vector with the names of the region sets kept. Default:
  `NULL`, all of them.

- FDR:

  Numeric value with the adjusted p-value below which a region counts as
  changed. Default: `NULL`, the threshold the test was run with.

- log2FC:

  Numeric value with the log2 fold change a region must reach. Default:
  `NULL`, the threshold the test was run with.

## Value

A data.frame with one row per occupancy class: the class itself, the
number of regions in it, how many of them came out `down`, `null` and
`up`, and the percentage that changed in either direction.

## Details

The table says where the differences sit, which is the question a
peak-based experiment usually starts from: whether the condition
rearranged the shared peaks or whether it mostly gained and lost whole
sites. A region called in one group only and coming out unchanged is
worth as much attention as the reverse, since it means the peak caller
drew a line the counts do not support.

What the table is not is a test of occupancy. Whether a peak is called
in a sample depends on the depth of that library and on the threshold of
the caller as much as on the chromatin, so the number of group-specific
regions is not an effect size and comparing those numbers between groups
is not evidence of anything. The evidence is in the columns beside them,
which come from the counts.

## See also

[`plotPeakOccupancy`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakOccupancy.md),
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
[`plotPeakUpset`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakUpset.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (requireNamespace("consensusRegions", quietly = TRUE)) {
  # AR binding without ligand and after 24 hours of R1881
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  sampleSheet <- dplyr::filter(sampleSheet, condition %in% c("DMSO", "R1881_24h"))

  consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)

  counts <- countReads(consensus, sampleSheet = sampleSheet, verbose = FALSE)
  counts <- countBackground(counts, binSize = 10000, verbose = FALSE)
  counts <- normalizeCounts(counts, method = "background", verbose = FALSE)

  fit <- fitRegions(counts, design = ~ condition, verbose = FALSE)
  results <- testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)

  peakOccupancyTable(results)

  # By number of samples calling a peak, whatever the group
  peakOccupancyTable(results, by = "samples")
}
#> calcNormFactors has been renamed to normLibSizes
#>   occupancy n.regions down null up percent.changed
#> 1         2        25    0   15 10            40.0
#> 2         3        61    0    8 53            86.9
#> 3         4         5    0    1  4            80.0
#> 4         5         2    0    0  2           100.0
#> 5         6         9    0    1  8            88.9
```
