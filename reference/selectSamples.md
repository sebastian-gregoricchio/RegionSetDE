# selectSamples

Restricts a `RegionSetDE.counts` object to a subset of its samples. The
conditions are written as they would be in
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html) and
are evaluated on the `colData`, so any column of the sample metadata can
be used. The regions are left untouched, only the columns change.

## Usage

``` r
selectSamples(
  counts,
  ...,
  samples = NULL,
  dropNormalization = TRUE,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- ...:

  Conditions evaluated on the `colData`, in the syntax of
  [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html),
  e.g. `mark == "H3K27ac"`. Several conditions are combined with AND.

- samples:

  Character vector with the names of the samples to keep, or a numeric
  vector of column positions. Applied before the conditions. Default:
  `NULL`, all the samples.

- dropNormalization:

  Logical value to indicate whether the normalisation stored in the
  object must be discarded. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.counts` object with the selected columns.

## Details

Marks, and more generally experiments run on different antibodies or
different assays, should not share a model. The dispersion, the dynamic
range and the meaning of the scaling factors all differ between them, so
a fit that pools them borrows information across rows that have nothing
to say about each other. The selection therefore belongs upstream of
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
rather than at test time, and this is why `dropNormalization` defaults
to `TRUE`: factors estimated over a set of samples that no longer exists
describe a library composition that no longer exists either. The raw
counts are never modified, so re-normalising costs one call.

The background bins stored in the metadata, when present, are subset
along with the regions.

## See also

[`splitSamples`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitSamples.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

brownNorway <- selectSamples(counts, condition == "BN", verbose = FALSE)
SummarizedExperiment::colData(brownNorway)
#> DataFrame with 2 rows and 7 columns
#>                                                 sample               bam.file
#>                                            <character>            <character>
#> lv-H3K4me3-BN-female-bio1-tech1 lv-H3K4me3-BN-female.. /home/s.gregoricchio..
#> lv-H3K4me3-BN-male-bio2-tech1   lv-H3K4me3-BN-male-b.. /home/s.gregoricchio..
#>                                 condition         sex biologicalReplicate
#>                                  <factor> <character>         <character>
#> lv-H3K4me3-BN-female-bio1-tech1        BN      female                bio1
#> lv-H3K4me3-BN-male-bio2-tech1          BN        male                bio2
#>                                 paired.end library.size
#>                                  <logical>    <numeric>
#> lv-H3K4me3-BN-female-bio1-tech1      FALSE       386378
#> lv-H3K4me3-BN-male-bio2-tech1        FALSE       400384

firstReplicates <- selectSamples(counts, biologicalReplicate == "bio2",
                                 verbose = FALSE)
SummarizedExperiment::colData(firstReplicates)
#> DataFrame with 2 rows and 7 columns
#>                                                sample               bam.file
#>                                           <character>            <character>
#> lv-H3K4me3-BN-male-bio2-tech1  lv-H3K4me3-BN-male-b.. /home/s.gregoricchio..
#> lv-H3K4me3-SHR-male-bio2-tech1 lv-H3K4me3-SHR-male-.. /home/s.gregoricchio..
#>                                condition         sex biologicalReplicate
#>                                 <factor> <character>         <character>
#> lv-H3K4me3-BN-male-bio2-tech1        BN         male                bio2
#> lv-H3K4me3-SHR-male-bio2-tech1       SHR        male                bio2
#>                                paired.end library.size
#>                                 <logical>    <numeric>
#> lv-H3K4me3-BN-male-bio2-tech1       FALSE       400384
#> lv-H3K4me3-SHR-male-bio2-tech1      FALSE       337600
```
