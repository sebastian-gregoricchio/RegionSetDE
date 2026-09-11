# splitSamples

Splits a `RegionSetDE.counts` object into a list of objects, one per
level of a column of the `colData`. Convenient when several marks or
assays have been counted together over the same regions and each of them
needs its own normalisation and its own fit.

## Usage

``` r
splitSamples(
  counts,
  by,
  dropNormalization = TRUE,
  minSamples = 1,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- by:

  String, or character vector, with the names of the `colData` columns
  defining the groups. When more than one is given the levels are
  combined.

- dropNormalization:

  Logical value to indicate whether the normalisation stored in the
  object must be discarded in each piece. Default: `TRUE`.

- minSamples:

  Numeric value with the minimum number of samples a group must contain
  to be returned. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A named list of `RegionSetDE.counts` objects.

## See also

[`selectSamples`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/selectSamples.md),
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

byStrain <- splitSamples(counts, by = "condition", verbose = FALSE)
names(byStrain)
#> [1] "BN"  "SHR"
byStrain[[1]]
#> class: RegionSetDE.counts 
#> dim: 3224 2 
#> metadata(2): signal.type background
#> assays(1): counts
#> rownames(3224): promoterNonCpG|region_00002 promoterNonCpG|region_00003
#>   ... promoterCpG|region_03797 promoterCpG|region_03798
#> rowData names(3): region.set region.id tile.id
#> colnames(2): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1
#> colData names(7): sample bam.file ... paired.end library.size
```
