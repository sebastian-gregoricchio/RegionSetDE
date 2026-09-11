# applyBlacklist

Removes from every region set the regions overlapping a blacklist, such
as the ENCODE blacklisted regions. The blacklist can be provided as a
BED-like file, a `GRanges` or a data.frame.

## Usage

``` r
applyBlacklist(
  regionSet,
  blacklist,
  overlapType = "any",
  minOverlapBp = 1,
  minOverlapFraction = 0,
  trimRegions = FALSE,
  ignoreStrand = TRUE,
  emptySets = "stop",
  verbose = TRUE
)
```

## Arguments

- regionSet:

  A `RegionSetDE` object, a `GRangesList`, a named list of `GRanges` or
  a single `GRanges`.

- blacklist:

  String indicating the path to a BED-like file, a `GRanges` or a
  data.frame with the regions to exclude.

- overlapType:

  String indicating the type of overlap required to blacklist a region,
  one among `"any"`, `"within"`, `"start"`, `"end"` or `"equal"`.
  Default: `"any"`.

- minOverlapBp:

  Numeric value indicating the minimum number of bases that must overlap
  the blacklist for a region to be removed. Default: `1`.

- minOverlapFraction:

  Numeric value between 0 and 1 indicating the minimum fraction of a
  region that must overlap the blacklist for it to be removed. Default:
  `0`, any overlap is sufficient.

- trimRegions:

  Logical value to indicate whether the blacklisted portion must be
  subtracted from the regions instead of removing them entirely. Notice
  that trimming collapses the regions overlapping each other within the
  same set. Default: `FALSE`.

- ignoreStrand:

  Logical value to indicate whether the strand must be ignored when
  computing the overlaps. Default: `TRUE`.

- emptySets:

  String indicating how to handle the sets left without any region, one
  among `"stop"`, `"remove"` or `"keep"`. Default: `"stop"`.

- verbose:

  Logical value to indicate whether the filtering messages must be
  printed. Default: `TRUE`.

## Value

An object of the same class as the input, with the blacklisted regions
removed. For a `RegionSetDE` object the blacklist and the filtering
counts are stored in the corresponding slots.

## See also

[`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md),
[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)
exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

regions <- splitLoadRegions(
  GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
  splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

# The example regions are shipped before the filtering, so this removes rows
regions <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)
regions
#> ### RegionSetDE object ###
#> Genome assembly:   rn4
#> Chromosome style:  UCSC
#> Region sets:       4
#> 
#>   promoterNonCpG  464 regions  (464,000 bp)
#>   intergenic      1,112 regions  (1,112,000 bp)
#>   geneBody        1,370 regions  (1,370,000 bp)
#>   promoterCpG     278 regions  (278,000 bp)
#> 
#> Blacklist:  applied (2,524 regions)
#> Whitelist:  not applied
#> 
#> Filtering steps: blacklist
#> (see the 'filtering.log' slot for the details)
```
