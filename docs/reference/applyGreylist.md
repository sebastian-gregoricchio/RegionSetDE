# applyGreylist

Removes from every region set the regions overlapping a greylist,
typically the one built from the input libraries by
[`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md),
or one exported by another tool as a BED-like file. It works as
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
does, but the step is recorded as a greylist and the blacklist already
stored in the object is left as it is.

## Usage

``` r
applyGreylist(
  regionSet,
  greylist,
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

- greylist:

  `GRanges` returned by
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md),
  string indicating the path to a BED-like file, or data.frame with the
  regions to exclude.

- overlapType:

  String indicating the type of overlap required to greylist a region,
  one among `"any"`, `"within"`, `"start"`, `"end"` or `"equal"`.
  Default: `"any"`.

- minOverlapBp:

  Numeric value indicating the minimum number of bases that must overlap
  the greylist for a region to be removed. Default: `1`.

- minOverlapFraction:

  Numeric value between 0 and 1 indicating the minimum fraction of a
  region that must overlap the greylist for it to be removed. Default:
  `0`, any overlap is sufficient.

- trimRegions:

  Logical value to indicate whether the greylisted portion must be
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

An object of the same class as `regionSet`. For a `RegionSetDE` object
the step is added to the `filtering.log` as `"greylist"`, and the number
of greylisted regions and the bases they cover are stored in
`parameters$greylist`.

## Details

A greylist removes what the inputs flag as artefacts, and for broad
marks some of what it flags is genuine signal: heterochromatin marks
such as H3K9me3 sit on satellites and repeats, where inputs pile up too.
`trimRegions = TRUE` cuts the greylisted stretch out of a broad domain
and keeps the rest of it, and `minOverlapFraction` removes only the
regions mostly covered by the greylist. The messages report how many
regions each set keeps, which is the number to look at before going
further.

## See also

[`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md),
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
[`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)

regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
                            splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

# The exclusion list shipped with the package stands in for a greylist
greylist <- loadExampleData("exclusionRegions", verbose = FALSE)

greylisted <- applyGreylist(regions, greylist = greylist)
#> Applying the greylist (2524 regions):
#>   promoterNonCpG: 464/498 regions retained (6.8% removed)
#>   intergenic: 1112/1500 regions retained (25.9% removed)
#>   geneBody: 1370/1500 regions retained (8.7% removed)
#>   promoterCpG: 278/303 regions retained (8.3% removed)
filteringLog(greylisted)
#>       step     region.set n.before n.after n.removed
#> 1 greylist promoterNonCpG      498     464        34
#> 2 greylist     intergenic     1500    1112       388
#> 3 greylist       geneBody     1500    1370       130
#> 4 greylist    promoterCpG      303     278        25

# Broad domains lose the greylisted stretch only
trimmed <- applyGreylist(regions, greylist = greylist, trimRegions = TRUE, verbose = FALSE)
```
