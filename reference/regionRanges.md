# regionRanges

Returns the region sets of a `RegionSetDE` object as a named
`GRangesList`, one element per set.

## Usage

``` r
regionRanges(object)

# S4 method for class 'RegionSetDE'
regionRanges(object)
```

## Arguments

- object:

  `RegionSetDE` object returned by
  [`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md),
  [`splitLoadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitLoadRegions.md)
  or
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).

## Value

A named `GRangesList` with the regions of every set, together with their
metadata columns.

## See also

[`regionSetNames`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionSetNames.md),
[`resultRanges`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultRanges.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)

regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
                            splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

lengths(regionRanges(regions))
#> promoterNonCpG     intergenic       geneBody    promoterCpG 
#>            498           1500           1500            303 
head(regionRanges(regions)$promoterCpG, 3)
#> GRanges object with 3 ranges and 1 metadata column:
#>       seqnames          ranges strand |     regionId
#>          <Rle>       <IRanges>  <Rle> |  <character>
#>   [1]    chr12 1019593-1020592      * | region_00090
#>   [2]    chr12 1291077-1292076      * | region_00114
#>   [3]    chr12 1365795-1366794      * | region_00119
#>   -------
#>   seqinfo: 1 sequence from rn4 genome; no seqlengths
```
