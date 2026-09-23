# filteringLog

Returns the record of the filters applied to a RegionSetDE object, with
the number of regions before and after each step. The record travels
with the regions into the counts, the fit and the results.

## Usage

``` r
filteringLog(object)

# S4 method for class 'RegionSetDE.provenance'
filteringLog(object)
```

## Arguments

- object:

  Any object of the package: `RegionSetDE`, `RegionSetDE.counts`,
  `RegionSetDE.fit`, `RegionSetDE.results`, `RegionSetDE.setResults` or
  `RegionSetDE.setScores`.

## Value

A data.frame with one row per filtering step, empty when no filter has
been applied.

## See also

[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
[`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md),
[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)
exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
                            splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
regions <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)

filteringLog(regions)
#>        step     region.set n.before n.after n.removed
#> 1 blacklist promoterNonCpG      498     464        34
#> 2 blacklist     intergenic     1500    1112       388
#> 3 blacklist       geneBody     1500    1370       130
#> 4 blacklist    promoterCpG      303     278        25
```
