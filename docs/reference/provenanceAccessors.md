# Provenance accessors

Return what every object of the package carries about the origin of its
regions: the lists of regions applied to them, the chromosome naming
style, and the parameters of every step that built the object. The
record travels from the regions into the counts, the fit and the
results, so the same accessors answer at every step.

## Usage

``` r
blacklistRanges(object)

# S4 method for class 'RegionSetDE.provenance'
blacklistRanges(object)

greylistRanges(object)

# S4 method for class 'RegionSetDE.provenance'
greylistRanges(object)

whitelistRanges(object)

# S4 method for class 'RegionSetDE.provenance'
whitelistRanges(object)

chromosomeStyle(object)

# S4 method for class 'RegionSetDE.provenance'
chromosomeStyle(object)

parameterLog(object)

# S4 method for class 'RegionSetDE.provenance'
parameterLog(object)
```

## Arguments

- object:

  Any object of the package: `RegionSetDE`, `RegionSetDE.counts`,
  `RegionSetDE.fit`, `RegionSetDE.results`, `RegionSetDE.setResults`,
  `RegionSetDE.setScores` or `RegionSetDE.geneResults`.

## Value

`blacklistRanges`, `greylistRanges` and `whitelistRanges` return a
`GRanges` with the list applied by
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
[`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md)
or
[`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md),
or by the matching argument of
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
and
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
and `NULL` when none was applied. `chromosomeStyle` returns a string
with the naming style of the chromosomes, such as `"UCSC"`, or `NA` when
the names were kept as they came. `parameterLog` returns a named list
with one element per step, named after the function, holding the
arguments it was called with.

## See also

[`filteringLog`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filteringLog.md),
[`genomeAssembly`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/genomeAssembly.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)
exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
    splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE
)
regions <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)

blacklistRanges(regions)
#> GRanges object with 2524 ranges and 0 metadata columns:
#>          seqnames            ranges strand
#>             <Rle>         <IRanges>  <Rle>
#>      [1]    chr12         2685-3206      *
#>      [2]    chr12         5789-5838      *
#>      [3]    chr12         9368-9417      *
#>      [4]    chr12       16603-17124      *
#>      [5]    chr12       23627-23676      *
#>      ...      ...               ...    ...
#>   [2520]    chr12 46653064-46653167      *
#>   [2521]    chr12 46676570-46676619      *
#>   [2522]    chr12 46698898-46698947      *
#>   [2523]    chr12 46755155-46756128      *
#>   [2524]    chr12 46757592-46759787      *
#>   -------
#>   seqinfo: 1 sequence from rn4 genome
is.null(greylistRanges(regions))
#> [1] TRUE
chromosomeStyle(regions)
#> [1] "UCSC"

# The steps that built the example counts
counts <- loadExampleData("counts", verbose = FALSE)
names(parameterLog(counts))
#> [1] "loadRegions"     "blacklist"       "countReads"      "countBackground"
```
