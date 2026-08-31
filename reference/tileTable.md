# tileTable

Returns the per-tile table of a `RegionSetDE.results` object, empty when
the counts were not tiled.

## Usage

``` r
tileTable(results)

# S4 method for class 'RegionSetDE.results'
tileTable(results)

# S4 method for class 'RegionSetDE.resultsList'
tileTable(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A data.frame with one row per tile.

## Author

Sebastian Gregoricchio

## Examples

``` r
# Tiling splits each region into fixed-width windows, kept alongside the region
bamFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")

exampleRegions <- GenomicRanges::GRanges(
  seqnames = rep(c("seq1", "seq2"), each = 2),
  ranges = IRanges::IRanges(start = rep(c(1, 800), 2), width = 400))

exampleRegions$setName <- rep(c("firstSet", "secondSet"), each = 2)

exampleSets <- splitLoadRegions(exampleRegions, splitBy = "setName",
                                seqlevelsStyle = NULL, verbose = FALSE)

tiledCounts <- countReads(exampleSets, bamFiles = bamFile,
                          sampleNames = "example", tileWidth = 100,
                          verbose = FALSE)

head(SummarizedExperiment::rowData(tiledCounts))
#> DataFrame with 6 rows and 3 columns
#>                               region.set     region.id   tile.id
#>                              <character>   <character> <integer>
#> firstSet|seq1:1-400|tile1       firstSet    seq1:1-400         1
#> firstSet|seq1:1-400|tile2       firstSet    seq1:1-400         2
#> firstSet|seq1:1-400|tile3       firstSet    seq1:1-400         3
#> firstSet|seq1:1-400|tile4       firstSet    seq1:1-400         4
#> firstSet|seq1:800-1199|tile1    firstSet seq1:800-1199         1
#> firstSet|seq1:800-1199|tile2    firstSet seq1:800-1199         2

if (FALSE) { # \dontrun{
# After testing, the per-tile statistics sit behind the combined ones
tiledResults <- testRegions(tiledFit, contrast = c("condition", "SHR", "BN"))
head(tileTable(tiledResults))
} # }
```
