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
# The peaks of one sample of the AR example, cut into tiles of 100 bp
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)

tiledCounts <- countReads(peakRegions, sampleSheet = sampleSheet, tileWidth = 100, verbose = FALSE)
tiledCounts <- normalizeCounts(tiledCounts, method = "TMM", verbose = FALSE)
tiledFit <- fitRegions(tiledCounts, design = ~ condition, verbose = FALSE)

# The tiles are combined into regions, and their own statistics sit behind the combined ones
tiledResults <- testRegions(tiledFit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)
head(tileTable(tiledResults), 3)
#>                                       log2FC average.signal      stat
#> peaks|chr19:46089737-46089961|tile1 2.594558       9.218571 0.9729137
#> peaks|chr19:46089737-46089961|tile2 3.470645       9.472299 1.9880706
#> peaks|chr19:46089737-46089961|tile3 3.233674       9.376084 1.6438457
#>                                     stat.distribution df1  df2   p.value
#> peaks|chr19:46089737-46089961|tile1                 f   1 2712 0.3240426
#> peaks|chr19:46089737-46089961|tile2                 f   1 2712 0.1586574
#> peaks|chr19:46089737-46089961|tile3                 f   1 2712 0.1999082
#>                                     average.signal.DMSO average.signal.R1881_4h
#> peaks|chr19:46089737-46089961|tile1            8.822609                9.565452
#> peaks|chr19:46089737-46089961|tile2            8.822609                9.805591
#> peaks|chr19:46089737-46089961|tile3            8.822609                9.690221
#>                                     average.signal.R1881_24h region.set
#> peaks|chr19:46089737-46089961|tile1                 9.148300      peaks
#> peaks|chr19:46089737-46089961|tile2                 9.565067      peaks
#> peaks|chr19:46089737-46089961|tile3                 9.439105      peaks
#>                                                   region.id tile.id
#> peaks|chr19:46089737-46089961|tile1 chr19:46089737-46089961       1
#> peaks|chr19:46089737-46089961|tile2 chr19:46089737-46089961       2
#> peaks|chr19:46089737-46089961|tile3 chr19:46089737-46089961       3
#>                                                        region.key seqnames
#> peaks|chr19:46089737-46089961|tile1 peaks|chr19:46089737-46089961    chr19
#> peaks|chr19:46089737-46089961|tile2 peaks|chr19:46089737-46089961    chr19
#> peaks|chr19:46089737-46089961|tile3 peaks|chr19:46089737-46089961    chr19
#>                                        start      end width
#> peaks|chr19:46089737-46089961|tile1 46089737 46089836   100
#> peaks|chr19:46089737-46089961|tile2 46089837 46089936   100
#> peaks|chr19:46089737-46089961|tile3 46089937 46089961    25
#>                                                           name score      V7
#> peaks|chr19:46089737-46089961|tile1 AR_R1881_24h_r1_peak_21400   159 9.71586
#> peaks|chr19:46089737-46089961|tile2 AR_R1881_24h_r1_peak_21400   159 9.71586
#> peaks|chr19:46089737-46089961|tile3 AR_R1881_24h_r1_peak_21400   159 9.71586
#>                                          V8      V9 V10       FDR
#> peaks|chr19:46089737-46089961|tile1 18.4847 15.9682 103 0.5503016
#> peaks|chr19:46089737-46089961|tile2 18.4847 15.9682 103 0.3900078
#> peaks|chr19:46089737-46089961|tile3 18.4847 15.9682 103 0.4290677
```
