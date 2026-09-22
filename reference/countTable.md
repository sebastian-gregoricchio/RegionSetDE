# countTable

Returns the values of a counts object as a table, one row per region or
per tile, with the coordinates and the annotation of each row next to
the values of the samples. The values can be raw or normalised, and come
in wide or long format, or as a plain matrix. On a tiled object the
tiles can be combined back into their region.

## Usage

``` r
countTable(object, ...)

# S4 method for class 'RegionSetDE.counts'
countTable(
  object,
  level = "region",
  normalized = FALSE,
  format = "wide",
  set = NULL,
  tileSummary = NULL,
  extraColumns = TRUE,
  verbose = TRUE
)

# S4 method for class 'RegionSetDE.fit'
countTable(object, ...)

# S4 method for class 'RegionSetDE.results'
countTable(object, ...)

# S4 method for class 'RegionSetDE.setResults'
countTable(object, ...)

# S4 method for class 'RegionSetDE.resultsList'
countTable(object, ...)

# S4 method for class 'RegionSetDE.setResultsList'
countTable(object, ...)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit`, `RegionSetDE.results` or
  `RegionSetDE.setResults` object, or either of the two list classes
  holding several contrasts. For the results the counts are the ones
  carried inside them, which requires the test to have been run with
  `carryCounts = TRUE`.

- ...:

  Arguments passed on to the method for `RegionSetDE.counts` objects.

- level:

  String indicating whether the table must have one row per region
  (`"region"`) or one row per tile (`"tile"`). On a tiled object
  `"region"` combines the tiles of each region into a single row.
  Default: `"region"`.

- normalized:

  Logical value to indicate whether the normalised values must be
  returned instead of the raw ones. Default: `FALSE`.

- format:

  String with the shape of the output: `"wide"` gives one column per
  sample, `"long"` one row per row of the object and sample, with the
  `colData` of the samples attached, and `"matrix"` a numeric matrix
  with the rows named after the regions. Default: `"wide"`.

- set:

  Character vector with the names of the region sets to keep. Default:
  `NULL`, all of them.

- tileSummary:

  String indicating how the tiles are combined into their region when
  `level = "region"`, one among `"sum"`, `"mean"`, `"max"` and `"min"`.
  Default: `NULL`, read from the way the object was counted.

- extraColumns:

  Annotation carried by the regions that must be added to the table,
  after the coordinates. Either `TRUE` for every column of the `rowData`
  beyond the ones the package writes itself, `FALSE` for none, or a
  character vector naming the ones wanted. Ignored when
  `format = "matrix"`. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

With `format = "wide"`, a data.frame with one row per region, or per
tile, described by the `region.set`, `region.id`, `tile.id` (tiles
only), `seqnames`, `start`, `end` and `width` columns, followed by the
annotation of the regions and by one column per sample. When the tiles
have been combined, `n.tiles` reports how many of them each region was
built from. With `format = "long"`, a data.frame with the same
description repeated for every sample, a `sample` column, the values in
a column named after the assay they were read from (`counts`, or
`norm.counts` for the normalised ones), and the `colData` of the
samples. With `format = "matrix"`, a numeric matrix with one column per
sample and the rows named `"set|id"`, or `"set|id|tileN"` for the tiles.

## Details

With `level = "region"` on a tiled object the tiles of each region are
combined into one row, spanning from the first to the last tile present.
Read counts are summed. Signal read from bigWig files follows the
`summaryFunction` used by
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md):
a sum stays a sum, a mean is averaged with the width of each tile as its
weight, so that a shorter trailing tile counts for the bases it covers,
and maxima and minima stay maxima and minima. `tileSummary` overrides
the rule, for instance for a matrix imported through
[`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)
that holds a mean signal rather than counts.

For bigWig signal the combination is exact, since the values are
integrated base by base. The one exception is a mean computed with
`missingAsZero = FALSE`, which is taken over the covered bases only,
while the tiles are weighted by their full width. For reads the
combination is not exact.
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
counts a fragment in every row it overlaps, so a fragment lying across
the border between two tiles is counted in both, and the sum over the
tiles is higher than the count of the same region taken whole. The
excess grows with the fragment length relative to the tile width: a few
percent with tiles of several kb, about double when the tiles are as
narrow as the fragments. It is similar across libraries with similar
fragment lengths, so the values still compare between samples, but the
number of fragments falling in a region can only come from counting the
regions without tiles.

[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
combines the tiles differently. It averages them
(`tileHandling = "collapse"`) so that every region weighs the same
within its set, whatever its width. A table of counts describes the
regions themselves, and the sum is the value closest to a region counted
in one piece.

The normalised assay holds the raw counts divided by the scaling factor
of each sample, so the region values built from it equal the combined
raw counts divided by the same factor. The exception is
`method = "loess"`, where every row carries its own offset: each tile is
corrected at its own abundance before the tiles are combined, which is
not what a loess fit on the region counts would return.

After
[`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
a region holds only the tiles that passed the filter, and its value
covers those tiles alone.

## See also

[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`resultCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultCounts.md),
[`fitCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitCounts.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

# Raw counts, one column per sample
head(countTable(counts), 3)
#>                                 region.set    region.id seqnames start   end
#> promoterNonCpG|region_00002 promoterNonCpG region_00002    chr12  4041  5040
#> promoterNonCpG|region_00003 promoterNonCpG region_00003    chr12  6212  7211
#> promoterNonCpG|region_00005 promoterNonCpG region_00005    chr12 10685 11684
#>                             width     regionId lv-H3K4me3-BN-female-bio1-tech1
#> promoterNonCpG|region_00002  1000 region_00002                               0
#> promoterNonCpG|region_00003  1000 region_00003                               0
#> promoterNonCpG|region_00005  1000 region_00005                               1
#>                             lv-H3K4me3-BN-male-bio2-tech1
#> promoterNonCpG|region_00002                             0
#> promoterNonCpG|region_00003                             0
#> promoterNonCpG|region_00005                             0
#>                             lv-H3K4me3-SHR-male-bio2-tech1
#> promoterNonCpG|region_00002                              0
#> promoterNonCpG|region_00003                              0
#> promoterNonCpG|region_00005                              0
#>                             lv-H3K4me3-SHR-male-bio3-tech1
#> promoterNonCpG|region_00002                              0
#> promoterNonCpG|region_00003                              0
#> promoterNonCpG|region_00005                              1

# Normalised values in long format, with the sample annotation attached for ggplot2
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes
longTable <- countTable(counts, normalized = TRUE, format = "long")
head(longTable, 3)
#>       region.set    region.id seqnames start   end width     regionId
#> 1 promoterNonCpG region_00002    chr12  4041  5040  1000 region_00002
#> 2 promoterNonCpG region_00003    chr12  6212  7211  1000 region_00003
#> 3 promoterNonCpG region_00005    chr12 10685 11684  1000 region_00005
#>                            sample norm.counts
#> 1 lv-H3K4me3-BN-female-bio1-tech1    0.000000
#> 2 lv-H3K4me3-BN-female-bio1-tech1    0.000000
#> 3 lv-H3K4me3-BN-female-bio1-tech1    1.840561
#>                                                                                                                     bam.file
#> 1 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#> 2 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#> 3 /home/s.gregoricchio/R/x86_64-pc-linux-gnu-library/4.6/chromstaRData/extdata/euratrans/lv-H3K4me3-BN-female-bio1-tech1.bam
#>   condition    sex biologicalReplicate paired.end library.size norm.factor
#> 1        BN female                bio1      FALSE       386378   0.7959318
#> 2        BN female                bio1      FALSE       386378   0.7959318
#> 3        BN female                bio1      FALSE       386378   0.7959318
#>   scaling.factor
#> 1      0.5433126
#> 2      0.5433126
#> 3      0.5433126

# A plain matrix of one set, for ComplexHeatmap or any other tool
countMatrix <- countTable(counts, normalized = TRUE, format = "matrix", set = "promoterCpG")
dim(countMatrix)
#> [1] 278   4

# On a tiled object the tiles can be kept, or summed back into their region
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

countTable(tiledCounts, level = "tile")
#>                               region.set     region.id tile.id seqnames start
#> firstSet|seq1:1-400|tile1       firstSet    seq1:1-400       1     seq1     1
#> firstSet|seq1:1-400|tile2       firstSet    seq1:1-400       2     seq1   101
#> firstSet|seq1:1-400|tile3       firstSet    seq1:1-400       3     seq1   201
#> firstSet|seq1:1-400|tile4       firstSet    seq1:1-400       4     seq1   301
#> firstSet|seq1:800-1199|tile1    firstSet seq1:800-1199       1     seq1   800
#> firstSet|seq1:800-1199|tile2    firstSet seq1:800-1199       2     seq1   900
#> firstSet|seq1:800-1199|tile3    firstSet seq1:800-1199       3     seq1  1000
#> firstSet|seq1:800-1199|tile4    firstSet seq1:800-1199       4     seq1  1100
#> secondSet|seq2:1-400|tile1     secondSet    seq2:1-400       1     seq2     1
#> secondSet|seq2:1-400|tile2     secondSet    seq2:1-400       2     seq2   101
#> secondSet|seq2:1-400|tile3     secondSet    seq2:1-400       3     seq2   201
#> secondSet|seq2:1-400|tile4     secondSet    seq2:1-400       4     seq2   301
#> secondSet|seq2:800-1199|tile1  secondSet seq2:800-1199       1     seq2   800
#> secondSet|seq2:800-1199|tile2  secondSet seq2:800-1199       2     seq2   900
#> secondSet|seq2:800-1199|tile3  secondSet seq2:800-1199       3     seq2  1000
#> secondSet|seq2:800-1199|tile4  secondSet seq2:800-1199       4     seq2  1100
#>                                end width example
#> firstSet|seq1:1-400|tile1      100   100      20
#> firstSet|seq1:1-400|tile2      200   100      58
#> firstSet|seq1:1-400|tile3      300   100     115
#> firstSet|seq1:1-400|tile4      400   100     154
#> firstSet|seq1:800-1199|tile1   899   100     151
#> firstSet|seq1:800-1199|tile2   999   100     164
#> firstSet|seq1:800-1199|tile3  1099   100     185
#> firstSet|seq1:800-1199|tile4  1199   100     181
#> secondSet|seq2:1-400|tile1     100   100      58
#> secondSet|seq2:1-400|tile2     200   100     105
#> secondSet|seq2:1-400|tile3     300   100     178
#> secondSet|seq2:1-400|tile4     400   100     200
#> secondSet|seq2:800-1199|tile1  899   100     178
#> secondSet|seq2:800-1199|tile2  999   100     198
#> secondSet|seq2:800-1199|tile3 1099   100     206
#> secondSet|seq2:800-1199|tile4 1199   100     213
countTable(tiledCounts, level = "region")
#> The tiles have been summed into their regions. A fragment overlapping two tiles is counted in both, so these values are higher than the counts of the same regions taken whole.
#>                         region.set     region.id seqnames start  end width
#> firstSet|seq1:1-400       firstSet    seq1:1-400     seq1     1  400   400
#> firstSet|seq1:800-1199    firstSet seq1:800-1199     seq1   800 1199   400
#> secondSet|seq2:1-400     secondSet    seq2:1-400     seq2     1  400   400
#> secondSet|seq2:800-1199  secondSet seq2:800-1199     seq2   800 1199   400
#>                         n.tiles example
#> firstSet|seq1:1-400           4     347
#> firstSet|seq1:800-1199        4     681
#> secondSet|seq2:1-400          4     541
#> secondSet|seq2:800-1199       4     795
```
