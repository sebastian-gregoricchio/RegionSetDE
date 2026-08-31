# loadExampleData

Loads one of the example objects installed with RegionSetDE. The objects
cover the whole path from a collection of genomic regions to a fitted
model, and are the ones used throughout the vignette and the examples.

## Usage

``` r
loadExampleData(dataset = "counts", verbose = TRUE)
```

## Arguments

- dataset:

  String indicating which object to load, one among `"counts"`, `"fit"`,
  `"regions"`, `"exclusionRegions"`, `"sampleSheet"` and
  `"buildMetadata"`. The string `"blacklist"` is accepted as a synonym
  of `"exclusionRegions"`. Default: `"counts"`.

- verbose:

  Logical value to indicate whether the loading message must be printed.
  Default: `TRUE`.

## Value

The requested object. `"counts"` returns a `RegionSetDE.counts` object,
unnormalised and unfiltered, with the background bins already stored in
its metadata. `"fit"` returns a `RegionSetDE.fit` object ready for
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
and
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md).
`"regions"` returns a data.frame with one row per region and its set
membership, `"exclusionRegions"` a `GRanges`, `"sampleSheet"` a
data.frame describing every library of the source dataset, and
`"buildMetadata"` a list with the parameters used to generate the
others.

## Details

The example data comes from the liver ChIP-seq libraries of the
EURATRANS project, distributed by the `chromstaRData` package and
aligned to rn4. The contrast is H3K4me3 in the spontaneously
hypertensive (SHR) rat against the Brown Norway (BN) strain, two
biological replicates each, restricted to chromosome 12.

The regions are one kilobase windows split into four sets ordered by the
amount of H3K4me3 they are expected to carry: promoters overlapping a
CpG island, promoters without one, positions inside gene bodies away
from any transcription start site, and intergenic positions that serve
as the low-signal control. The sets are disjoint, and a window claimed
by an earlier set is never reused by a later one.

Only the counts are stored on disk. Asking for `"fit"` runs
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
and
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
on them, with the parameters recorded in `"buildMetadata"`, and keeps
the result for the rest of the session. A fitted model carries the
internals of its engine, and those change between releases of `edgeR`
and `limma`, so a serialised fit would break as soon as the engine moved
underneath it.

rn4 has no curated blacklist, so `"exclusionRegions"` is assembled from
the UCSC assembly gap track and from bins carrying implausible coverage
in the input libraries. It is not an ENCODE-grade exclusion list and
should not be reused outside these examples.

The script that generated the stored objects is installed with the
package, at
`system.file("scripts", "make-data.R", package = "RegionSetDE")`.

## See also

[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md),
[`splitLoadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitLoadRegions.md),
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts")
#> Loaded the 'counts' example object (RegionSetDE.counts).
counts
#> class: RegionSetDE.counts 
#> dim: 3224 4 
#> metadata(2): signal.type background
#> assays(1): counts
#> rownames(3224): promoterNonCpG|region_00002 promoterNonCpG|region_00003
#>   ... promoterCpG|region_03797 promoterCpG|region_03798
#> rowData names(3): region.set region.id tile.id
#> colnames(4): lv-H3K4me3-BN-female-bio1-tech1
#>   lv-H3K4me3-BN-male-bio2-tech1 lv-H3K4me3-SHR-male-bio2-tech1
#>   lv-H3K4me3-SHR-male-bio3-tech1
#> colData names(7): sample bam.file ... paired.end library.size

fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
topRegions(results, n = 5, FDR = 1)
#>       region.set    region.id seqnames    start      end width    log2FC
#> 1 promoterNonCpG region_02996    chr12 36842295 36843294  1000 -5.203835
#> 2     intergenic region_03590    chr12 44174500 44175499  1000 -2.887100
#> 3 promoterNonCpG region_00212    chr12  2500829  2501828  1000 -3.222630
#> 4       geneBody region_02435    chr12 29881730 29882729  1000 -2.778908
#> 5       geneBody region_02220    chr12 27481625 27482624  1000 -2.281658
#>   average.signal     stat      p.value          FDR diff.status
#> 1       5.159079 84.13617 1.148685e-08 2.176757e-05        down
#> 2       5.237329 43.04577 2.743053e-06 2.599042e-03        down
#> 3       4.816977 34.29977 1.370844e-05 6.573186e-03        down
#> 4       5.406565 36.15528 1.387480e-05 6.573186e-03        down
#> 5       5.277404 29.11255 3.347081e-05 1.268544e-02        down
```
