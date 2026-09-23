# loadBlacklist

Returns a blacklist shipped with the package, ready for
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
or for the `excludeRegions` argument of
[`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).
The files travel with the package, so nothing is downloaded and nothing
depends on a hub being reachable.

## Usage

``` r
loadBlacklist(
  genome,
  assay = NULL,
  source = NULL,
  seqlevelsStyle = "UCSC",
  verbose = TRUE
)
```

## Arguments

- genome:

  String with the genome assembly, such as `"hg38"`, `"mm10"` or `"hs1"`
  for T2T-CHM13v2.0. The usual aliases are understood, `"GRCh38"` and
  `"T2T"` for instance. Run
  [`availableRegionLists`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/availableRegionLists.md)
  for the whole list.

- assay:

  String with the assay the list was built for, `"cutrun"` or
  `"cuttag"`. Default: `NULL`, the ENCODE blacklist, which is not tied
  to an assay.

- source:

  String with the laboratory or project the list comes from, as
  [`availableRegionLists`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/availableRegionLists.md)
  prints it, needed only when a genome carries more than one list for
  the same assay. Default: `NULL`.

- seqlevelsStyle:

  String with the chromosome naming style of the output, one among
  `"UCSC"` (chr1), `"Ensembl"` (1) and `"NCBI"`, or `NULL` to keep the
  names of the file. Default: `"UCSC"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `GRanges` with the regions of the list and their name, carrying the
assembly in its `genome` and the source, the version and the reference
in its metadata.

## Details

Without an assay the list is the one meant for ChIP-seq and ATAC-seq,
the regions of anomalous coverage found across many experiments: the
ENCODE blacklist version 2 for the assemblies it covers, and for
T2T-CHM13v2.0, which ENCODE never covered, the set the excluderanges
authors built by running the same software on that assembly. Naming an
assay returns instead the high signal regions of the CUT&RUN greenlist
paper, built from hundreds of CUT&RUN or CUT&Tag libraries, one list per
assay. They are not interchangeable: the first kind is about the genome,
the other two are about what a given protocol does to it.

The files are the published ones, converted to gzipped BED and nothing
else. `inst/extdata/regionLists/SOURCES.md` records where each of them
comes from, the version, the date it was taken and the licence.

## See also

[`loadGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md),
[`availableRegionLists`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/availableRegionLists.md),
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
blacklist <- loadBlacklist("hg38")
#> The ENCODE blacklist v2 for hg38: 636 regions covering 227.2 Mb.
head(blacklist, 3)
#> GRanges object with 3 ranges and 1 metadata column:
#>       seqnames            ranges strand |               name
#>          <Rle>         <IRanges>  <Rle> |        <character>
#>   [1]    chr10           1-45700      * |    Low Mappability
#>   [2]    chr10 38481301-38596500      * | High Signal Region
#>   [3]    chr10 38782601-38967900      * | High Signal Region
#>   -------
#>   seqinfo: 24 sequences from hg38 genome; no seqlengths

# The regions the CUT&RUN protocol piles reads on, which the ENCODE list does not cover
cutrunBlacklist <- loadBlacklist("hg38", assay = "cutrun")
#> The deMello blacklist v1 for hg38: 832 regions covering 10.1 Mb.

# T2T-CHM13v2.0, under any of the names it goes by
t2tBlacklist <- loadBlacklist("T2T")
#> The excluderanges blacklist v1 for hs1: 3565 regions covering 275.5 Mb.
```
