# loadRegions

Imports a collection of genomic regions provided as BED-like files
and/or GRanges objects and returns them as a single `RegionSetDE`
object. Files, ranges and data.frames can be mixed within the same list.

## Usage

``` r
loadRegions(
  regions,
  regionNames = NULL,
  keepMetadata = TRUE,
  sortRegions = TRUE,
  reduceRegions = FALSE,
  removeDuplicatedRegions = TRUE,
  duplicatedSets = "stop",
  seqlevelsStyle = "UCSC",
  genomeAssembly = NULL,
  skipInvalid = FALSE,
  outputFormat = "RegionSetDE",
  verbose = TRUE
)
```

## Arguments

- regions:

  List in which each element is either a string indicating the path to a
  BED-like file (bed/narrowPeak/broadPeak, gzipped files allowed), a
  `GRanges` object, or a data.frame containing chromosome/start/end
  columns. A single path or a single `GRanges` is accepted as well.

- regionNames:

  Character vector with the names to assign to each element of the list.
  Default: `NULL`, meaning that the names of the input list are used
  and, when missing, the file base names.

- keepMetadata:

  Logical value to indicate whether the metadata columns must be kept.
  Default: `TRUE`.

- sortRegions:

  Logical value to indicate whether the regions must be sorted by
  coordinate. Default: `TRUE`.

- reduceRegions:

  Logical value to indicate whether the overlapping regions within the
  same element must be collapsed
  ([`IRanges::reduce`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)).
  Notice that metadata are lost upon reduction. Default: `FALSE`.

- removeDuplicatedRegions:

  Logical value to indicate whether the regions with identical
  coordinates within the same element must be collapsed to a single
  entry. Only the metadata of the first occurrence are kept. Ignored
  when `reduceRegions = TRUE`. Default: `TRUE`.

- duplicatedSets:

  String indicating how to handle two or more elements containing
  exactly the same regions: `"stop"` interrupts the loading, `"remove"`
  keeps only the first occurrence, `"keep"` disables the check. Notice
  that identical sets are perfectly correlated and inflate the number of
  tests in the multiple-testing correction. Default: `"stop"`.

- seqlevelsStyle:

  String indicating the chromosome naming style to apply to all the
  elements, one among `"UCSC"` (chr1), `"Ensembl"` (1) or `"NCBI"`. When
  set to `NULL` the names are left as they are and the loading is
  interrupted if the elements use different styles. Default: `"UCSC"`.

- genomeAssembly:

  String indicating the genome assembly to store in the ranges seqinfo,
  e.g. `"hg38"`. Default: `NULL`, no assembly is assigned.

- skipInvalid:

  Logical value to indicate whether the elements that cannot be loaded
  (missing files, malformed tables, unsupported classes) must be skipped
  with a warning rather than interrupting the loading. Default: `FALSE`.

- outputFormat:

  String indicating the class of the returned object, one among
  `"RegionSetDE"`, `"GRangesList"` or `"list"`. Default:
  `"RegionSetDE"`.

- verbose:

  Logical value to indicate whether the loading messages must be
  printed. Default: `TRUE`.

## Value

A `RegionSetDE` object or, depending on `outputFormat`, a named
`GRangesList` or a named list of `GRanges`.

## See also

[`splitLoadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitLoadRegions.md),
[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
[`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)

regionList <- split(regionTable[, c("seqnames", "start", "end")],
                    regionTable$setName)

regions <- loadRegions(regionList, genomeAssembly = "rn4", verbose = FALSE)
regions
#> ### RegionSetDE object ###
#> Genome assembly:   rn4
#> Chromosome style:  UCSC
#> Region sets:       4
#> 
#>   geneBody        1,500 regions  (1,500,000 bp)
#>   intergenic      1,500 regions  (1,500,000 bp)
#>   promoterCpG     303 regions  (303,000 bp)
#>   promoterNonCpG  498 regions  (498,000 bp)
#> 
#> Blacklist:  not applied
#> Whitelist:  not applied

if (FALSE) { # \dontrun{
regions <- loadRegions(list(promoters = "peaks/promoters.bed",
                            enhancers = grEnhancers,
                            "peaks/CTCF.narrowPeak"),
                       genomeAssembly = "hg38")
} # }
```
