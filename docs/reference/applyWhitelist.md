# applyWhitelist

Restricts every region set to the regions overlapping a whitelist, for
instance a set of accessible or mappable regions. The whitelist can be
provided as a BED-like file, a `GRanges` or a data.frame.

## Usage

``` r
applyWhitelist(
  regionSet,
  whitelist,
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

- whitelist:

  String indicating the path to a BED-like file, a `GRanges` or a
  data.frame with the regions to retain.

- overlapType:

  String indicating the type of overlap required to retain a region, one
  among `"any"`, `"within"`, `"start"`, `"end"` or `"equal"`. Default:
  `"any"`.

- minOverlapBp:

  Numeric value indicating the minimum number of bases that must overlap
  the whitelist for a region to be retained. Default: `1`.

- minOverlapFraction:

  Numeric value between 0 and 1 indicating the minimum fraction of a
  region that must overlap the whitelist for it to be retained. Default:
  `0`, any overlap is sufficient.

- trimRegions:

  Logical value to indicate whether the regions must be clipped at the
  whitelist boundaries instead of being retained entirely. A region
  spanning two whitelisted blocks is split accordingly. Default:
  `FALSE`.

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

An object of the same class as the input, restricted to the whitelisted
regions. For a `RegionSetDE` object the whitelist and the filtering
counts are stored in the corresponding slots.

## See also

[`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
[`loadRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)

regions <- splitLoadRegions(
  GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
  splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

# Restrict the analysis to the first half of the chromosome
whitelist <- GenomicRanges::GRanges(
  seqnames = "chr12",
  ranges = IRanges::IRanges(start = 1, end = 25e6))

regions <- applyWhitelist(regions, whitelist = whitelist, verbose = FALSE)
regions
#> ### RegionSetDE object ###
#> Genome assembly:   rn4
#> Chromosome style:  UCSC
#> Region sets:       4
#> 
#>   promoterNonCpG  319 regions  (319,000 bp)
#>   intergenic      781 regions  (781,000 bp)
#>   geneBody        768 regions  (768,000 bp)
#>   promoterCpG     160 regions  (160,000 bp)
#> 
#> Blacklist:  not applied
#> Whitelist:  applied (1 regions)
#> 
#> Filtering steps: whitelist
#> (see the 'filtering.log' slot for the details)
```
