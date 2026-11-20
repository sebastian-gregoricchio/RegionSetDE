# RegionSetDE class

S4 class collecting a group of genomic region sets together with the
filters applied to them. The regions are stored as a `GRangesList`, so
that any Bioconductor operation remains available through the `regions`
slot.

## Usage

``` r
# S4 method for class 'RegionSetDE'
show(object)
```

## Arguments

- object:

  A `RegionSetDE` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `regions`:

  `GRangesList` containing the region sets, one element per set.

- `consensus`:

  List with the consensus data when the regions were built from peaks by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md):
  the consensus of every group, the peaks of every sample after the
  exclusion, the total consensus, the table of the samples and the
  sample sheet. Empty for regions loaded any other way. Read it with
  [`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md).

## Author

Sebastian Gregoricchio

## Examples

``` r
regionTable <- loadExampleData("regions", verbose = FALSE)

regions <- loadRegions(split(regionTable[, c("seqnames", "start", "end")], regionTable$setName),
    genomeAssembly = "rn4", verbose = FALSE
)
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
#> Greylist:   not applied
#> Whitelist:  not applied
regionSetNames(regions)
#> [1] "geneBody"       "intergenic"     "promoterCpG"    "promoterNonCpG"
```
