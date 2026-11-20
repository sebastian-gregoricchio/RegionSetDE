# RegionSetDE.provenance class

Virtual class collecting the filters and the parameters shared by all
the RegionSetDE objects, so that the origin of the regions survives
every step of the analysis. Not meant to be instantiated directly.

## Slots

- `blacklist`:

  `GRanges` with the regions removed by
  [`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
  `NULL` when no blacklist has been applied.

- `greylist`:

  `GRanges` with the regions removed by
  [`applyGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md),
  by the `greylist` argument of
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  or by
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md),
  `NULL` when no greylist has been applied. Several greylists applied
  one after the other are merged.

- `whitelist`:

  `GRanges` with the regions used by
  [`applyWhitelist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md)
  to restrict the sets, `NULL` when no whitelist has been applied.

- `genome.assembly`:

  String with the genome assembly of the regions, `NULL` when not
  declared.

- `seqlevels.style`:

  String with the chromosome naming style shared by all the sets.

- `filtering.log`:

  `DataFrame` collecting the number of regions before and after each
  filtering step. Read it with
  [`filteringLog`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filteringLog.md).

- `parameters`:

  List with the arguments used at each step.

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

# Every object of the package carries the provenance slots
is(counts, "RegionSetDE.provenance")
#> [1] TRUE
genomeAssembly(counts)
#> [1] "rn4"
filteringLog(counts)
#>        step     region.set n.before n.after n.removed
#> 1 blacklist promoterNonCpG      498     464        34
#> 2 blacklist     intergenic     1500    1112       388
#> 3 blacklist       geneBody     1500    1370       130
#> 4 blacklist    promoterCpG      303     278        25
```
