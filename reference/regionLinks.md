# regionLinks

Returns the region-to-gene links stored in an object by
[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md),
one row per link. Also sets them, which is how a contact table computed
elsewhere is attached to an object that never went through a `TxDb`.

## Usage

``` r
regionLinks(object)

# S4 method for class 'RegionSetDE.provenance'
regionLinks(object)

regionLinks(object) <- value

# S4 method for class 'RegionSetDE.provenance'
regionLinks(object) <- value

# S4 method for class 'RegionSetDE.geneResults'
regionLinks(object)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit`, `RegionSetDE.results` or
  `RegionSetDE.setResults` object. Not a `RegionSetDE` object: the links
  live beside the counted rows, so the regions have to be counted first.

- value:

  Data.frame or `GRanges` with the links, for the replacement form.

## Value

A data.frame with the columns `region` and `gene.id`, together with
whatever else the links carried: `gene.symbol`, `distance.to.tss`,
`link.weight`, `weight.source` and `link.type`. `NULL` when the object
has not been annotated.

## Details

A region linked to no gene is kept in the table with a missing
`gene.id`, which is how an unassigned region stays distinguishable from
a region nothing was computed for. Count the genes with `na.omit` rather
than on the rows.

The replacement form runs the same checks
[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md)
does and refuses a table that would be ambiguous, so attaching links by
hand is as safe as building them.

## See also

[`annotateRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md),
[`geneSetsToRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

regionTable <- as.data.frame(SummarizedExperiment::rowRanges(counts))[1:4, 1:3]
regionLinks(counts) <- data.frame(regionTable, gene.id = c("11287", "66968", "73824", "20683"))

regionLinks(counts)
#>                                        region gene.id
#> promoterNonCpG|region_00002   chr12:4041-5040   11287
#> promoterNonCpG|region_00003   chr12:6212-7211   66968
#> promoterNonCpG|region_00005 chr12:10685-11684   73824
#> promoterNonCpG|region_00006 chr12:13540-14539   20683
```
