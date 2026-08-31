# filterRegions

Removes the rows of a `RegionSetDE.counts` object that carry too little
signal to say anything about a contrast. The decision is taken on the
average abundance alone, never on the variance or on a fold change, so
that the rows kept are independent of the comparison that will be run on
them afterwards.

## Usage

``` r
filterRegions(
  counts,
  method = "background",
  foldChange = 2,
  minCount = 10,
  proportion = 0.5,
  design = NULL,
  group = NULL,
  keep = NULL,
  byWidth = TRUE,
  referenceWidth = NULL,
  bySet = TRUE,
  widthStrata = 5,
  wholeRegion = FALSE,
  minTiles = 1,
  assay = "counts",
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- method:

  String with the criterion, one of `"background"`, `"abundance"`,
  `"proportion"`, `"byExpr"` and `"manual"`. Default: `"background"`.

- foldChange:

  Numeric value with the enrichment over the background a region must
  reach, on the linear scale. Only for `method = "background"`. Default:
  `2`.

- minCount:

  Numeric value with the number of reads a region of reference width
  must carry, on average across the samples. Only for
  `method = "abundance"`. Default: `10`.

- proportion:

  Numeric value between 0 and 1 with the fraction of rows to keep. Only
  for `method = "proportion"`. Default: `0.5`.

- design:

  Design matrix or formula passed to
  [`edgeR::filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).
  Only for `method = "byExpr"`. Default: `NULL`.

- group:

  String with the name of a `colData` column holding the experimental
  groups, passed to
  [`edgeR::filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).
  Only for `method = "byExpr"`. Default: `NULL`.

- keep:

  Logical vector with one value per row, or a vector of row positions.
  Only for `method = "manual"`. Default: `NULL`.

- byWidth:

  Logical value to indicate whether the abundance must be brought to a
  common width before being compared to the threshold. Default: `TRUE`.

- referenceWidth:

  Numeric value with the width the abundances are scaled to. Default:
  `NULL`, the median width of the rows.

- bySet:

  Logical value to indicate whether the threshold of
  `method = "proportion"` must be computed inside each region set rather
  than over the whole object. Default: `TRUE`.

- widthStrata:

  Numeric value with the number of width strata used by
  `method = "proportion"` when `byWidth = TRUE`. Default: `5`.

- wholeRegion:

  Logical value to indicate whether a tiled region must be kept in full
  as soon as one of its tiles passes. Ignored when the counts were not
  tiled. Default: `FALSE`.

- minTiles:

  Numeric value with the number of tiles a region must retain to
  survive. Ignored when the counts were not tiled. Default: `1`.

- assay:

  String with the name of the assay holding the values used to compute
  the abundance. Default: `"counts"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.counts` object with fewer rows.

## Details

The filter runs before
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
because the dispersion trend is fitted across the rows of the object:
leave in a few thousand regions that see two reads each and the trend
describes them rather than the regions being tested. It also feeds the
multiple testing correction, since every row that survives costs power
at the adjustment step.

Width is the reason a plain count threshold does not work on arbitrary
region sets. A 40 kb domain accumulates more reads than a 400 bp
promoter window at the same signal density, so a threshold in reads
keeps every broad region and drops every narrow one, whatever the
biology. With `byWidth = TRUE` the abundance of each row is divided by
its width and brought back to `referenceWidth`, and the threshold then
reads as "this many reads in a region of that width". This matters here
more than in a window based analysis, where every window has the same
size by construction.

The methods differ in what they compare against. `"background"` needs
the bins from
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)
and keeps the regions whose abundance exceeds, by `foldChange`, what a
stretch of genome of the same width would carry: an absolute statement,
in the sense that it does not depend on which other regions were loaded.
`"abundance"` is the simpler version of the same idea, with the
threshold given directly in reads. `"proportion"` keeps the strongest
fraction of the rows and is relative by construction, which is why
`bySet` defaults to `TRUE`: filtering the sets together lets a uniformly
weak set be removed entirely, and a set that no longer has rows cannot
be tested at the set level later. `"byExpr"` hands the decision to
[`edgeR::filterByExpr`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html),
which reads the group sizes from the design.

On a tiled object the filter applies to the tiles. A region that loses
every tile disappears, and the count of regions lost this way is
reported. `wholeRegion = TRUE` keeps all the tiles of a region as soon
as one of them passes, which preserves the span of the region at the
cost of carrying the empty tiles through the fit; the Simes combination
in
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
pays for those tiles in multiplicity, so the default leaves them out.

## See also

[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

filtered <- filterRegions(counts, foldChange = 2, verbose = FALSE)

nrow(counts)
#> [1] 3224
nrow(filtered)
#> [1] 1895

# The intergenic control set is the one that loses most of its rows
table(SummarizedExperiment::rowData(filtered)$region.set)
#> 
#>       geneBody     intergenic    promoterCpG promoterNonCpG 
#>            909            440            269            277 
```
