# topRegions

Returns the regions that respond most strongly to a contrast, optionally
restricted to one region set or to one direction.

## Usage

``` r
topRegions(
  results,
  n = Inf,
  set = NULL,
  contrast = NULL,
  sortBy = "FDR",
  direction = "both",
  FDR = NULL,
  log2FC = NULL,
  level = "region"
)
```

## Arguments

- results:

  `RegionSetDE.results` or `RegionSetDE.resultsList` object.

- n:

  Numeric value with the number of regions to return. When `set` holds
  more than one name, `n` regions are returned for each of them.
  Default: `Inf`, every region passing the thresholds.

- set:

  Character vector with the names of the region sets to consider.
  Default: `NULL`, all of them pooled together.

- contrast:

  String with the name of the contrast to take, or its position, when
  `results` holds several of them. Default: `NULL`.

- sortBy:

  String with the column driving the ranking, one of `"FDR"`,
  `"p.value"`, `"log2FC"` and `"stat"`. Default: `"FDR"`.

- direction:

  String restricting the output to one direction of change, one of
  `"both"`, `"up"` and `"down"`. Default: `"both"`.

- FDR:

  Numeric value with an adjusted p-value cut-off applied before the
  ranking. Default: `NULL`, the threshold stored in the object.

- log2FC:

  Numeric value with an absolute log2 fold change cut-off applied before
  the ranking. Default: `NULL`, the threshold stored in the object.

- level:

  String indicating whether the regions (`"region"`) or the tiles
  (`"tile"`) must be ranked. Default: `"region"`.

## Value

A data.frame with the selected rows.

## Details

Ranking by `"log2FC"` sorts on the effect size alone and returns
whatever passes the cut-offs, which on a small object is often a handful
of low-count regions with a large and badly estimated fold change. Keep
an FDR cut-off in place when doing so.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`plotVolcano`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

# By default only the regions passing the thresholds come back
topRegions(results, n = 10)
#>        region.set    region.id seqnames    start      end width    log2FC
#> 1  promoterNonCpG region_02996    chr12 36842295 36843294  1000 -5.203835
#> 2      intergenic region_03590    chr12 44174500 44175499  1000 -2.887100
#> 3  promoterNonCpG region_00212    chr12  2500829  2501828  1000 -3.222630
#> 4        geneBody region_02435    chr12 29881730 29882729  1000 -2.778908
#> 5        geneBody region_02220    chr12 27481625 27482624  1000 -2.281658
#> 6      intergenic region_03406    chr12 42124500 42125499  1000  3.180337
#> 7      intergenic region_02700    chr12 33474500 33475499  1000 -2.202950
#> 8        geneBody region_00572    chr12  7295412  7296411  1000 -2.343642
#> 9  promoterNonCpG region_01472    chr12 17981162 17982161  1000  2.864652
#> 10    promoterCpG region_01273    chr12 15719347 15720346  1000 -1.771521
#>    average.signal     stat      p.value          FDR diff.status
#> 1        5.159079 84.13617 1.148685e-08 2.176757e-05        down
#> 2        5.237329 43.04577 2.743053e-06 2.599042e-03        down
#> 3        4.816977 34.29977 1.370844e-05 6.573186e-03        down
#> 4        5.406565 36.15528 1.387480e-05 6.573186e-03        down
#> 5        5.277404 29.11255 3.347081e-05 1.268544e-02        down
#> 6        4.777653 26.92056 4.205306e-05 1.328176e-02          up
#> 7        5.128097 25.79674 6.712520e-05 1.683932e-02        down
#> 8        6.847950 25.58223 7.108947e-05 1.683932e-02        down
#> 9        5.445639 23.21950 1.531036e-04 3.223681e-02          up
#> 10       5.933614 21.00985 2.058661e-04 3.901163e-02        down

# FDR = 1 ranks everything instead, which is what you want when power is low
topRegions(results, n = 10, FDR = 1)
#>        region.set    region.id seqnames    start      end width    log2FC
#> 1  promoterNonCpG region_02996    chr12 36842295 36843294  1000 -5.203835
#> 2      intergenic region_03590    chr12 44174500 44175499  1000 -2.887100
#> 3  promoterNonCpG region_00212    chr12  2500829  2501828  1000 -3.222630
#> 4        geneBody region_02435    chr12 29881730 29882729  1000 -2.778908
#> 5        geneBody region_02220    chr12 27481625 27482624  1000 -2.281658
#> 6      intergenic region_03406    chr12 42124500 42125499  1000  3.180337
#> 7      intergenic region_02700    chr12 33474500 33475499  1000 -2.202950
#> 8        geneBody region_00572    chr12  7295412  7296411  1000 -2.343642
#> 9  promoterNonCpG region_01472    chr12 17981162 17982161  1000  2.864652
#> 10    promoterCpG region_01273    chr12 15719347 15720346  1000 -1.771521
#>    average.signal     stat      p.value          FDR diff.status
#> 1        5.159079 84.13617 1.148685e-08 2.176757e-05        down
#> 2        5.237329 43.04577 2.743053e-06 2.599042e-03        down
#> 3        4.816977 34.29977 1.370844e-05 6.573186e-03        down
#> 4        5.406565 36.15528 1.387480e-05 6.573186e-03        down
#> 5        5.277404 29.11255 3.347081e-05 1.268544e-02        down
#> 6        4.777653 26.92056 4.205306e-05 1.328176e-02          up
#> 7        5.128097 25.79674 6.712520e-05 1.683932e-02        down
#> 8        6.847950 25.58223 7.108947e-05 1.683932e-02        down
#> 9        5.445639 23.21950 1.531036e-04 3.223681e-02          up
#> 10       5.933614 21.00985 2.058661e-04 3.901163e-02        down

# One set at a time, ranked rather than filtered
topRegions(results, n = 5, set = "promoterCpG", FDR = 1)
#>    region.set    region.id seqnames    start      end width    log2FC
#> 1 promoterCpG region_01273    chr12 15719347 15720346  1000 -1.771521
#> 2 promoterCpG region_02452    chr12 30097777 30098776  1000  1.513446
#> 3 promoterCpG region_01859    chr12 22815633 22816632  1000  1.445776
#> 4 promoterCpG region_00824    chr12 10369814 10370813  1000  1.635046
#> 5 promoterCpG region_00836    chr12 10545505 10546504  1000 -1.223234
#>   average.signal      stat      p.value        FDR diff.status
#> 1       5.933614 21.009851 0.0002058661 0.03901163        down
#> 2       5.117243 10.076467 0.0049984001 0.28702934        null
#> 3       5.016526  8.582683 0.0086069702 0.36371630        null
#> 4       4.658080  7.080136 0.0163032446 0.45433307        null
#> 5       5.879394  6.508768 0.0207077166 0.49672308        null

# Sorting by effect size instead of significance
topRegions(results, n = 5, FDR = 1, sortBy = "log2FC", direction = "down")
#>       region.set    region.id seqnames    start      end width    log2FC
#> 1 promoterNonCpG region_02996    chr12 36842295 36843294  1000 -5.203835
#> 2 promoterNonCpG region_00212    chr12  2500829  2501828  1000 -3.222630
#> 3     intergenic region_01115    chr12 13844500 13845499  1000 -3.154327
#> 4       geneBody region_03415    chr12 42209041 42210040  1000 -3.020137
#> 5     intergenic region_03590    chr12 44174500 44175499  1000 -2.887100
#>   average.signal      stat      p.value          FDR diff.status
#> 1       5.159079 84.136173 1.148685e-08 2.176757e-05        down
#> 2       4.816977 34.299774 1.370844e-05 6.573186e-03        down
#> 3       3.819372 15.362838 9.731331e-04 9.705722e-02        null
#> 4       3.189783  8.417927 8.828997e-03 3.637163e-01        null
#> 5       5.237329 43.045774 2.743053e-06 2.599042e-03        down
```
