# contrastInfo

Describes the contrasts of a results object in one table: the engine
that tested them, the two groups compared and how many samples each
holds, and the distribution the test statistic follows with its degrees
of freedom. Together with the `stat` column of
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md)
it is what a power or sample size analysis needs, for instance with
power4peaks.

## Usage

``` r
contrastInfo(results)
```

## Arguments

- results:

  `RegionSetDE.results` object, or a `RegionSetDE.resultsList` holding
  several contrasts.

## Value

A data.frame with one row per contrast:

- `contrast`: the name of the contrast in the list, or its description
  for a single one.

- `contrast.description`: what the contrast compares.

- `engine`: `"edgeR"`, `"voom"`, `"limma"`, `"dream"` or `"deseq2"`.

- `column`, `group1`, `group2`: the column of the sample table the
  contrast separates and its two levels, the fold change being `group1`
  over `group2`.

- `n.group1`, `n.group2`: the number of samples in each group.

- `stat.distribution`: the distribution of `stat` under the null
  hypothesis, `"f"`, `"chisq"`, `"t"` or `"norm"`.

- `df1`, `df2`: its degrees of freedom, the median over the regions when
  they differ between regions, as the prior degrees of freedom of a
  robust fit make them do.

- `n.regions`: the number of rows tested.

The group columns are `NA` when the contrast is not a difference between
two levels of one column, a numeric vector over the coefficients for
instance, and `stat.distribution` is `NA` for the threshold tests run
with `lfcThreshold > 0`.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md),
[`contrastName`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastName.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

contrastInfo(results)
#>               contrast contrast.description engine    column group1 group2
#> 1 condition: SHR vs BN condition: SHR vs BN  edgeR condition    SHR     BN
#>   n.group1 n.group2 stat.distribution df1     df2 n.regions
#> 1        2        2                 f   1 18.9037      1895

# The statistics themselves, one per region
head(resultsTable(results)$stat)
#> [1] 0.166290413 2.337588301 0.010028735 2.084865356 0.007114484 2.273415136
```
