# testSetContrast

Asks whether a contrast affects one region set differently from another.
This is the comparison behind questions of the kind "does the treatment
reduce the mark more at Polycomb promoters than at active enhancers",
and it is the one claim that a global normalisation error cannot
manufacture, since a scaling factor that is wrong for one set is wrong
for the other in the same way.

## Usage

``` r
testSetContrast(
  fit,
  contrast,
  set1 = NULL,
  set2 = NULL,
  interRegionCor = NULL,
  useRanks = FALSE,
  sharedRegions = "drop",
  FDR = 0.05,
  adjustMethod = "BH",
  carryCounts = TRUE,
  verbose = TRUE
)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrast:

  Contrast to test, in the syntax accepted by
  [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
  or a named list of contrasts to run on the same fit.

- set1:

  Character vector with the name, or names, of the first region set.
  Default: `NULL`, every pair of sets is tested.

- set2:

  Character vector with the name, or names, of the second region set.
  Default: `NULL`.

- interRegionCor:

  Numeric value with the correlation between the regions of a set.
  Default: `NULL`, estimated from the residuals.

- useRanks:

  Logical value to indicate whether the test must work on the ranks
  rather than on the statistics. Default: `FALSE`.

- sharedRegions:

  String with what to do with the regions belonging to both sets, either
  `"drop"` or `"stop"`. Default: `"drop"`.

- FDR:

  Numeric value with the adjusted p-value cut-off reported in the
  output. Default: `0.05`.

- adjustMethod:

  String with the multiple testing correction across the pairs. Default:
  `"BH"`.

- carryCounts:

  Logical value to indicate whether the counts must travel inside the
  result. Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.setResults` object with one row per pair of sets, or a
`RegionSetDE.setResultsList` when `contrast` is a named list.

## Details

The test restricts the universe to the two sets and runs the competitive
test of
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
on the first of them, which is exactly a comparison of the first set
against the second. The effect size is the difference between the two
mean log2 fold changes, with a confidence interval carrying the variance
inflation of both sets. A region that belongs to both sets carries the
same reads into both sides of the comparison and pulls the difference
towards zero. Those regions are removed by default and the number
removed is reported; `sharedRegions = "stop"` refuses to run instead,
which is the safer setting when the overlap is unexpected.

## See also

[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`plotSetEffect`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

# Do the CpG island promoters respond differently from the intergenic control?
setContrast <- testSetContrast(fit,
                               contrast = c("condition", "SHR", "BN"),
                               set1 = "promoterCpG",
                               set2 = "intergenic",
                               verbose = FALSE)
setContrast
#> An object of class 'RegionSetDE.setResults'
#>   test            : setContrast 
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   methods         : camera 
#>   universe        : pairedSet  
#> 
#>        set.1      set.2 delta.log2FC CI.lower CI.upper camera.FDR
#>  promoterCpG intergenic         -1.1    -2.73    0.535      0.316
```
