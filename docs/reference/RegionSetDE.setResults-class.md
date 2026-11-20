# RegionSetDE.setResults class

S4 class holding the outcome of a set level test. The table in the
`results` slot has one row per region set and carries the effect size
with its confidence interval next to the p-values, since with tens of
thousands of regions in a set the p-value stops discriminating long
before the effect size does.

## Usage

``` r
# S4 method for class 'RegionSetDE.setResults'
show(object)
```

## Arguments

- object:

  `RegionSetDE.setResults` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Slots

- `results`:

  `DataFrame` with one row per region set, or per pair of sets for a set
  contrast. Read it with
  [`resultsTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md).

- `regionStats`:

  `DataFrame` with the per-region statistics the set level test was
  computed on.

- `contrast`:

  String describing the contrast that was tested.

- `contrast.groups`:

  List with the `column` of the `colData` the contrast separates and the
  two `groups` it compares, the first one being the level the fold
  change is positive for. Empty when the contrast is not a difference
  between two levels of one variable.

- `test`:

  String with the kind of test, either `"set"` or `"setContrast"`.

- `methods`:

  Character vector with the tests that were run.

- `universe`:

  `RegionSetDE.universe` object with the rows every set was compared
  against. Empty when the test was self-contained only.

- `engine`:

  String with the engine that produced the per-region statistics.

- `thresholds`:

  List with the FDR cut-off and the correction method.

- `counts`:

  `RegionSetDE.counts` object the test was computed on, so that the
  signal behind a set level claim can be drawn without carrying a second
  object around. Empty when `carryCounts` was `FALSE`.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

setResults <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
setResults
#> An object of class 'RegionSetDE.setResults'
#>   test            : set 
#>   contrast        : condition: SHR vs BN 
#>   engine          : edgeR 
#>   methods         : camera, fry 
#>   universe        : otherSets (matched on width and abundance) 
#> 
#>      region.set n.regions mean.log2FC delta.log2FC CI.lower CI.upper CI.type
#>     promoterCpG       269     -0.8340      -0.9310   -1.930    0.561  sample
#>        geneBody       909      0.2910       0.4070   -0.333    0.919  sample
#>      intergenic       440      0.2650       0.2060   -0.171    0.367  sample
#>  promoterNonCpG       277     -0.0209      -0.0787   -0.263    0.251  sample
#>  camera.FDR fry.FDR inter.region.cor.universe
#>       0.938   0.479                     0.378
#>       0.938   0.479                     0.453
#>       0.938   0.669                     0.382
#>       0.938   0.669                     0.350
```
