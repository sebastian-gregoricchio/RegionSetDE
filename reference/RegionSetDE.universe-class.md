# RegionSetDE.universe class

S4 class holding, for every region set, the rows that make up the
universe of its competitive test: the set itself together with the rows
it is compared against. The positions refer to the object the universe
was built from, so the same universe serves every contrast run on that
object.

## Usage

``` r
# S4 method for class 'RegionSetDE.universe'
show(object)
```

## Arguments

- object:

  `RegionSetDE.universe` object.

## Value

[`show()`](https://rdrr.io/r/methods/show.html) prints a summary of the
object to the console.

## Details

The universe answers the question "more than what". A competitive test
asks whether the regions of a set respond more strongly than the rows
around them, and the answer depends entirely on which rows those are, so
the choice is stored rather than assumed. This has nothing to do with
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
which counts genome bins for the normalisation.

## Slots

- `index`:

  List with one element per region set, each a vector of row positions
  forming its universe.

- `type`:

  String with the way the universe was built, one of `"otherSets"` and
  `"supplied"`.

- `matching`:

  Character vector with the covariates the comparison rows were matched
  on.

- `comparison.sets`:

  Character vector with the names of the region sets the comparison rows
  were drawn from. A competitive p-value is relative to these and to
  nothing else, so they are stored beside the result rather than left to
  be reconstructed.

- `diagnostics`:

  `DataFrame` with, for every set, the median width and abundance of the
  set and of the rows it is compared against.

- `n.rows`:

  Numeric value with the number of rows of the object the positions
  refer to.

## See also

[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)

universe <- makeSetUniverse(fit, verbose = FALSE)
universe
#> An object of class 'RegionSetDE.universe'
#>   type            : otherSets 
#>   matched on      : width, abundance 
#>   compared against: promoterNonCpG, intergenic, geneBody, promoterCpG 
#>   sets            : 4 
#> 
#>      region.set n.regions n.comparison median.width median.width.comparison
#>  promoterNonCpG       277         1378         1000                    1000
#>      intergenic       440         1354         1000                    1000
#>        geneBody       909          986         1000                    1000
#>     promoterCpG       269          314         1000                    1000
#>  median.abundance median.abundance.comparison
#>              3.77                        3.78
#>              3.44                        3.76
#>              3.60                        4.08
#>              9.81                        3.75
```
