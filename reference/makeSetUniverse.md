# makeSetUniverse

Builds, for every region set, the universe of its competitive test: the
set itself together with the rows it will be compared against. The
comparison rows come from the other sets of the object, or from an index
of your own, and they can be matched to each set on width and on
baseline abundance so that the comparison is not driven by the sets
simply being made of different kinds of intervals.

Calling this is optional.
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
builds a matched universe and keeps it in the fit, and `matchOn` and
`universeRatio` cover the common adjustments there. Come here to reach
`strata`, or to supply an index of your own through `type = "supplied"`.

## Usage

``` r
makeSetUniverse(
  object,
  type = "otherSets",
  match = c("width", "abundance"),
  ratio = 5,
  strata = 5,
  regionSets = NULL,
  index = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.fit` or `RegionSetDE.counts` object.

- type:

  String with the source of the comparison rows, either `"otherSets"`
  (the rows of every other set) or `"supplied"`. Default: `"otherSets"`.

- match:

  Character vector with the covariates the comparison rows are matched
  on, among `"width"` and `"abundance"`. An empty vector takes every
  other row as it is. Default: `c("width", "abundance")`.

- ratio:

  Numeric value with the number of comparison rows drawn per region of
  the set, when matching. Default: `5`.

- strata:

  Numeric value with the number of strata used per covariate. Default:
  `5`.

- regionSets:

  Character vector with the names of the sets to build a universe for.
  Default: `NULL`, all of them.

- index:

  List with one vector of row positions per set, holding the comparison
  rows. Only for `type = "supplied"`. Default: `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.universe` object.

## Details

Whether a set responds "more than the rest" depends entirely on what the
rest is, so the choice is worth making explicitly rather than inheriting
it from whatever happened to be loaded. The comparison rows are drawn
from the other sets in the object, which asks whether a set responds
differently from the others it was loaded alongside. That is usually the
interesting claim when the sets were chosen to be compared with each
other, and it is the only claim available when the object holds nothing
else.

Matching on width and abundance is what keeps the answer from being
about the intervals rather than the biology. A set of 40 kb domains has
more reads per region than a set of 400 bp promoter windows, and more
reads mean a tighter fold change estimate, so an unmatched competitive
test can separate the two sets on precision alone. The rows are binned
on the quantiles of each covariate and the comparison is drawn within
the bins, at `ratio` rows per region of the set. When a stratum runs out
of eligible rows the whole stratum is taken and the shortfall shows up
in the diagnostics, where the medians of the set and of its comparison
should sit close together.

The abundance used for the matching is the width-adjusted one, the same
quantity
[`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
thresholds on, so a universe matched on abundance is also matched on
signal density rather than on total signal.

The draw inside a stratum is deterministic. The candidates are ordered
on the covariate being matched and taken evenly spaced across that
order, which spreads the comparison over the stratum, gives the same
universe on every call, and leaves the random number generator of the
session alone.

A single region set leaves nothing to compare against and this function
stops. Either load the sets that make the comparison interesting, or add
the genome bins as a set of their own before counting.

## See also

[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`plotUniverseMatching`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md),
and
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
which is a different thing entirely: it counts genome bins for the
normalisation and has nothing to do with the comparison built here.

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
