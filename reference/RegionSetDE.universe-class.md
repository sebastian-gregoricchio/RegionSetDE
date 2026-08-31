# RegionSetDE.universe class

S4 class holding, for every region set, the rows that make up the
universe of its competitive test: the set itself together with the rows
it is compared against. The positions refer to the object the universe
was built from, so the same universe serves every contrast run on that
object.

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

- `diagnostics`:

  Data.frame with, for every set, the median width and abundance of the
  set and of the rows it is compared against.

- `n.rows`:

  Numeric value with the number of rows of the object the positions
  refer to.

## See also

[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)

## Author

Sebastian Gregoricchio
