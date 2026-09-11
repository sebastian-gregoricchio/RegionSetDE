# .combineTiles

Combines the tile level statistics into one row per region, through
[`csaw::combineTests`](https://rdrr.io/pkg/csaw/man/combineTests.html).

## Usage

``` r
.combineTiles(
  tileTable,
  tileRanges,
  extraColumns = character(0),
  method = "simes",
  lfcThreshold = 0,
  adjustMethod = "BH",
  verbose = TRUE
)
```

## Arguments

- tileTable:

  Data.frame with one row per tile, as returned by the engine specific
  test.

- tileRanges:

  `GRanges` with the coordinates of the tiles.

- extraColumns:

  Character vector with the annotation columns carried over from the
  tiles.

- method:

  String with the combination method.

- lfcThreshold:

  Numeric value used to count the tiles moving in each direction.

- adjustMethod:

  String with the multiple testing correction.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with the `results` data.frame and the `regions` `GRanges`.

## Author

Sebastian Gregoricchio
