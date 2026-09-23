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
  directionFDR = 0.05,
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

- directionFDR:

  Numeric value with the false discovery rate, within each region, below
  which a tile counts as moving up or down in `n.tiles.up` and
  `n.tiles.down`. It is the `fc.threshold` of csaw, which despite its
  name is not a fold change. Default: `0.05`.

- adjustMethod:

  String with the multiple testing correction.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with the `results` data.frame and the `regions` `GRanges`.

## Author

Sebastian Gregoricchio
