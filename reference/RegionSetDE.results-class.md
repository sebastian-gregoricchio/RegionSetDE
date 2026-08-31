# RegionSetDE.results class

S4 class holding the outcome of a per-region contrast. The table in the
`results` slot has one row per region and follows the coordinates in the
`regions` slot, so the two can be bound together without any matching
step. When the counts were tiled, the tile level statistics are kept in
the `tiles` slot next to the combined ones, since a region that changes
only over part of its length is easier to read from the tiles than from
the summary.

## Slots

- `results`:

  Data.frame with one row per region.

- `tiles`:

  Data.frame with one row per tile, empty when the counts were not
  tiled.

- `regions`:

  `GRanges` with the coordinates of the rows of `results`.

- `contrast`:

  String describing the contrast that was tested.

- `contrast.vector`:

  Numeric vector with the coefficients of the contrast over the columns
  of the design.

- `engine`:

  String with the engine that produced the statistics.

- `counting.level`:

  String indicating whether the model was fitted on regions or on tiles.

- `combination`:

  List with the method used to combine the tiles and whether it was
  applied.

- `thresholds`:

  List with the FDR and log2 fold change cut-offs used to fill the
  `diff.status` column.

- `counts`:

  `RegionSetDE.counts` object the contrast was computed on, so that the
  values behind a result can be drawn without carrying a second object
  around. Empty when `testRegions` was called with
  `carryCounts = FALSE`.

## Author

Sebastian Gregoricchio
