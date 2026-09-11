# .collapseTileStats

Averages the tiles of a region back into a single row, so that the
region rather than the tile is the unit a set is built from.

## Usage

``` r
.collapseTileStats(regionStats)
```

## Arguments

- regionStats:

  Data.frame returned by `.setStatistics`, one row per tile.

## Value

A list with the collapsed `stats` table and `map`, an integer vector
giving the collapsed row of every tile.

## Details

A set assembled from tiles gives every region a weight equal to the
number of tiles it was cut into, which turns the mean over the set into
a mean over base pairs. Averaging first restores the region as the unit.
The averaged statistic is deliberately the plain mean of the per-tile
statistics rather than a combined one: it keeps the scale of the values
comparable between a region cut into two tiles and one cut into forty,
which is what the competitive test ranks them on.

## Author

Sebastian Gregoricchio
