# .representativeRow

Finds, for every region of a result table, the row of the counts that
carries its statistic. On a tiled object that is the tile the
combination picked, and on a non-tiled one it is the region itself.

## Usage

``` r
.representativeRow(rowTable, topTable, results)
```

## Arguments

- rowTable:

  Data.frame with the rows of the counts and their positions.

- topTable:

  Data.frame with the regions to draw.

- results:

  `RegionSetDE.results` object.

## Value

An integer vector with one position per region, `NA` when the region is
absent from the counts.

## Author

Sebastian Gregoricchio
