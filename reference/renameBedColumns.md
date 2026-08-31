# renameBedColumns

Assigns the standard BED column names to the first columns of a
data.frame, leaving the remaining ones untouched. Useful to convert a
headerless table, read with generic `V1`, `V2`, `V3` names, into a
data.frame ready for
[`GenomicRanges::makeGRangesFromDataFrame`](https://rdrr.io/pkg/GenomicRanges/man/makeGRangesFromDataFrame.html).

## Usage

``` r
renameBedColumns(table, bedFormat = 3)
```

## Arguments

- table:

  A data.frame whose first columns hold the genomic coordinates.

- bedFormat:

  Numeric value indicating how many columns must be renamed, one among
  `3` (seqnames, start, end), `6` (adding name, score, strand), `9`
  (adding thickStart, thickEnd, itemRgb) or `12` (adding blockCount,
  blockSizes, blockStarts). Default: `3`.

## Value

The input data.frame with the first columns renamed.

## Author

Sebastian Gregoricchio

## Examples

``` r
bedTable <- data.frame(V1 = "chr12",
                       V2 = c(1000, 5000),
                       V3 = c(2000, 6000),
                       V4 = c("peak_1", "peak_2"),
                       V5 = c(120, 340),
                       V6 = c("+", "-"))

renameBedColumns(bedTable, bedFormat = 6)
#>   seqnames start  end   name score strand
#> 1    chr12  1000 2000 peak_1   120      +
#> 2    chr12  5000 6000 peak_2   340      -

# Only the first three columns are renamed, the rest keep their names
renameBedColumns(bedTable, bedFormat = 3)
#>   seqnames start  end     V4  V5 V6
#> 1    chr12  1000 2000 peak_1 120  +
#> 2    chr12  5000 6000 peak_2 340  -
```
