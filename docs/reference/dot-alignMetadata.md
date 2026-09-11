# .alignMetadata

Builds the metadata table of one region set in the shape every set has
to share, filling with `NA` the columns that set does not carry.

## Usage

``` r
.alignMetadata(gr, plan)
```

## Arguments

- gr:

  `GRanges` of one region set.

- plan:

  List returned by `.metadataPlan`.

## Value

A `DataFrame` with one row per region and one column per entry of the
plan.

## Author

Sebastian Gregoricchio
