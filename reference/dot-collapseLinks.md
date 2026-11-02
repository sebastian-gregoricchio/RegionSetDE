# .collapseLinks

Picks one link per region, for the columns written next to the regions.
The closest transcription start site when the distances are known, the
heaviest link when the weights are, and the first one otherwise. It is a
view for reading and for labelling a plot, never the quantity anything
is computed from: that is the full table, which
[`regionLinks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md)
returns.

## Usage

``` r
.collapseLinks(linkTable, regionKeys)
```

## Arguments

- linkTable:

  Data.frame with the validated links.

- regionKeys:

  Character vector with the key of every region of the object.

## Value

A data.frame with one row per region, in the order of `regionKeys`.

## Author

Sebastian Gregoricchio
