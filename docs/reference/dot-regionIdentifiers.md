# .regionIdentifiers

Returns the identifier of every region, the `region.id` metadata column
when the regions carry one and their names otherwise.

Returns the identifier of every region of a set: the names of the
ranges, a metadata column when one is asked for, or the coordinates when
there is nothing else.

## Usage

``` r
.regionIdentifiers(gr, regionId = NULL, setName = "")

.regionIdentifiers(gr, regionId = NULL, setName = "")
```

## Arguments

- gr:

  `GRanges` of one region set.

- regionId:

  String with the name of a metadata column, or `NULL`.

- setName:

  String with the name of the set, used in the messages.

- regions:

  `GRanges` with the regions of the object.

## Value

A character vector with one identifier per region.

A character vector with one identifier per region.

## Author

Sebastian Gregoricchio
