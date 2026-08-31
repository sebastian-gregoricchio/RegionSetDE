# Coerce a counts object to a DGEList

Registers the `as(x, "DGEList")` idiom, which calls
[`asDGEList`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md)
with its defaults.

## Value

A `DGEList` holding the counts, the sample metadata, the region
annotation and the normalisation as offsets.

## See also

[`asDGEList`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (FALSE) { # \dontrun{
dgeList <- as(counts, "DGEList")
} # }
```
