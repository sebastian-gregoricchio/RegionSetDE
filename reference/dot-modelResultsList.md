# .modelResultsList

Collects the contrasts whose statistics the model brackets of
[`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
are written from: the single contrast of a result, every contrast of a
list, or the one named in `contrast`.

## Usage

``` r
.modelResultsList(object, contrast = NULL)
```

## Arguments

- object:

  Object handed to
  [`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md).

- contrast:

  String with the name of a contrast, or its position, or `NULL`.
  Default: `NULL`.

## Value

A named list of `RegionSetDE.results` objects.

## Author

Sebastian Gregoricchio
