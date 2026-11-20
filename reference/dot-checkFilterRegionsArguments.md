# .checkFilterRegionsArguments

Checks the arguments of
[`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkFilterRegionsArguments(counts, method, proportion, referenceWidth, assay)
```

## Arguments

- counts:

  See
  [`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md).

- method:

  See
  [`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md).

- proportion:

  See
  [`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md).

- referenceWidth:

  See
  [`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md).

- assay:

  See
  [`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md).

## Value

A list with the arguments as checked and completed: `countMatrix`,
`librarySizes`, `referenceWidth`, `rowTable`, `rowWidths`.

## Author

Sebastian Gregoricchio
