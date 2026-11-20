# .checkMakeGreylistArguments

Checks the arguments of
[`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkMakeGreylistArguments(
  inputFiles,
  inputNames,
  binSize,
  quantile,
  maxGap,
  minInputs,
  excludeChromosomes,
  pairedEnd
)
```

## Arguments

- inputFiles:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- inputNames:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- binSize:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- quantile:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- maxGap:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- minInputs:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- excludeChromosomes:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

- pairedEnd:

  See
  [`makeGreylist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md).

## Value

A list with the arguments as checked and completed: `binSize`,
`inputFiles`, `inputNames`, `maxGap`, `minInputs`, `pairedEnd`.

## Author

Sebastian Gregoricchio
