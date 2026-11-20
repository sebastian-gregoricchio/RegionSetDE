# .checkCountBigwigArguments

Checks the arguments of
[`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkCountBigwigArguments(
  regionSet,
  bigwigFiles,
  sampleSheet,
  sampleNames,
  sampleMetadata,
  summaryFunction,
  countLike,
  roundValues,
  BPPARAM,
  verbose
)
```

## Arguments

- regionSet:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- bigwigFiles:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- sampleSheet:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- sampleNames:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- sampleMetadata:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- summaryFunction:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- countLike:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- roundValues:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- BPPARAM:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

- verbose:

  See
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md).

## Value

A list with the arguments as checked and completed: `bigwigFiles`,
`parallelParam`, `regionSet`, `sampleMetadata`, `sampleNames`,
`summaryFunction`.

## Author

Sebastian Gregoricchio
