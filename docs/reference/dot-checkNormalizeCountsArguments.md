# .checkNormalizeCountsArguments

Checks the arguments of
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkNormalizeCountsArguments(
  counts,
  method,
  scalingFactors,
  factorType,
  spikeInCounts,
  greenlistEstimator,
  backgroundHoldout
)
```

## Arguments

- counts:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

- method:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

- scalingFactors:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

- factorType:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

- spikeInCounts:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

- greenlistEstimator:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

- backgroundHoldout:

  See
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md).

## Value

A list with the arguments as checked and completed: `backgroundHoldout`,
`countMatrix`, `factorType`, `greenlistEstimator`, `method`,
`sampleNumber`.

## Author

Sebastian Gregoricchio
