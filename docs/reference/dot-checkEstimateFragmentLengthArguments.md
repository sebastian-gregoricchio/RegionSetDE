# .checkEstimateFragmentLengthArguments

Checks the arguments of
[`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkEstimateFragmentLengthArguments(
  bamFiles,
  sampleSheet,
  sampleNames,
  pairedEnd,
  maxDistance,
  discardRegions,
  maxReads,
  verbose
)
```

## Arguments

- bamFiles:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- sampleSheet:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- sampleNames:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- pairedEnd:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- maxDistance:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- discardRegions:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- maxReads:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

- verbose:

  See
  [`estimateFragmentLength`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md).

## Value

A list with the arguments as checked and completed: `bamFiles`,
`maxDistance`, `pairedEnd`, `sampleNames`.

## Author

Sebastian Gregoricchio
