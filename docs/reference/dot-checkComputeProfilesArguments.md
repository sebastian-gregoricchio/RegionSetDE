# .checkComputeProfilesArguments

Checks the arguments of
[`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkComputeProfilesArguments(
  object,
  contrast,
  samples,
  direction,
  maxRegions,
  groupBy,
  signalFiles,
  distance,
  binWidth,
  centre
)
```

## Arguments

- object:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- contrast:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- samples:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- direction:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- maxRegions:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- groupBy:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- signalFiles:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- distance:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- binWidth:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

- centre:

  See
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md).

## Value

A list with the arguments as checked and completed: `allCounts`,
`binWidth`, `binsPerSide`, `centre`, `counts`, `direction`, `distance`,
`results`, `sampleIndex`, `sampleTable`, `signalFiles`.

## Author

Sebastian Gregoricchio
