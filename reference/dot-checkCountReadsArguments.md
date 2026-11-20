# .checkCountReadsArguments

Checks the arguments of
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
and brings them to the form the rest of the function works with.

## Usage

``` r
.checkCountReadsArguments(
  regionSet,
  bamFiles,
  sampleSheet,
  sampleNames,
  sampleMetadata,
  tileWidth,
  pairedEnd,
  fragmentLength,
  maxFragmentLength,
  excludeChromosomes,
  greylist,
  discardListedReads,
  discardRegions,
  fullLibrarySize,
  countInput,
  summits,
  summitSource,
  verbose
)
```

## Arguments

- regionSet:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- bamFiles:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- sampleSheet:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- sampleNames:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- sampleMetadata:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- tileWidth:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- pairedEnd:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- fragmentLength:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- maxFragmentLength:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- excludeChromosomes:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- greylist:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- discardListedReads:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- discardRegions:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- fullLibrarySize:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- countInput:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- summits:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- summitSource:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- verbose:

  See
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

## Value

A list with the arguments as checked and completed: `bamFiles`,
`bamSeqlevels`, `bamTargets`, `excludeChromosomes`, `pairedEnd`,
`regionSet`, `sampleMetadata`, `sampleNames`, `summitSource`.

## Author

Sebastian Gregoricchio
