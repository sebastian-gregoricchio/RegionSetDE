# .resolveFragmentLength

Works out the length every single-end sample is extended to, from a
single value, one value per sample, a column of the sample annotation,
or an estimate made on the reads.

## Usage

``` r
.resolveFragmentLength(
  fragmentLength,
  sampleTable,
  pairedEnd,
  bamFiles,
  regions,
  minMapq,
  removeDuplicates,
  discardRegions,
  nThreads,
  verbose
)
```

## Arguments

- fragmentLength:

  Value given to `countReads`.

- sampleTable:

  Data.frame with one row per sample, as built by `.buildSampleTable`.

- pairedEnd:

  Logical vector with one value per sample.

- bamFiles:

  Character vector with the paths of the BAM files.

- regions:

  `GRanges` with the regions, named after the chromosomes of the BAM
  files, used by the estimate.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, or `NULL`.

- nThreads:

  Number of threads.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with `values`, a numeric vector with one length per sample, `NA`
for the paired-end ones, and `source`, a string saying where the lengths
came from.

## Author

Sebastian Gregoricchio
