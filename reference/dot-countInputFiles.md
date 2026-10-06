# .countInputFiles

Counts every distinct input file once over the regions, with the read
filters of the samples, and spreads the counts over the samples they
serve.

## Usage

``` r
.countInputFiles(
  inputPaths,
  inputIds,
  sampleFragmentLength,
  fragmentLength,
  bamTargets,
  regions,
  maxFragmentLength,
  minMapq,
  removeDuplicates,
  excludeChromosomes,
  discardRegions,
  fullLibrarySize,
  nThreads,
  progressBar = FALSE,
  verbose
)
```

## Arguments

- inputPaths:

  Character vector with the input of every sample, `NA` where there is
  none.

- inputIds:

  Character vector with the name of the input of every sample, or `NULL`
  to name them after the files.

- sampleFragmentLength:

  Numeric vector with the fragment length of every sample, `NA` for the
  paired-end ones.

- fragmentLength:

  Value given to `countReads`, used to tell a single length from lengths
  varying by sample.

- bamTargets:

  Named vector with the chromosome lengths of the sample BAM files,
  whose names the regions are written in.

- regions:

  `GRanges` with the regions, named after the chromosomes of the BAM
  files.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- excludeChromosomes:

  Character vector with the chromosomes left out of the library sizes.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, or `NULL`.

- fullLibrarySize:

  Logical value indicating whether the library sizes cover every
  chromosome.

- nThreads:

  Number of threads.

- progressBar:

  Logical value to indicate whether a progress bar must be drawn while
  the inputs are read. Default: `FALSE`.

- verbose:

  Logical value to indicate whether the messages must be printed.

## Value

A list with `counts`, a matrix with one column per sample, `NA` for the
samples without input; `library.size`, the library size of the input of
every sample; and `input.id`, its name.

## Author

Sebastian Gregoricchio
