# .bamSummits

Finds, in every BAM file, the summit of the fragment pileup within each
region: the middle of the highest stretch of coverage, and its height.

## Usage

``` r
.bamSummits(
  bamFiles,
  regions,
  pairedEnd,
  fragmentLength,
  maxFragmentLength,
  minMapq,
  removeDuplicates,
  discardRegions = NULL,
  progressBar = FALSE,
  nThreads = 1
)
```

## Arguments

- bamFiles:

  Character vector with the paths of the BAM files.

- regions:

  `GRanges` with the regions, named after the chromosomes of the BAM
  files.

- pairedEnd:

  Logical vector with one value per BAM file.

- fragmentLength:

  Numeric vector with the length to which single-end reads are extended,
  one value per BAM file.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, or `NULL`.

- progressBar:

  Logical value to indicate whether a progress bar must be drawn, one
  step for every file. Default: `FALSE`.

- nThreads:

  Number of threads, one file per thread. Default: `1`.

## Value

A list with `position` and `height`, two matrices with one row per
region and one column per file, and `fragments`, the number of fragments
read in each file.

## Author

Sebastian Gregoricchio
