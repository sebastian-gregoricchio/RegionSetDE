# .readSummits

Places the summit of every region from the fragment pileup of the
samples: the summit of each sample, averaged with weights proportional
to the height of its pileup over its depth.

## Usage

``` r
.readSummits(
  bamFiles,
  regions,
  pairedEnd,
  fragmentLength,
  maxFragmentLength,
  minMapq,
  removeDuplicates,
  discardRegions,
  nThreads
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

  Numeric vector with the fragment length of every sample, `NA` for the
  paired-end ones.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, or `NULL`.

- nThreads:

  Number of threads.

## Value

An integer vector with the summit of every region, `NA` where no sample
has a fragment.

## Author

Sebastian Gregoricchio
