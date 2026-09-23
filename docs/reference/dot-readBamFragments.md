# .readBamFragments

Reads the fragments of one BAM file over a group of pieces of genome and
counts them over the ranges of the same chromosomes. It is the job run
by every thread of `.countBamFragments`.

## Usage

``` r
.readBamFragments(
  job,
  bamFiles,
  pairedEnd,
  fragmentLength,
  maxFragmentLength,
  minMapq,
  removeDuplicates,
  countMode
)
```

## Arguments

- job:

  List describing the job: `file.index`, the `pieces` table, the
  `ranges` to count with their `range.chromosome` and `range.index`, and
  the `discard` regions of its chromosomes.

- bamFiles:

  Character vector with the paths of the BAM files.

- pairedEnd:

  Logical vector with one value per BAM file.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

- countMode:

  String, either `"overlap"` or `"bin"`, see `.countBamFragments`.

## Value

A list with the `file.index` and `range.index` of the job, the `counts`
of its ranges, the `total.fragments` that went through the filters and
`mate.mapq.found`.

## Author

Sebastian Gregoricchio
