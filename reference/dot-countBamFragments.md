# .countBamFragments

Counts the fragments of a group of BAM files over a set of ranges. The
chromosomes are cut into pieces of at most 50 Mb, and the pieces of all
the files are shared among the threads, so that even a single file keeps
every thread busy. Paired-end fragments are rebuilt from the first mate
of each proper pair, whose position, mate position and template length
(TLEN) give the start and the width of the fragment: the two reads never
have to be matched. Single-end reads are extended to the fragment length
from their 5' end.

## Usage

``` r
.countBamFragments(
  bamFiles,
  ranges,
  pairedEnd,
  fragmentLength = 150,
  maxFragmentLength = 1000,
  minMapq = 20,
  removeDuplicates = TRUE,
  restrictChromosomes = NULL,
  discardRegions = NULL,
  fullLibrarySize = TRUE,
  countMode = "overlap",
  pieceLength = 5e+07,
  nThreads = 1
)
```

## Arguments

- bamFiles:

  Character vector with the paths of the BAM files, all sharing the same
  header.

- ranges:

  `GRanges` with the ranges to count, named after the chromosomes of the
  BAM files. The strand is ignored.

- pairedEnd:

  Logical vector with one value per BAM file.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.
  Default: `150`.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.
  Default: `1000`.

- minMapq:

  Numeric value with the minimum mapping quality of a read. For
  paired-end data the mate is checked through its `MQ` tag, when the
  file carries it. Default: `20`.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded. Default: `TRUE`.

- restrictChromosomes:

  Character vector with the chromosomes to read, named as in the BAM
  files. Default: `NULL`, all of them.

- discardRegions:

  `GRanges` with the regions whose reads must be ignored, named after
  the chromosomes of the BAM files. A fragment is dropped when one of
  its reads starts inside them. Default: `NULL`.

- fullLibrarySize:

  Logical value: `TRUE` reads every chromosome to compute the library
  sizes, `FALSE` only the chromosomes carrying ranges, which is faster
  but leaves the library sizes partial. Default: `TRUE`.

- countMode:

  String with the way a fragment is assigned to the ranges: `"overlap"`
  counts it in every range it overlaps, `"bin"` counts it once, at its
  centre for paired-end data and at the 5' end of the read for
  single-end data, as csaw does for genome wide bins. Default:
  `"overlap"`.

- pieceLength:

  Numeric value with the maximum length of the stretch of genome read by
  a single job, in base pairs. Default: `5e7`.

- nThreads:

  Number of threads. Default: `1`.

## Value

A list with three elements: `counts`, an integer matrix with one row per
range and one column per file; `library.size`, the number of fragments
that went through the filters on the chromosomes read;
`mate.mapq.found`, telling for each paired-end file whether the `MQ` tag
was found (`NA` for single-end files).

## Author

Sebastian Gregoricchio
