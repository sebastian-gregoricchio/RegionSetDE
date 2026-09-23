# countGreenlist

Counts the reads falling in the regions of a CUT&RUN or CUT&Tag
greenlist, the places where the background of the protocol is
reproducible enough to measure how much material was sequenced. The
counts are stored in the metadata of the object, where
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
picks them up with `method = "greenlist"`.

## Usage

``` r
countGreenlist(
  counts,
  greenlist,
  bamFiles = NULL,
  minCount = 1,
  pairedEnd = NULL,
  fragmentLength = NULL,
  maxFragmentLength = NULL,
  minMapq = NULL,
  removeDuplicates = NULL,
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object returned by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

- greenlist:

  `GRanges` with the greenlist regions, typically from
  [`loadGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md),
  or the path to a BED file holding them.

- bamFiles:

  Character vector with the paths of the BAM files, in the same order as
  the samples of `counts`. Default: `NULL`, the files recorded by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  are reused.

- minCount:

  Numeric value with the minimum total count required to keep a
  greenlist region. Default: `1`.

- pairedEnd:

  Logical value, or one logical value per BAM file, indicating whether
  the reads must be counted as proper pairs. Default: `NULL`, the
  layouts resolved at the counting step.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.
  Default: `NULL`, the value used at the counting step.

- maxFragmentLength:

  Numeric value with the maximum insert size accepted for a pair.
  Default: `NULL`, the value used at the counting step.

- minMapq:

  Numeric value with the minimum mapping quality of a read. Default:
  `NULL`, the value used at the counting step.

- removeDuplicates:

  Logical value indicating whether the duplicated reads must be
  discarded. Default: `NULL`, the value used at the counting step.

- nThreads:

  Number of threads. The files are cut into pieces of at most 50 Mb,
  shared among the threads. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The input `RegionSetDE.counts` object with the greenlist counts stored
as a `RangedSummarizedExperiment` in `metadata(counts)$greenlist`.

## Details

The read filters are taken from
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
unless they are given here, for the same reason as in
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md):
a reference counted with another mapping quality or duplicate policy
describes a library that is not the one under study.

Each fragment is counted once, in the region holding its centre, as the
background bins do, so the totals stay a share of the library and two
neighbouring regions never claim the same fragment. The list is merged
beforehand for the same reason. Greenlist regions lying on chromosomes
absent from the BAM files are dropped, and how many were is reported,
which is what catches a list built for another assembly before it
quietly halves the counts.

## See also

[`loadGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (FALSE) { # \dontrun{
greenlist <- loadGreenlist("hg38", assay = "cutrun")

counts <- countGreenlist(counts, greenlist = greenlist, nThreads = 4)
counts <- normalizeCounts(counts, method = "greenlist")

SummarizedExperiment::colData(counts)$scaling.factor
} # }
```
