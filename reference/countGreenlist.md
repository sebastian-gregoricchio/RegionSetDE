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
  excludeCounted = TRUE,
  bamFiles = NULL,
  minCount = 1,
  pairedEnd = NULL,
  fragmentLength = NULL,
  maxFragmentLength = NULL,
  minMapq = NULL,
  removeDuplicates = NULL,
  nThreads = 1,
  progressBar = interactive(),
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

- excludeCounted:

  Logical value to indicate whether the greenlist regions overlapping
  the regions of `counts` must be left out, since the reads there carry
  the signal under study rather than the background of the protocol.
  Default: `TRUE`.

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

  Numeric value with the length to which single-end reads are extended,
  or one value per BAM file. Default: `NULL`, the lengths used at the
  counting step, sample by sample.

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

- progressBar:

  Logical value to indicate whether a progress bar must be drawn while
  the files are read. It advances with the pieces of the files as the
  threads hand them back, and it is drawn only when `verbose = TRUE`.
  Default: [`interactive()`](https://rdrr.io/r/base/interactive.html),
  which keeps it out of scripts and rendered documents.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The input `RegionSetDE.counts` object with the greenlist counts stored
as a `RangedSummarizedExperiment` in `metadata(counts)$greenlist`. Its
`colData` describes every library over the list: `totals`, the fragments
counted, `regions.covered`, the regions holding at least one of them,
and `library.fraction`, the share of the library they represent.

## Details

The read filters are taken from
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
unless they are given here, for the same reason as in
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md):
a reference counted with another mapping quality or duplicate policy
describes a library that is not the one under study. The reads
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
discarded, those of the blacklist, of the greylist and of
`discardRegions`, are left out here too.

Each fragment is counted once, in the region holding its centre, as the
background bins do, so the totals stay a share of the library and two
neighbouring regions never claim the same fragment. The list is merged
beforehand for the same reason. This is the quantification of the
greenlist paper, which counted the lists with
`multiBamSummary --centerReads`. Greenlist regions lying on chromosomes
absent from the BAM files are dropped, and how many were is reported,
which is what catches a list built for another assembly before it
quietly halves the counts.

The checks follow the conditions under which de Mello *et al.* built and
validated the lists. The regions were chosen at least 5 kb away from
genes so that no target binds there; a region of this experiment that
overlaps one of them means the target does bind there, in these cells,
and its reads would carry the biology into the factors, which is why
`excludeCounted` drops it. The libraries used to build the lists had at
least 1.5 million aligned reads for human CUT&RUN, 1 million for mouse
CUT&RUN and 500,000 for human CUT&Tag, and a library below that depth,
counted here in fragments, is reported. So is a library covering fewer
than half as many greenlist regions as the median library, whose factor
rests on a small part of the list. The paper sets no rule per group of
samples: what matters for the factors is how many regions carry reads in
every library, which
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
reports when it computes them.

## See also

[`loadGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)

# A published list, loadGreenlist("hg38", assay = "cuttag"), holds too few reads in these small
# slices of chromosome 19, so a few stretches of the window stand in for it here
greenlist <- GenomicRanges::GRanges("19", IRanges::IRanges(start = seq(46.5e6, 57.5e6, by = 1e6), width = 2e5))

counts <- countGreenlist(counts, greenlist = greenlist, verbose = FALSE)
counts <- normalizeCounts(counts, method = "greenlist", verbose = FALSE)
#> Warning: Only 4 greenlist regions carry a read in every sample, the factors rest on those alone.

SummarizedExperiment::colData(counts)$scaling.factor
#> [1] 0.8202319 1.2421577 0.8396263 0.8382576 0.9573529 1.2259310 1.1672375
#> [8] 0.9747252 0.9344800
```
