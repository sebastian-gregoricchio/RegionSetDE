# computeProfiles

Computes the signal of every sample in bins around the centre of a group
of regions, by default the regions found to change in a contrast split
into those going up and those going down. The signal is read from the
BAM files the counts came from, with the same read filters and scaled by
the normalisation of the object, or from bigWig files. The result feeds
[`plotProfile`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotProfile.md).

## Usage

``` r
computeProfiles(
  object,
  regions = NULL,
  contrast = NULL,
  set = NULL,
  samples = NULL,
  direction = "both",
  FDR = NULL,
  log2FC = NULL,
  maxRegions = 1000,
  blacklist = NULL,
  whitelist = NULL,
  groupBy = NULL,
  signal = "auto",
  signalFiles = NULL,
  distance = 1500,
  binWidth = 50,
  centre = "summit",
  useOffsets = TRUE,
  nThreads = 1,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.results` or `RegionSetDE.resultsList` object carrying its
  counts, whose differential regions are profiled, or a
  `RegionSetDE.counts` or `RegionSetDE.fit` object, whose region sets
  are.

- regions:

  Regions to profile instead of those picked from `object`: a `GRanges`,
  a named list of `GRanges` or a `GRangesList`, one group of rows per
  element, or a `RegionSetDE` object, one group per set. Default:
  `NULL`.

- contrast:

  String with the name of a contrast, or its position, when `object`
  holds several of them. Default: `NULL`.

- set:

  Character vector with the names of the region sets used. Default:
  `NULL`, all of them.

- samples:

  Samples profiled: a character vector with their names, a numeric
  vector with their positions, or a logical vector with one value per
  sample. The scaling factors of the whole analysis are kept. Default:
  `NULL`, all of them.

- direction:

  String with the differential regions profiled, one among `"both"`,
  which draws the regions going up and those going down as two groups of
  rows, `"up"` and `"down"`. Only for results. Default: `"both"`.

- FDR:

  Numeric value with the adjusted p-value cut-off defining a
  differential region. Default: `NULL`, the one used by the test.

- log2FC:

  Numeric value with the absolute log2 fold change cut-off defining a
  differential region. Default: `NULL`, the one used by the test.

- maxRegions:

  Numeric value with the maximum number of regions per group of rows:
  the most significant ones for results, the strongest ones for the
  region sets of a counts object, and regions spread evenly along the
  list for regions given through `regions`. It is applied after
  `blacklist` and `whitelist`. Default: `1000`.

- blacklist:

  Regions whose signal must stay out of the figure: a row is left out
  when its drawn window, `distance` on each side of the centre, overlaps
  them. A `GRanges`, a path to a BED-like file, a data.frame, a list of
  them, or `TRUE` for the blacklist stored in the object by
  [`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
  or
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).
  Default: `NULL`.

- whitelist:

  Regions the rows are restricted to: a row is kept when its region
  overlaps them. Accepts the same forms as `blacklist`, except `TRUE`.
  Default: `NULL`.

- groupBy:

  String with the name of a `colData` column whose groups are averaged
  into one profile each, for instance `"condition"`. Default: `NULL`,
  one profile per sample.

- signal:

  String with the files the signal is read from: `"bam"`, the BAM files
  of
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
  `"bigwig"`, the bigWig files of the sample sheet or of
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md),
  or `"auto"`, the BAM files when they are known and the bigWig files
  otherwise. Default: `"auto"`.

- signalFiles:

  Character vector with one BAM or bigWig file per sample, in the order
  of the samples, overriding the ones recorded in the object. Default:
  `NULL`.

- distance:

  Numeric value with the number of base pairs drawn on each side of the
  centre. Default: `1500`.

- binWidth:

  Numeric value with the width of the bins the signal is averaged over,
  in base pairs. Default: `50`.

- centre:

  String with the point the regions are aligned on: `"summit"`, the
  `summit` column written by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  with `summits`, falling back to the midpoint where there is none, or
  `"midpoint"`. Default: `"summit"`.

- useOffsets:

  Logical value to indicate whether the signal read from BAM files must
  be scaled by the scaling factors of
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
  `FALSE` scaling it by the library sizes alone. Ignored for bigWig
  files, whose signal is taken as it is. Default: `TRUE`.

- nThreads:

  Number of threads, one file per thread. Default: `1`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with `profiles`, a named list with one matrix per sample, or per
group with `groupBy`, holding one row per region and one column per bin;
`regions`, a `GRanges` with the window of every row, its `row.group` and
its `region.key`; `bins`, the distance of the centre of every bin from
the centre of the regions; `samples`, the annotation of the samples; and
`parameters`.

## Details

The windows span `distance` base pairs on each side of the centre and
are cut into bins of `binWidth`, the signal of a bin being the mean
coverage over its bases. Regions on the minus strand are read from right
to left, so that a set of promoters shows the transcription start site
facing the same way. Windows running past the end of a chromosome are
left out, and how many were is reported.

From BAM files the coverage is built from the fragments as
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
counts them: paired-end fragments from their pairs, single-end reads
extended to the fragment length of their sample, with the same mapping
quality, duplicate and discarded region filters. It is then divided by
the scaling factor of the sample and expressed per million fragments of
the average library, so that the profiles of different samples compare
the way the normalised counts do. A normalisation with one offset per
region, as `method = "loess"` gives, has no factor to divide by, and the
library sizes are used instead. This is the counterpart of
`DiffBind::dba.plotProfile`, which also reads the reads of the analysis
and applies its normalisation.

A bigWig file is read as it is, so it has to be normalised already, and
to the same scale for every sample, for its profiles to be compared.

A region cleaned of blacklisted regions can still sit next to one, and
the window drawn around it then shows the artefact. `blacklist`
therefore acts on the whole window rather than on the region, while
`whitelist` asks the region itself to overlap the list, which is how a
profile is restricted to promoters, enhancers or any other class of
sites. How many rows each of them left out is reported.

## See also

[`plotProfile`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotProfile.md),
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)

counts <- countReads(peakRegions, sampleSheet = sampleSheet, summits = 200, verbose = FALSE)
counts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)

peakProfiles <- computeProfiles(counts, groupBy = "condition", distance = 1000, verbose = FALSE)
names(peakProfiles$profiles)
#> [1] "DMSO"      "R1881_4h"  "R1881_24h"
dim(peakProfiles$profiles$DMSO)
#> [1] 101  40
```
