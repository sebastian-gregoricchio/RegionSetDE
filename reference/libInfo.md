# libInfo

Summarises the libraries of a counts object in one table: how many reads
each BAM file holds, how many fragments went through the read filters,
how many of those fall in the regions, and the fraction of reads in
regions (FRiP) that follows from the two.

## Usage

``` r
libInfo(counts, annotationColumns = NULL, bamFiles = NULL)
```

## Arguments

- counts:

  `RegionSetDE.counts` object returned by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
  or a `RegionSetDE.fit`, whose counts are used.

- annotationColumns:

  Character vector with the columns of the `colData` to add after the
  sample names, for instance `"condition"`. Default: `NULL`, none.

- bamFiles:

  Character vector with the BAM files, in the order of the samples.
  Default: `NULL`, the files recorded by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md).

## Value

A data.frame with one row per sample: `sample`, the `annotationColumns`,
`paired.end`, `fragment.length` (the length single-end reads were
extended to, `NA` for paired-end samples), `bam.reads` (every record of
the BAM file), `bam.mapped` (the mapped ones), `library.size` (the
fragments that went through the filters of the counting),
`reads.in.regions` and `FRiP`, the ratio of the last two. When
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
discarded the reads of a blacklist, of a greylist or of other regions,
`discarded.reads` follows `library.size` with the fragments each sample
lost to them. When the inputs were counted, `input.id`,
`input.library.size`, `input.in.regions` and `input.FRiP` give the same
numbers for the input of each sample.

## Details

The three counts measure different things and are not expected to agree.
`bam.reads` and `bam.mapped` come from the index of each file, so they
are read in an instant, and they count alignment records: a paired-end
fragment is two of them. `library.size` and `reads.in.regions` come from
the counting, where a paired-end fragment counts once and only after the
mapping quality, duplicate and proper pair filters. On paired-end data
`bam.mapped` is therefore about twice `library.size`, less what the
filters removed.

`reads.in.regions` counts a fragment once per region it overlaps, the
way the counting does. A region shared by several sets is counted once,
but a fragment lying across two neighbouring regions, or two tiles of
the same region, is counted in both, so on a tiled object the FRiP comes
out slightly high.

The FRiP is computed on `library.size`, the reads that went through the
same filters as the counts. Computed on `bam.reads` it would mix
fragments with alignment records and mapped with filtered reads.

`discarded.reads` counts the fragments that passed every other filter
and lie on the discarded regions, on the chromosomes entering the
library sizes. They are in neither `library.size` nor
`reads.in.regions`, and a sample losing a much larger share of its
library than the others is worth a look.

The FRiP of the input is the share of the input library falling in the
regions, which is what the regions would collect with no enrichment at
all. A sample whose FRiP sits close to the one of its input carries
little signal in the regions, however many reads it has.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`countTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countTable.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)

libInfo(counts, annotationColumns = "condition")
#>            sample condition paired.end fragment.length bam.reads bam.mapped
#> 1      AR_DMSO_r1      DMSO       TRUE              NA      7913       7913
#> 2      AR_DMSO_r2      DMSO       TRUE              NA     12038      12038
#> 3      AR_DMSO_r3      DMSO       TRUE              NA      9185       9185
#> 4  AR_R1881_4h_r1  R1881_4h       TRUE              NA      7836       7836
#> 5  AR_R1881_4h_r2  R1881_4h       TRUE              NA     11205      11205
#> 6  AR_R1881_4h_r3  R1881_4h       TRUE              NA     12693      12693
#> 7 AR_R1881_24h_r1 R1881_24h       TRUE              NA     14080      14080
#> 8 AR_R1881_24h_r2 R1881_24h       TRUE              NA     10793      10793
#> 9 AR_R1881_24h_r3 R1881_24h       TRUE              NA     10332      10332
#>   library.size reads.in.regions   FRiP input.id input.library.size
#> 1         3588              116 0.0323    input               5353
#> 2         5469              161 0.0294    input               5353
#> 3         4183              131 0.0313    input               5353
#> 4         3582              558 0.1558    input               5353
#> 5         5121              881 0.1720    input               5353
#> 6         5806              718 0.1237    input               5353
#> 7         6326             1139 0.1801    input               5353
#> 8         4889             1158 0.2369    input               5353
#> 9         4651              839 0.1804    input               5353
#>   input.in.regions input.FRiP
#> 1               26     0.0049
#> 2               26     0.0049
#> 3               26     0.0049
#> 4               26     0.0049
#> 5               26     0.0049
#> 6               26     0.0049
#> 7               26     0.0049
#> 8               26     0.0049
#> 9               26     0.0049
```
