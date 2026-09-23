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
`paired.end`, `bam.reads` (every record of the BAM file), `bam.mapped`
(the mapped ones), `library.size` (the fragments that went through the
filters of the counting), `reads.in.regions` and `FRiP`, the ratio of
the last two.

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

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`countTable`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countTable.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (FALSE) { # \dontrun{
sampleSheet <- loadExampleData("peakSheet")
consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl")
counts <- countReads(consensus, sampleSheet = sampleSheet)

libInfo(counts, annotationColumns = "condition")
} # }
```
