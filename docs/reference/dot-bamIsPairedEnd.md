# .bamIsPairedEnd

Tells whether BAM files hold paired-end reads, from the flags of their
first records. In a paired-end library every read carries the paired
flag, so the head of the file answers for all of it.
[`Rsamtools::testPairedEndBam`](https://rdrr.io/pkg/Rsamtools/man/testPairedEndBam.html)
stops at the first paired read as well, but on a single-end file it
finds none and reads to the last record, a million at a time, printing
the running total.

## Usage

``` r
.bamIsPairedEnd(bamFiles, nRecords = 1e+05)
```

## Arguments

- bamFiles:

  Character vector with the paths of the BAM files.

- nRecords:

  Numeric value with the number of records read from the head of each
  file. Default: `1e5`.

## Value

A logical vector with one value per file, `FALSE` for a file without
records.

## Author

Sebastian Gregoricchio
