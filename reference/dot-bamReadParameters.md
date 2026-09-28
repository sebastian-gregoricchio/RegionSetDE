# .bamReadParameters

Builds the `ScanBamParam` shared by every function reading fragments
from a BAM file, so that counting, summits and profiles apply the same
filters. Of a proper pair only the first mate is read, since its
position, the position of its mate and the template length already
describe the fragment.

## Usage

``` r
.bamReadParameters(which, isPairedEnd, minMapq, removeDuplicates)
```

## Arguments

- which:

  `GRanges` with the stretches of genome to read.

- isPairedEnd:

  Logical value, `TRUE` for a paired-end file.

- minMapq:

  Numeric value with the minimum mapping quality of a read.

- removeDuplicates:

  Logical value indicating whether the reads flagged as duplicates must
  be discarded.

## Value

A `ScanBamParam` object.

## Author

Sebastian Gregoricchio
