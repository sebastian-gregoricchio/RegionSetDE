# .fragmentsFromReads

Turns the records returned by `scanBam` into fragments. Paired-end
fragments are rebuilt from the first mate of each proper pair,
single-end reads are extended from their 5' end to the fragment length.
The records of a pair longer than the maximum, of a mate below the
mapping quality, or starting in a discarded region are dropped.

## Usage

``` r
.fragmentsFromReads(
  reads,
  keep,
  isPairedEnd,
  fragmentLength,
  maxFragmentLength,
  minMapq,
  chromosomeLength,
  discard = NULL
)
```

## Arguments

- reads:

  List returned by `scanBam` for one stretch of genome, or several of
  them pasted together.

- keep:

  Logical vector with one value per record, telling which records enter
  at all.

- isPairedEnd:

  Logical value, `TRUE` for a paired-end file.

- fragmentLength:

  Numeric value with the length to which single-end reads are extended.

- maxFragmentLength:

  Numeric value with the maximum length of a paired-end fragment.

- minMapq:

  Numeric value with the minimum mapping quality of a read, applied to
  the mate through its `MQ` tag.

- chromosomeLength:

  Numeric value, or one value per record, with the length of the
  chromosome the fragments are clipped to.

- discard:

  `IRanges` with the discarded regions of the chromosome, or `NULL`.
  Only for records of a single chromosome.

## Value

A list with the `start`, the `end` and the `point` of every fragment,
the point being its centre for paired-end data and the 5' end of the
read for single-end data, the `index` of the record each fragment comes
from, and `mate.mapq.found`.

## Author

Sebastian Gregoricchio
