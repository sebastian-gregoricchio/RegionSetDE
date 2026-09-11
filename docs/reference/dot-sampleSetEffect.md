# .sampleSetEffect

Computes one set score per biological sample and runs it through the
design, so that the interval on the effect rests on the samples rather
than on the genomic loci.

## Usage

``` r
.sampleSetEffect(
  expressionMatrix,
  setIndex,
  backgroundIndex,
  design,
  contrastVector,
  level = 0.95
)
```

## Arguments

- expressionMatrix:

  Numeric matrix of log2 values, one row per region and one column per
  sample.

- setIndex:

  Integer vector with the rows of the set.

- backgroundIndex:

  Integer vector with the rows the set is compared against.

- design:

  Design matrix of the fit.

- contrastVector:

  Numeric vector with the contrast, in the columns of the design.

- level:

  Numeric value with the confidence level. Default: `0.95`.

## Value

A list with the per-sample scores, the estimate of the contrast on them,
its standard error, degrees of freedom, p-value and the bounds of the
interval.

## Details

The score of a sample is the mean signal over the set minus the mean
signal over its comparison, inside that library. Taking the difference
within the sample removes anything that scales the whole library, the
sequencing depth and the normalisation factor included, which is what
makes the score comparable across samples in the first place. Those
scores are then a single response fitted on the design of the
experiment, and the interval that comes out has as many degrees of
freedom as the design leaves, whether the set holds twenty regions or
thirty thousand.

This is a different quantity from the interval built on the per-region
fold changes, and it is the one that answers the question a reader takes
a set-level confidence interval to be answering. The regions still
contribute, through the precision of each sample's score, but they are
not counted as replicates.

## Author

Sebastian Gregoricchio
