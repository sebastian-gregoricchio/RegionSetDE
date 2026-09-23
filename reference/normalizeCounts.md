# normalizeCounts

Estimates the scaling factors of a `RegionSetDE.counts` object and
stores them together with a normalised assay. The factors can be
computed from the counts themselves, from the background bins collected
by
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
from the greenlist regions collected by
[`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md),
from a spike-in, or supplied by the user.

## Usage

``` r
normalizeCounts(
  counts,
  method = "TMM",
  scalingFactors = NULL,
  factorType = "division",
  spikeInCounts = NULL,
  greenlistCounts = NULL,
  greenlistEstimator = "medianRatio",
  useRegionSets = NULL,
  minCount = 1,
  referenceSample = NULL,
  backgroundHoldout = 0,
  normalizedAssay = "norm.counts",
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object returned by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
  or
  [`loadCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md).

- method:

  String with the method used to estimate the factors, one among
  `"TMM"`, `"TMMwsp"`, `"RLE"`, `"upperQuartile"`, `"librarySize"`,
  `"readsInRegions"`, `"background"`, `"greenlist"`, `"loess"`,
  `"spikeIn"`, `"manual"` or `"none"`. Default: `"TMM"`.

- scalingFactors:

  Numeric vector with one factor per sample, required by the `"manual"`
  method. Named vectors are matched to the sample names and may cover
  samples absent from the object, unnamed ones must follow the column
  order. Default: `NULL`.

- factorType:

  String declaring how the values of `scalingFactors` must be applied,
  either `"division"` when the counts have to be divided by them, as for
  the size factors of DESeq2, or `"multiplication"` when they have to be
  multiplied, as for the scale factors of deeptools and of most spike-in
  protocols. Default: `"division"`.

- spikeInCounts:

  Numeric vector with the number of reads assigned to the exogenous
  genome in each sample, required by the `"spikeIn"` method. Default:
  `NULL`.

- greenlistCounts:

  Counts over the greenlist regions when they were collected outside the
  package, either a numeric vector with one total per sample or a matrix
  with one row per region and one column per sample. Default: `NULL`,
  the counts stored by
  [`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md).

- greenlistEstimator:

  String with the way the greenlist counts become factors, one among
  `"medianRatio"`, the median of the ratios to the geometric mean of
  each region, which is the size factor of DESeq2 and what the greenlist
  paper used, `"TMM"`, and `"sum"`, the total signal over the list.
  Default: `"medianRatio"`.

- useRegionSets:

  Character vector with the names of the region sets used to estimate
  the factors. Default: `NULL`, all of them.

- minCount:

  Numeric value with the minimum total count required for a region to
  take part in the estimation. Default: `1`.

- referenceSample:

  String or numeric position of the sample used as reference by the
  `"TMM"` and `"TMMwsp"` methods. Default: `NULL`, chosen by edgeR.

- backgroundHoldout:

  Numeric value between 0 and 1 with the fraction of the background bins
  kept out of the estimation of the factors, so that
  [`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
  and
  [`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
  can work on bins the normalisation has never seen. Only for
  `method = "background"`. Default: `0`, every bin is used.

- normalizedAssay:

  String with the name given to the normalised assay. Default:
  `"norm.counts"`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

The input `RegionSetDE.counts` object with three additions: the
`norm.factor` and `scaling.factor` columns in the `colData`, the
normalised values in the assay named after `normalizedAssay`, and, for
the `"loess"` method, the log offsets in the `offset` assay. With
`backgroundHoldout` above zero the positions of the held-out bins are
stored in the `background.holdout` entry of the metadata.

## Details

The raw counts are never overwritten, the normalisation only adds an
assay and the factors beside it, so that the testing functions can keep
working on the counts and their offsets.

Whatever the method, `scaling.factor` always holds a divisor: the
normalised values are the raw counts divided by it, and a sample
sequenced more deeply than the others therefore gets a factor above one.
This is the convention of the size factors of DESeq2 and the opposite of
the scale factors written by `bamCoverage` or by the spike-in pipelines,
which are meant to multiply the signal. Factors coming from those tools
must be declared with `factorType = "multiplication"` and are inverted
on the way in. The factors are centred so that their mean is one, which
keeps the normalised values on the scale of the raw counts instead of
collapsing them to fractions.

A calibration check is only worth as much as the independence behind it.
The background bins serve twice, once here to estimate the factors and
once in
[`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
to estimate the dispersion, and a bin used in both has contributed to
the model it is later asked to test. Splitting the bins here, through
`backgroundHoldout`, keeps a fraction of them out of the factors, and
[`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
then picks up exactly that split rather than making its own. The
held-out bins are outside the whole preprocessing chain, which is the
version of the check that means what it appears to mean. Without it the
check still works, but on bins that are holdouts of the dispersion
alone.

The choice of the method matters more than usual on region sets. `"TMM"`
and `"RLE"` assume that most of the regions do not change, which is
reasonable for a catalogue of thousands of peaks but not for a handful
of hand-picked ones, and not for a mark that is globally redistributed
by the treatment. `"background"` sidesteps that assumption by estimating
the factors on the genome wide bins, where the signal of the experiment
is diluted, and is the safest option when a global shift is expected.
`"spikeIn"` relies on the exogenous genome alone and ignores the regions
altogether. `"loess"` corrects a bias that changes with the abundance,
which no single factor per sample can describe, so it returns a matrix
of offsets rather than a vector and leaves `scaling.factor` empty.

`"greenlist"` is the spike-in free reference of CUT&RUN and CUT&Tag. The
greenlist regions carry the background of the protocol rather than the
factor being mapped, reproducibly enough between experiments that the
reads landing there measure how much material was sequenced.
[`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md)
counts the libraries over them and this method turns those counts into
factors, by default with the median of ratios, which is the size factor
the greenlist paper computed with DESeq2. Unlike the endogenous methods
it says nothing about the regions under study, and unlike a spike-in it
needs nothing added to the experiment, but it does assume the list
applies to the assay and the assembly at hand.

`"librarySize"` and `"readsInRegions"` are the two depth-based options,
and they differ in what they call depth. The first takes the whole
library, so a sample where the immunoprecipitation worked better keeps
that advantage, which is what makes it the honest choice when the
question is how much signal each sample carries. The second takes the
reads collected in the regions, the quantity DiffBind calls reads in
peaks, so the samples are put on the same total enrichment before
anything is compared. That assumes the fraction of the library sitting
in the regions is a property of the protocol rather than of the biology,
and a treatment that genuinely redistributes a mark violates it in the
direction that hides the effect.

## See also

[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md),
[`countGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)

# Scaling factors from the background bins, which is what they are counted for
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes
SummarizedExperiment::colData(counts)
#> DataFrame with 4 rows and 9 columns
#>                                                 sample               bam.file
#>                                            <character>            <character>
#> lv-H3K4me3-BN-female-bio1-tech1 lv-H3K4me3-BN-female.. /home/s.gregoricchio..
#> lv-H3K4me3-BN-male-bio2-tech1   lv-H3K4me3-BN-male-b.. /home/s.gregoricchio..
#> lv-H3K4me3-SHR-male-bio2-tech1  lv-H3K4me3-SHR-male-.. /home/s.gregoricchio..
#> lv-H3K4me3-SHR-male-bio3-tech1  lv-H3K4me3-SHR-male-.. /home/s.gregoricchio..
#>                                 condition         sex biologicalReplicate
#>                                  <factor> <character>         <character>
#> lv-H3K4me3-BN-female-bio1-tech1       BN       female                bio1
#> lv-H3K4me3-BN-male-bio2-tech1         BN         male                bio2
#> lv-H3K4me3-SHR-male-bio2-tech1        SHR        male                bio2
#> lv-H3K4me3-SHR-male-bio3-tech1        SHR        male                bio3
#>                                 paired.end library.size norm.factor
#>                                  <logical>    <numeric>   <numeric>
#> lv-H3K4me3-BN-female-bio1-tech1      FALSE       386378    0.795932
#> lv-H3K4me3-BN-male-bio2-tech1        FALSE       400384    0.803658
#> lv-H3K4me3-SHR-male-bio2-tech1       FALSE       337600    0.921449
#> lv-H3K4me3-SHR-male-bio3-tech1       FALSE       780222    1.696608
#>                                 scaling.factor
#>                                      <numeric>
#> lv-H3K4me3-BN-female-bio1-tech1       0.543313
#> lv-H3K4me3-BN-male-bio2-tech1         0.568472
#> lv-H3K4me3-SHR-male-bio2-tech1        0.549586
#> lv-H3K4me3-SHR-male-bio3-tech1        2.338629

# TMM estimates them from the regions themselves instead
tmmCounts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
SummarizedExperiment::colData(tmmCounts)$scaling.factor
#> lv-H3K4me3-BN-female-bio1-tech1   lv-H3K4me3-BN-male-bio2-tech1 
#>                       0.7502489                       0.7881769 
#>  lv-H3K4me3-SHR-male-bio2-tech1  lv-H3K4me3-SHR-male-bio3-tech1 
#>                       0.6487672                       1.8128070 

if (FALSE) { # \dontrun{
# Factors coming from a spike-in, or from any external estimate
counts <- normalizeCounts(counts, scalingFactors = spikeInFactors)
} # }
```
