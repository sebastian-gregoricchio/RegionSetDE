# .fryValues

Returns the values the self-contained test is run on: the log values of
the fit for the linear engines, and the z-scores of the counts under the
null model of the fit for `edgeR` and `DESeq2`.

## Usage

``` r
.fryValues(fit, contrastVector, fryInput = "auto")
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrastVector:

  Numeric vector with the contrast, in the columns of the design.

- fryInput:

  String, either `"auto"` or `"logcpm"`. Default: `"auto"`.

## Value

A list with `values`, a numeric matrix with one row per region and one
column per sample, `standardize`, the value to pass to
[`limma::fry`](https://rdrr.io/pkg/limma/man/roast.html), and `input`, a
label of what the values are.

## Details

The z-scores follow `edgeR`: the null model is the design with the
contrast taken out, fitted with the offsets and the dispersion of every
region, and each count is turned into the normal deviate of its negative
binomial probability under that model
([`edgeR::zscoreNBinom`](https://rdrr.io/pkg/edgeR/man/zscoreNBinom.html)).
They are computed once for the contrast and shared by all the sets,
where `fry` on a `DGEList` would refit the whole null model for every
set. For `DESeq2` the offsets are the logarithm of its normalisation or
size factors and the dispersions its final ones, so the null model is
the one `DESeq2` fitted, estimated here by
[`edgeR::glmFit`](https://rdrr.io/pkg/edgeR/man/glmfit.html).

## Author

Sebastian Gregoricchio
