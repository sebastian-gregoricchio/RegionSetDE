# .bplapplySameLibraries

Runs
[`BiocParallel::bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
with socket workers that load the packages from the same libraries as
the calling session. A socket worker is a new R process: the functions
it receives refer to the RegionSetDE namespace by name, and the worker
loads that namespace from the libraries it starts with, which are not
the ones the session may have added with
[`.libPaths()`](https://rdrr.io/r/base/libPaths.html). When the two
disagree, as during `R CMD build`, where the package being built sits in
a temporary library, the worker runs new code against an older installed
namespace, or finds no package at all. Forked and serial workers share
the session and need nothing.

## Usage

``` r
.bplapplySameLibraries(X, FUN, ..., BPPARAM)
```

## Arguments

- X:

  List or vector to iterate over.

- FUN:

  Function applied to each element.

- ...:

  Further arguments of `FUN`.

- BPPARAM:

  `BiocParallelParam` object.

## Value

The list returned by
[`BiocParallel::bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html).

## Author

Sebastian Gregoricchio
