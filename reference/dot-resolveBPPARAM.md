# .resolveBPPARAM

Turns the `BPPARAM` argument of the user into a `BiocParallel` back end.
A `BiocParallelParam` object is returned untouched, a number of cores
builds the back end matching the operating system, and `NULL` gives a
serial one.

## Usage

``` r
.resolveBPPARAM(BPPARAM = NULL, tasks = 0L, progressBar = FALSE)
```

## Arguments

- BPPARAM:

  A `BiocParallelParam` object, a number of cores, or `NULL`. Default:
  `NULL`.

- tasks:

  Number of tasks the work is split into when the back end is built
  here, see
  [`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html).
  Setting it to the number of jobs hands the jobs out one at a time, as
  the workers become free. Default: `0`, one task per worker.

- progressBar:

  Logical value to indicate whether a back end built here must draw its
  progress bar, which advances by one step for every task that comes
  back. Default: `FALSE`.

## Value

A `BiocParallelParam` object.

## Author

Sebastian Gregoricchio
