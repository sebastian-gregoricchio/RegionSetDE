# .makeParallelParam

Builds the `BiocParallel` back end matching the number of requested
threads and the operating system.

## Usage

``` r
.makeParallelParam(nThreads = 1, tasks = 0L)
```

## Arguments

- nThreads:

  Number of threads. Default: `1`.

- tasks:

  Number of tasks the work is split into, see
  [`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html).
  Setting it to the number of jobs hands the jobs out one at a time, as
  the threads become free. Default: `0`, one task per thread.

## Value

A `BiocParallelParam` object.

## Author

Sebastian Gregoricchio
