# .makeParallelParam

Builds the `BiocParallel` back end matching the number of requested
threads and the operating system.

## Usage

``` r
.makeParallelParam(nThreads = 1)
```

## Arguments

- nThreads:

  Number of threads. Default: `1`.

## Value

A `BiocParallelParam` object.

## Author

Sebastian Gregoricchio
