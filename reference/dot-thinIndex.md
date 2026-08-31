# .thinIndex

Returns the positions of a regularly spaced subset of a vector, used to
keep the point clouds drawable. The thinning is deterministic, so the
same object always gives the same picture.

## Usage

``` r
.thinIndex(n, maxPoints)
```

## Arguments

- n:

  Number of available elements.

- maxPoints:

  Maximum number of elements to keep.

## Value

An integer vector of positions.

## Author

Sebastian Gregoricchio
