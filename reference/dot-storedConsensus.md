# .storedConsensus

Internal accessor of the consensus data of a `RegionSetDE` object.
Unlike
[`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md)
it does not stop when the regions were not built from peaks, and it
answers an empty list for any other object.

## Usage

``` r
.storedConsensus(object)
```

## Arguments

- object:

  Any object.

## Value

The list returned by
[`consensusData`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md),
or an empty list.

## Author

Sebastian Gregoricchio
