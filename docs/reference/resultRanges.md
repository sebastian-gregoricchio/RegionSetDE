# resultRanges

Returns the coordinates of a `RegionSetDE.results` object with the
statistics attached as metadata columns, ready to be exported as a
BED-like file.

## Usage

``` r
resultRanges(results)

# S4 method for class 'RegionSetDE.results'
resultRanges(results)

# S4 method for class 'RegionSetDE.resultsList'
resultRanges(results)
```

## Arguments

- results:

  `RegionSetDE.results` object.

## Value

A `GRanges` with one element per region.

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

resultRanges(results)
#> GRanges object with 1895 ranges and 15 metadata columns:
#>                               seqnames            ranges strand |
#>                                  <Rle>         <IRanges>  <Rle> |
#>   promoterNonCpG|region_00012    chr12       26988-27987      * |
#>   promoterNonCpG|region_00017    chr12       39449-40448      * |
#>   promoterNonCpG|region_00019    chr12       44116-45115      * |
#>   promoterNonCpG|region_00020    chr12       46527-47526      * |
#>   promoterNonCpG|region_00026    chr12       89229-90228      * |
#>                           ...      ...               ...    ... .
#>      promoterCpG|region_03781    chr12 46595293-46596292      * |
#>      promoterCpG|region_03789    chr12 46646773-46647772      * |
#>      promoterCpG|region_03792    chr12 46684188-46685187      * |
#>      promoterCpG|region_03797    chr12 46721630-46722629      * |
#>      promoterCpG|region_03798    chr12 46728313-46729312      * |
#>                                   region.set    region.id   tile.id     log2FC
#>                                  <character>  <character> <integer>  <numeric>
#>   promoterNonCpG|region_00012 promoterNonCpG region_00012      <NA> -0.4744875
#>   promoterNonCpG|region_00017 promoterNonCpG region_00017      <NA> -1.6911510
#>   promoterNonCpG|region_00019 promoterNonCpG region_00019      <NA>  0.1328365
#>   promoterNonCpG|region_00020 promoterNonCpG region_00020      <NA>  1.3000608
#>   promoterNonCpG|region_00026 promoterNonCpG region_00026      <NA> -0.0791789
#>                           ...            ...          ...       ...        ...
#>      promoterCpG|region_03781    promoterCpG region_03781      <NA>  -0.912532
#>      promoterCpG|region_03789    promoterCpG region_03789      <NA>  -1.066691
#>      promoterCpG|region_03792    promoterCpG region_03792      <NA>  -0.999014
#>      promoterCpG|region_03797    promoterCpG region_03797      <NA>  -1.278943
#>      promoterCpG|region_03798    promoterCpG region_03798      <NA>  -0.638963
#>                               average.signal average.signal.BN
#>                                    <numeric>         <numeric>
#>   promoterNonCpG|region_00012        3.09870           3.15461
#>   promoterNonCpG|region_00017        3.14197           3.51271
#>   promoterNonCpG|region_00019        3.27334           3.14980
#>   promoterNonCpG|region_00020        3.49069           2.93099
#>   promoterNonCpG|region_00026        3.22088           3.15220
#>                           ...            ...               ...
#>      promoterCpG|region_03781       10.11622          10.31085
#>      promoterCpG|region_03789       10.15290          10.41427
#>      promoterCpG|region_03792       10.14625          10.37415
#>      promoterCpG|region_03797       10.17313          10.51121
#>      promoterCpG|region_03798        7.03435           7.03222
#>                               average.signal.SHR       stat stat.distribution
#>                                        <numeric>  <numeric>       <character>
#>   promoterNonCpG|region_00012            3.07110 0.16629041                 f
#>   promoterNonCpG|region_00017            2.77439 2.33758830                 f
#>   promoterNonCpG|region_00019            3.39054 0.01002874                 f
#>   promoterNonCpG|region_00020            3.85473 2.08486536                 f
#>   promoterNonCpG|region_00026            3.29321 0.00711448                 f
#>                           ...                ...        ...               ...
#>      promoterCpG|region_03781            9.89195    0.95576                 f
#>      promoterCpG|region_03789            9.83439    1.26881                 f
#>      promoterCpG|region_03792            9.87635    1.13227                 f
#>      promoterCpG|region_03797            9.73201    1.80852                 f
#>      promoterCpG|region_03798            7.03908    1.80154                 f
#>                                     df1       df2   p.value       FDR
#>                               <numeric> <numeric> <numeric> <numeric>
#>   promoterNonCpG|region_00012         1   18.9454  0.688001  0.920086
#>   promoterNonCpG|region_00017         1   18.7356  0.142993  0.799720
#>   promoterNonCpG|region_00019         1   18.4527  0.921310  0.980729
#>   promoterNonCpG|region_00020         1   20.6420  0.163773  0.799720
#>   promoterNonCpG|region_00026         1   20.5652  0.933596  0.980729
#>                           ...       ...       ...       ...       ...
#>      promoterCpG|region_03781         1   18.9037  0.340604   0.79972
#>      promoterCpG|region_03789         1   18.9037  0.274087   0.79972
#>      promoterCpG|region_03792         1   18.9037  0.300696   0.79972
#>      promoterCpG|region_03797         1   18.9037  0.194600   0.79972
#>      promoterCpG|region_03798         1   18.8988  0.195429   0.79972
#>                               diff.status     regionId
#>                                  <factor>  <character>
#>   promoterNonCpG|region_00012        null region_00012
#>   promoterNonCpG|region_00017        null region_00017
#>   promoterNonCpG|region_00019        null region_00019
#>   promoterNonCpG|region_00020        null region_00020
#>   promoterNonCpG|region_00026        null region_00026
#>                           ...         ...          ...
#>      promoterCpG|region_03781        null region_03781
#>      promoterCpG|region_03789        null region_03789
#>      promoterCpG|region_03792        null region_03792
#>      promoterCpG|region_03797        null region_03797
#>      promoterCpG|region_03798        null region_03798
#>   -------
#>   seqinfo: 1 sequence from rn4 genome
```
