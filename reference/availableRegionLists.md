# availableRegionLists

Lists the blacklists and greenlists shipped with the package, with the
assembly, the assay, the source and the size of each of them.

## Usage

``` r
availableRegionLists(type = NULL, genome = NULL)
```

## Arguments

- type:

  String with the lists returned, either `"blacklist"` or `"greenlist"`.
  Default: `NULL`, both.

- genome:

  String with an assembly the lists are restricted to. Default: `NULL`,
  all of them.

## Value

A data.frame with one row per list: the `type`, the `genome`, the
`assay` it was built for or `"any"`, the `source`, the `version`, the
number of regions, the base pairs they cover and the reference.

## See also

[`loadBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadBlacklist.md),
[`loadGreenlist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
availableRegionLists()
#>         type genome  assay        source version n.regions covered.bp
#> 1  blacklist   ce10    any        ENCODE      v2       100    2205200
#> 2  blacklist   ce11    any        ENCODE      v2        97     728500
#> 3  blacklist    dm3    any        ENCODE      v2       271    2689400
#> 4  blacklist    dm6    any        ENCODE      v2       182    3740300
#> 5  blacklist   hg19    any        ENCODE      v2       834  274970000
#> 6  blacklist   hg38    any        ENCODE      v2       636  227162400
#> 7  blacklist   hg38 cutrun       deMello      v1       832   10133493
#> 8  blacklist   hg38 cuttag       deMello      v1      2020    9183890
#> 9  blacklist    hs1    any excluderanges      v1      3565  275454700
#> 10 blacklist   mm10    any        ENCODE      v2      3435  238977200
#> 11 blacklist   mm39 cutrun       deMello      v1      1452   22330024
#> 12 greenlist   hg38 cutrun       deMello      v1       869    1725126
#> 13 greenlist   hg38 cuttag       deMello      v1      2767    3811775
#> 14 greenlist   mm39 cutrun       deMello      v1      1694    6493000
#>                                            reference
#> 1   Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 2   Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 3   Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 4   Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 5   Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 6   Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 7  de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 8  de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 9   Dozmorov et al. (2023) Bioinformatics 39:btad198
#> 10  Amemiya, Kundaje and Boyle (2019) Sci Rep 9:9354
#> 11 de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 12 de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 13 de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 14 de Mello et al. (2024) Brief Bioinform 25:bbad538

availableRegionLists(type = "greenlist")
#>        type genome  assay  source version n.regions covered.bp
#> 1 greenlist   hg38 cutrun deMello      v1       869    1725126
#> 2 greenlist   hg38 cuttag deMello      v1      2767    3811775
#> 3 greenlist   mm39 cutrun deMello      v1      1694    6493000
#>                                           reference
#> 1 de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 2 de Mello et al. (2024) Brief Bioinform 25:bbad538
#> 3 de Mello et al. (2024) Brief Bioinform 25:bbad538
```
