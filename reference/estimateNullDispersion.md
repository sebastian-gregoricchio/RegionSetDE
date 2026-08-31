# estimateNullDispersion

Estimates the biological variation between samples from a collection of
rows assumed not to respond to the contrast, so that a design with no
replicates has a dispersion to be tested against. The rows are usually
the background bins, which cover the genome and should carry no
treatment effect, but any region set believed to be invariant works the
same way.

## Usage

``` r
estimateNullDispersion(
  counts,
  source = "background",
  regionSets = NULL,
  index = NULL,
  samples = NULL,
  minCount = 10,
  maxRows = 50000,
  holdout = 0.5,
  subset = NULL,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- source:

  String with where the null rows come from, one of `"background"` (the
  bins stored by
  [`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)),
  `"regionSet"` (a set of the object) and `"supplied"`. Default:
  `"background"`.

- regionSets:

  Character vector with the names of the sets used as null rows. Only
  for `source = "regionSet"`. Default: `NULL`.

- index:

  Integer vector with the positions of the null rows. Only for
  `source = "supplied"`. Default: `NULL`.

- samples:

  Character vector with the samples the estimate is computed on.
  Default: `NULL`, all of them.

- minCount:

  Numeric value with the average count a null row must carry to be used.
  Default: `10`.

- maxRows:

  Numeric value with the number of null rows kept, drawn
  deterministically when there are more. Default: `50000`.

- holdout:

  Numeric value between 0 and 1 with the fraction of the null rows left
  out of the estimate, so that
  [`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
  has rows the dispersion has not already been fitted to. Default:
  `0.5`.

- subset:

  Integer vector restricting the null rows to a subset of the ones
  `source` selects. Default: `NULL`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A list with the `dispersion`, its square root as `bcv`, the `source`,
the number of rows it was computed on, the samples used, and
`holdout.index`, the rows kept aside for the calibration check.

## Details

Without replicates nothing in the data measures how much two libraries
differ for reasons that have nothing to do with the treatment, and every
engine either refuses to fit or invents an answer. The way out is to
assume that some rows do not respond, and to read the variation across
those rows as the variation that would be seen between replicates. The
estimate is a common dispersion fitted under an intercept-only model,
which is exactly what `edgeR` recommends for an experiment with no
replication, with the housekeeping genes replaced here by the background
bins.

That assumption is the whole estimate, so it is worth being deliberate
about it. Background bins are the safest choice, since a treatment that
changed the genome-wide average would have broken the normalisation long
before it reached this point. A region set chosen because it looked flat
in the data is the unsafe one: picking rows for their small fold change
and then measuring the spread of those fold changes gives a dispersion
biased towards zero, and p-values that follow it down.

A plausible number is not a replicate. Everything downstream stays
conditional on this estimate being right, which is why
[`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
exists and why it should be run before any of the output is believed.
Checking the estimate against the rows it was fitted to would say
nothing, so half the null rows are held out by default and travel back
in `holdout.index` for that check to use.

## See also

[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
[`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md),
[`countBackground`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

# The background bins cover the genome and carry no strain effect
nullDispersion <- estimateNullDispersion(counts, source = "background",
                                         verbose = FALSE)
nullDispersion
#> $dispersion
#> [1] 0.06820971
#> 
#> $bcv
#> [1] 0.2611699
#> 
#> $source
#> [1] "background"
#> 
#> $n.rows
#> [1] 560
#> 
#> $holdout.index
#>   [1]   11   27   31   37   46   48   53   55   66   87   90   92   94   97   99
#>  [16]  104  106  108  111  113  115  117  119  121  123  125  127  131  133  135
#>  [31]  137  139  141  143  145  163  165  167  169  171  173  175  177  180  183
#>  [46]  185  187  190  195  198  200  203  207  209  211  213  215  217  219  222
#>  [61]  224  226  228  230  233  235  237  239  241  243  245  247  249  251  253
#>  [76]  255  258  261  263  265  267  269  271  273  275  277  279  281  283  285
#>  [91]  287  290  292  294  296  299  302  304  306  308  311  318  320  322  330
#> [106]  332  334  336  338  340  343  346  348  355  357  361  363  365  367  369
#> [121]  372  375  377  379  381  383  386  389  392  395  397  400  402  404  408
#> [136]  410  412  414  416  419  421  423  425  427  429  431  433  439  442  447
#> [151]  451  455  457  459  461  463  465  467  470  472  474  477  481  483  486
#> [166]  488  490  492  494  496  498  500  502  504  506  508  510  512  514  516
#> [181]  518  520  522  524  526  528  530  532  535  537  540  542  544  546  548
#> [196]  550  553  555  557  559  561  565  567  569  571  580  583  592  601  604
#> [211]  606  609  615  617  620  662  664  666  668  670  685  687  689  691  693
#> [226]  696  699  702  704  706  708  710  712  714  716  718  720  722  724  726
#> [241]  728  730  733  735  737  739  741  743  745  747  750  754  756  758  760
#> [256]  762  764  766  768  770  772  774  777  779  781  789  791  793  797  804
#> [271]  807  811  813  816  818  821  825  827  829  831  832  834  840  842  844
#> [286]  846  848  850  852  854  856  858  862  864  866  868  871  874  877  885
#> [301]  894  899  903  905  907  909  911  913  915  917  919  921  923  925  928
#> [316]  930  934  936  938  940  942  946  950  952  954  956  958  960  962  964
#> [331]  968  971  974  978  981  990  997  999 1001 1005 1025 1028 1035 1043 1049
#> [346] 1061 1067 1086 1088 1090 1092 1094 1096 1098 1100 1102 1104 1106 1108 1110
#> [361] 1112 1114 1116 1118 1120 1122 1124 1126 1128 1130 1132 1134 1136 1138 1141
#> [376] 1143 1145 1147 1149 1151 1153 1155 1157 1159 1161 1163 1165 1167 1169 1171
#> [391] 1173 1175 1177 1179 1181 1183 1185 1187 1189 1191 1193 1195 1197 1199 1201
#> [406] 1203 1205 1207 1209 1211 1213 1215 1217 1219 1221 1223 1225 1227 1230 1232
#> [421] 1234 1237 1239 1242 1244 1252 1255 1257 1259 1262 1264 1266 1270 1272 1274
#> [436] 1276 1278 1280 1282 1284 1288 1291 1293 1296 1298 1300 1303 1305 1307 1309
#> [451] 1312 1314 1316 1319 1322 1324 1326 1328 1330 1334 1336 1338 1340 1343 1346
#> [466] 1354 1356 1358 1360 1363 1365 1367 1369 1371 1374 1381 1387 1391 1395 1397
#> [481] 1399 1401 1403 1405 1409 1415 1418 1420 1422 1424 1429 1431 1433 1435 1437
#> [496] 1439 1441 1443 1445 1447 1449 1451 1453 1455 1457 1459 1461 1463 1465 1467
#> [511] 1469 1471 1473 1475 1477 1479 1481 1483 1485 1487 1489 1491 1494 1496 1503
#> [526] 1505 1507 1509 1511 1513 1515 1519 1521 1523 1525 1528 1530 1532 1534 1536
#> [541] 1538 1540 1542 1544 1546 1548 1550 1552 1554 1556 1558 1561 1563 1565 1567
#> [556] 1569 1571 1573 1575 1577
#> 
#> $samples
#> [1] "lv-H3K4me3-BN-female-bio1-tech1" "lv-H3K4me3-BN-male-bio2-tech1"  
#> [3] "lv-H3K4me3-SHR-male-bio2-tech1"  "lv-H3K4me3-SHR-male-bio3-tech1" 
#> 

# Any set believed to be invariant works the same way, as long as it has signal
fromRegions <- estimateNullDispersion(counts,
                                      source = "regionSet",
                                      regionSets = "geneBody",
                                      verbose = FALSE)
fromRegions
#> $dispersion
#> [1] 0.09171431
#> 
#> $bcv
#> [1] 0.3028437
#> 
#> $source
#> [1] "regionSet"
#> 
#> $n.rows
#> [1] 92
#> 
#> $holdout.index
#>  [1]   38   61   82   85  102  109  124  128  132  152  163  193  196  199  218
#> [16]  222  242  246  257  265  279  315  331  337  340  358  393  405  408  418
#> [31]  450  459  469  478  528  555  561  563  592  616  619  646  655  660  666
#> [46]  672  682  689  763  773  775  782  805  817  846  891  907  926  931  938
#> [61]  947  963  981  985  990 1009 1018 1020 1030 1038 1047 1059 1073 1104 1110
#> [76] 1129 1168 1178 1184 1195 1199 1227 1229 1240 1242 1253 1258 1261 1282 1291
#> [91] 1300 1330
#> 
#> $samples
#> [1] "lv-H3K4me3-BN-female-bio1-tech1" "lv-H3K4me3-BN-male-bio2-tech1"  
#> [3] "lv-H3K4me3-SHR-male-bio2-tech1"  "lv-H3K4me3-SHR-male-bio3-tech1" 
#> 
```
