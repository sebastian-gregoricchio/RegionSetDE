# .resolveSheetPaths

Turns the relative paths of a sample sheet into absolute ones, leaving
alone the missing values, the absolute paths and the addresses starting
with a protocol.

## Usage

``` r
.resolveSheetPaths(paths, basePath)
```

## Arguments

- paths:

  Character vector with the paths.

- basePath:

  String with the directory the relative paths refer to.

## Value

A character vector with the resolved paths.

## Author

Sebastian Gregoricchio
