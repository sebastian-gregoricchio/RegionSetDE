# .elapsedTime

Writes the time gone by since a starting point in the unit that reads
best: seconds, minutes or hours.

## Usage

``` r
.elapsedTime(startTime)
```

## Arguments

- startTime:

  Value returned by [`Sys.time()`](https://rdrr.io/r/base/Sys.time.html)
  when the work started.

## Value

A string such as `"42 s"`, `"3.5 min"` or `"1.2 h"`.

## Author

Sebastian Gregoricchio
