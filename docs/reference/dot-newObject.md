# .newObject

Creates an object of one of the classes of the package, as
[`methods::new`](https://rdrr.io/r/methods/new.html) does, converting
the tables given as data.frames for the slots holding a `DataFrame`.

## Usage

``` r
.newObject(Class, ...)
```

## Arguments

- Class:

  String with the class of the object.

- ...:

  Values of the slots, named after them. An unnamed value, such as the
  `RangedSummarizedExperiment` a counts object extends, is passed on as
  it is.

## Value

The new object.

## Author

Sebastian Gregoricchio
