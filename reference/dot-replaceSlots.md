# .replaceSlots

Writes new values into the slots of an object, the way
`object@slot <- value` would, with the class of every value checked
against its slot. The tables given as data.frames are converted for the
slots holding a `DataFrame`. Validity is left to the caller, since
several slots may have to change before the object is consistent again.

## Usage

``` r
.replaceSlots(object, ...)
```

## Arguments

- object:

  An object of one of the classes of the package.

- ...:

  New values, named after their slots.

## Value

The object with the new values.

## Author

Sebastian Gregoricchio
