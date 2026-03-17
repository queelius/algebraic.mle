# Collect a field from a list of MLEs, returning NULL if any are NULL.

Collect a field from a list of MLEs, returning NULL if any are NULL.

## Usage

``` r
collect_or_null(mles, extractor)
```

## Arguments

- mles:

  List of mle_fit objects.

- extractor:

  Function that extracts the field from a single MLE.

## Value

A list of field values, or NULL if any value is NULL.
