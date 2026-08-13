# hmagmaR 0.99.3

* Added this NEWS file.

# hmagmaR 0.99.2

* Fixed `DESCRIPTION`/`NAMESPACE` inconsistencies: removed unused `data.table`
  and `dplyr` imports, added proper imports for `org.Hs.eg.db` and `stats`.
* Removed a duplicate data file that redefined the `snpgeneexon` object,
  which was causing an installation-time conflict.
* Documented the `regulatoryRegions` object and moved a raw reference table
  out of `data/` and into `inst/extdata/`.
* Added runnable `@examples` to all exported functions.
* `SampledDownAnnotation()` no longer leaks changes to the caller's random
  number generator state after it returns.
* Removed `LazyData: true`, per Bioconductor guidance.

# hmagmaR 0.99.1

* Fixed missing `IRanges` and `S4Vectors` entries in `DESCRIPTION` Imports.

# hmagmaR 0.99.0

* Initial Bioconductor submission.
