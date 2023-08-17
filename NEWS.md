# ibdsim2 2.0.0

## Breaking changes

* As of version 2.0.0, the main data structure for IBD segments includes both megabase (MB) and centiMorgan (CM) coordinates. This allows the user to choose length unit in all downstream analyses. However, the new format is not compatible with previous versions of the package.

* `ibdsim()` has a new argument `simplify1`, by default TRUE. This means that `ibdsim(..., N = 1)` now simply returns a matrix, without the outer list layer. This typically the desired behaviour in interactive use, especially when piping. To enforce a list output, add `simplify1 = FALSE`.


## New features

* All downstream functions depending on segment lengths (e.g. `segmentStats()` and `plotSegmentDistribution()`) have a new argument `unit`, allowing the user to choose between "cm" (centiMorgan) and  "mb" (megabases).

* `haploDraw()` has nicer default colours and automatically produces sensible margins.

* New (experimental) function `karyoHaploid()` for visualising IBD segments in karyogram plots.

* `convertPos()` has been rewritten using Rcpp, making it much more efficient.

* General overhaul of documentation, examples and README.


# ibdsim2 1.5.0

## New features

* X-chromosomal simulations are now implemented.

* `haploDraw()` now handles and displays X-chromosomal simulations.

## Other

* Fix labelling bug in `haploDraw()`.

* `profileSimIBD()` has been overhauled, fixing several glitches and with significant speed improvements. Note: Simulations with a given `seed` may differ across versions.

* The internal dataset `decode19` has been recompiled, updating some attributes. This should not affect regular users.


# ibdsim2 1.4.0

## New features

* New function `segmentStats()` for summarising the segments identified by `findPattern()`.

* `findPattern()` has a new argument `cutoff` for excluding short segments.

* `findPattern()` can detect more patterns: the argument `pattern` now accepts entries `autozygous` and `heterozygous`.

* `findPattern()` is much faster now, mostly due to a new implementation of the internal `mergeSegments()`

* New function `extractIds()` replaces the previous (non-exported) `extractIdsFromSegmentSummary()`.

* Updated README, including links to Shiny app.


## Bug fixes

* Fixed wrong length attribute of `loadMap()` when `uniform = F` and `sexAverage = T`.


# ibdsim2 1.3.0

## New features

* `ibdsim()` now allows `map` to be a list of chromosome maps.

* New function `findPattern()` for identifying IBD patterns in simulation outputs.

* New function `convertPos()` for converting between genetic and physical positions.

* New functions `mapLen()` and `physRange()` for retrieving map info.

* `profileSimIBD()` now sorts the genotypes before returning.

* README is updated and expanded.


## Bug fixes

* Fix bug in `estimateTwoLocusInbreeding()`

* Fix issue with negative segment lengths


# ibdsim2 1.2

## Breaking changes

* The previous built-in recombination map (based on Kong et al., 2010) is 
replaced with one based on Halldorsson et al. (2019), using GRCh38 coordinates. 
For reasons of speed and efficiency, the built-in map is a thinned version of 
the published map, keeping around 60 000 data points.

## New features

* New function `haploDraw()` for visualising IBD patterns

* The internal data structure for genetic maps has been rewritten. 

* New functions `loadMap()` and `customMap()` for loading built-in and 
user-prepared maps.

## Bug fixes

* Fixed CRAN warnings about self-assignment in Rcpp code.


# ibdsim2 1.1

* Initial CRAN release
