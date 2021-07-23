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
