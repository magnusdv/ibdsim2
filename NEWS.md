# ibdsim2 1.3

## New features

* In `ibdsim()` allow `map` to be a list of chromosome maps.

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
