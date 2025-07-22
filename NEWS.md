# ibdsim2 2.3.0

## Breaking changes

* Update the built-in recombination map "decode19". This is still a thinned version of the map from Halldorsson et al. (2019), using GRCh38 coordinates. The new version includes physical chromosome endpoints, which were previously partially missing. It also uses a better thinning strategy, which allowed a reduction of data points from ~38 000 to ~14 000 without losing accuracy.

## App changes

* The app now uses the updated built-in map "decode19" (see above). This may lead to slightly different results when the length unit is set to ""

## New features

* The `realised...()` functions gain argument `simplify1` for simplifying the output when `N = 1`. This is useful for interactive use.

* Fixed a bug in `profileSimIBD()`, which sometimes affected markers in the telomeric regions.
* Added internal helper functions for taking unions and intersections of IBD segments.

# ibdsim2 2.2.0

## App updates

* `launchApp()` now opens the app in the default browser.
* The app now stops automatically when the browser is closed (locally).
* Minor appearance tweaks.

## Other

* Added `merge = TRUE` argument to `realisedKappa()` and other `realised...` functions. Previously, IBD segment merging was handled inconsistently, sometimes causing unexpected results. (Thanks, @mkruijver.)
* Added all required packages for the app to 'Suggests'.
* Remove unused `shinyBS` import


# ibdsim2 2.1.1

## App updates

* Tweaked the labels of the built-in pedigrees (some were confusingly named). 
* Added a few more inbred built-in pedigrees.
* Revised the default individuals to be selected for each built-in pedigree.
* Changed the pedigree shown at startup to half siblings.

## Package updates

* Fixed a bug in `haploDraw()` causing an additional empty plot in some situations.


# ibdsim2 2.1.0

This version includes a major update of the shiny app frontend to **ibdsim2**. Previously developed in a separate repository (accessible at https://github.com/magnusdv/ibdsim2-legacy/), the app is now included as part of the ibdsim2 package, and can be run locally with `ibdsim2::launchApp()`. The live version is available at https://magnusdv.shinyapps.io/ibdsim2-shiny/.

## New app features

* The simulations are now much faster than before. As a result, the default number of simulations has been increased from 50 to 500.

* X-chromosomal IBD simulations are now supported.

* The user can now choose the length unit for IBD segments; either centiMorgan ("cM", default) or megabytes ("Mb").

* Fixed buggish unit conversion in the previous version: In some cases the segments were measured in Mb while the plot labels said "cM".

* The random number seed can now be selected by the user.

* More coherent layout and better pedigree plots.

* New function `launchApp()` for running the app from within R.

## Other 

* `findPattern()` now works as intended also for X-chromosomal simulations. 
* Various updates of docs and examples.


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
