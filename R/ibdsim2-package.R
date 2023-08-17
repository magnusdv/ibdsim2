#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import pedtools
#' @importFrom stats rpois runif

#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib ibdsim2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onUnload <- function (libpath) {
  library.dynam.unload('ibdsim2', libpath)
}