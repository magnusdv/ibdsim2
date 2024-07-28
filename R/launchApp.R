#' Launch the ibdsim2 app
#'
#' This launches the Shiny app for simulating IBD segment distributions.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#'
#' \dontrun{
#' launchApp()
#' }
#'
#' @export
launchApp = function() {
  packages = c("shiny", "shinyjs", "lubridate")
  for(pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      msg = sprintf("Package '%s' is required but not installed.\nPlease run `install.packages('%s')` and try again", pkg, pkg)
      stop2(msg)
    }
  }
  shiny::runApp(system.file("shiny", package = "ibdsim2"))
}
