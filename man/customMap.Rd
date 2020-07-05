% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_utils.R
\name{customMap}
\alias{customMap}
\title{Custom recombination map}
\usage{
customMap(x)
}
\arguments{
\item{x}{A data frame (or matrix), whose column names must include "chrom",
"mb", and either "cm" or both "male" and "female". (See Examples.)}
}
\value{
An object of class \code{genomeMap}
}
\description{
Create custom recombination maps for use in \code{\link[=ibdsim]{ibdsim()}}.
}
\examples{
# A map including two chromosomes.
df1 = data.frame(chrom = c(1, 1, 2, 2),
                 mb = c(0, 2, 0, 5),
                 cm = c(0, 3, 0, 6))
map1 = customMap(df1)
map1

# Use columns "male" and "female" to make sex specific maps
df2 = data.frame(chrom = c(1, 1, 2, 2),
                 mb = c(0, 2, 0, 5),
                 male = c(0, 3, 0, 6),
                 female = c(0, 4, 0, 7))
map2 = customMap(df2)
map2


}
\seealso{
\code{\link[=uniformMap]{uniformMap()}}, \code{\link[=loadMap]{loadMap()}}
}