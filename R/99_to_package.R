library(devtools)
library(roxygen2)
# library(testthat)
package.path <- './package/x2y'
# create(package.path)

roxygenize(package.path)
load_all(package.path)
document(package.path)
