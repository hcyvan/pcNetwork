library(devtools)
library(pryr)
library(roxygen2)
# library(testthat)
package.path <- './package/x2y'
# create(package.path)

roxygenize(package.path)
load_all(package.path)
document(package.path)


############################################################
setwd('/home/c509/文档/bio/pcNetwork')
source('./lib/globals.R')
source('./lib/helpers.R')
#prad.rna.count <- helper.get.data.count()
#prad.rna.fpkm <- helper.get.fpkm()

setwd('/home/c509/文档/bio/pcProfile')
#devtools::use_data(prad.rna.count, overwrite = TRUE)
#devtools::use_data(prad.rna.fpkm, overwrite = TRUE)
load_all('.')

