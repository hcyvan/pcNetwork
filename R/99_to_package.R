library(devtools)
library(pryr)
library(roxygen2)
# library(testthat)

# create('./package/x2y')
roxygenize('./package/x2y')
load_all('./package/x2y')
document('./package/x2y')


############################################################
setwd('/home/c509/文档/bio/pcNetwork')
source('./lib/globals.R')
source('./lib/helpers.R')

library(devtools)
library(pryr)
library(roxygen2)

roxygenize('../pcProfile')
load_all('../pcProfile')
document('../pcProfile')

#prad.rna.count <- helper.get.data.count()
#prad.rna.fpkm <- helper.get.fpkm()

setwd('/home/c509/文档/bio/pcProfile')
#devtools::use_data(prad.rna.count, overwrite = TRUE)
#devtools::use_data(prad.rna.fpkm, overwrite = TRUE)
load_all('.')


install_github('hcyvan/pcProfile')
