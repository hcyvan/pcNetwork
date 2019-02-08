
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


setwd('../pcProfile')
tf2gene.jaspar <- fimo.tss.set.2
tf2gene.gtrd <- gtrd.set.3
tf2gene.trrust <- trrust.set.2
#devtools::use_data(prad.rna.count, overwrite = TRUE)
#devtools::use_data(prad.rna.fpkm, overwrite = TRUE)
devtools::use_data(tf2gene.jaspar, overwrite = TRUE)
devtools::use_data(tf2gene.trrust, overwrite = TRUE)
devtools::use_data(tf2gene.gtrd, overwrite = TRUE)



roxygenize('.')
load_all('.')


install_github('hcyvan/pcProfile')
