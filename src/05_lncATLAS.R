source('./lib/globals.R')
source('./lib/helpers.R')

# Data.Type
# ratio2=log2(cytosol/nucleus)
#

atlas.raw <- read.csv('./data/2018-12-17_lncATLAS_all_data.csv')
atlas <- atlas.raw%>%filter(Data.Type=='ratio2')

atlas.ratio2.mean <- sapply(split(atlas,as.vector(atlas$ENSEMBL.ID)), function(g){
  mean(g$Value,na.rm=T)
})
atlas.ratio2.mean<-atlas.ratio2.mean[!is.na(atlas.ratio2.mean)]

loc <- data.frame(id=names(atlas.ratio2.mean), rci=atlas.ratio2.mean)
saveRDS(loc,file = './cache/cell.location.rds')
