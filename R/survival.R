library(stringr)

source('./lib/globals.R')
source('./lib/helpers.R')
data.fpkm <- helper.get.fpkm.count()

slc45a2 <- 'ENSG00000164175'
unc5a <- 'ENSG00000113763'
tumor <- data.fpkm[,1:499]
tumor.2 <- tumor[match(c(slc45a2,unc5a),rownames(tumor)),]
final <- data.frame(t(tumor.2))

apply(str_split_fixed(rownames(final), '\\.', 4)[,1:3], 1,function(x){str_c(x, collapse = '-')})->sample

final <- data.frame(case=sample, final)

write.csv(final, file = 'two.gene.fpkm.csv', row.names = F)
