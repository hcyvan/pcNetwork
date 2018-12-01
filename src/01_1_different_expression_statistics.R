source('./lib/globals.R')
source('./lib/helpers.R')

data.fpkm <- helper.get.fpkm.count()


diff.gene <- helper.get.diff.gene()
pcg <- helper.get.PCG()
lncRNA <- helper.get.lncRNA()
dim(diff.gene)
dim(pcg)
dim(lncRNA)
table(diff.gene$GeneType)
config$lncRNA
config$PCGs

show()

