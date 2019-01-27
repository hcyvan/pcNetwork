source('./lib/helpers.R')
source('./R/lib.R')
devtools::load_all('./package/x2y/')


############################################### GET x.2.y
# plot.lncRNA.2.PCG()
# plot.tf.2.PCG()
################################################# tf-pcg
# get.tf.2.PCG.from.enricher <- function() {
#   tf.enricher <- read.delim('./data/enricher.all.bg0.new.csv', sep = ',', stringsAsFactors = F)
#   tf.enricher.list <- setNames(split(tf.enricher, seq(nrow(tf.enricher))), tf.enricher$detail.ID)
#
#   load('./cache/biomart.symbol.biotype.rda')
#
#   lapply(tf.enricher.list, function(tf){
#     symbol <- str_split(tf$detail.geneID, '/')[[1]]
#     filter(biomart.symbol.biotype, hgnc_symbol%in%symbol)$ensembl_gene_id
#   })
# }

#############################################################


pcg <- pf.get.diff('pcg')
lncRNA <- pf.get.diff('lncRNA')

biomart <- pf.get.biomart()
cor.pairs <- readRDS('./cache/cor.pearson.pairs.rds')

t1 <- biomart[match(cor.pairs$v1, biomart$ensembl_gene_id), 'gene_biotype']
t2 <- biomart[match(cor.pairs$v2, biomart$ensembl_gene_id), 'gene_biotype']
pcg2lncRNA <- filter(data.frame(cor.pairs,t1=t1,t2=t2), (t1%in%pv.pcg & t2%in%pv.lncRNA))[,1:2] # n:0
lncRNA2pcg.pairs <- filter(data.frame(cor.pairs,t1=t1,t2=t2), (t1%in%pv.lncRNA & t2%in%pv.pcg), FDR<0.05, abs(r) >= 0.3)

fimo.gss.list <- readRDS('./cache/fimo.gss.list.rds')
fimo.tss.pairs <- readRDS('./cache/fimo.tss.set.rds')
trrust.pairs <- readRDS('./cache/trrust.set.rds')
gtrd.pairs <- readRDS('./cache/gtrd.set.rds')


getAdjust <- function(xa, xb, background=NULL) {
  a.m <- x2yMatrix(xa, b=background)
  b.m <- x2yMatrix(xb, b=background)
  fix <- x2yMatrixAdjust(a.m, b.m, y.names = pcg$GeneID)
  fix
}


fix.fimo.gss <-getAdjust(lncRNA2pcg.pairs, fimo.gss.list, pcg$GeneID)
fix.fimo.tss <-getAdjust(lncRNA2pcg.pairs, fimo.tss.pairs, pcg$GeneID)
fix.trrust <-getAdjust(lncRNA2pcg.pairs, trrust.pairs, pcg$GeneID)
fix.gtrd <-getAdjust(lncRNA2pcg.pairs, gtrd.pairs, pcg$GeneID)

system.time(lnc2tf.fimo.gss <- xyCor(fix.fimo.gss$a, fix.fimo.gss$b, cores=6))#4
saveRDS(lnc2tf.fimo.gss, file = './cache/lnc2tf.fimo.gss.rds')
system.time(lnc2tf.fimo.tss <- xyCor(fix.fimo.tss$a, fix.fimo.tss$b, cores=6))#3
saveRDS(lnc2tf.fimo.tss, file = './cache/lnc2tf.fimo.tss.rds')
system.time(lnc2tf.trrust <- xyCor(fix.trrust$a, fix.trrust$b, cores=6))#4
saveRDS(lnc2tf.trrust, file = './cache/lnc2tf.trrust.rds')
system.time(lnc2tf.gtrd <- xyCor(fix.gtrd$a, fix.gtrd$b, cores=6))#9
saveRDS(lnc2tf.gtrd, file = './cache/lnc2tf.gtrd.rds')


lnc2tf.fimo.gss <- readRDS('./cache/lnc2tf.fimo.gss.rds')
lnc2tf.fimo.tss <- readRDS('./cache/lnc2tf.fimo.tss.rds')
lnc2tf.trrust <- readRDS('./cache/lnc2tf.trrust.rds')
lnc2tf.gtrd <- readRDS('./cache/lnc2tf.gtrd.rds')

lncTP.fimo.gss <- new('XY2Z', raw=lnc2tf.fimo.gss, x=fix.fimo.gss$a, y=fix.fimo.gss$b)
lncTP.fimo.tss <- new('XY2Z', raw=lnc2tf.fimo.tss, x=fix.fimo.gss$a, y=fix.fimo.tss$b)
lncTP.trrust <- new('XY2Z', raw=lnc2tf.trrust, x=fix.trrust$a, y=fix.trrust$b)
lncTP.gtrd <- new('XY2Z', raw=lnc2tf.gtrd, x=fix.gtrd$a, y=fix.gtrd$b)

#filter(lncTP.gtrd@detail, b=='MYC', a%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
#########################################

getFix <- function(tf, s=0) {
  if (tf=='enricher') {
    tf.2.PCG <- get.tf.2.PCG.from.enricher()
  } else {
    tf.2.PCG <- get.tf.2.PCG.from.fimo()
  }
  lncRNA.2.PCG <- get.lncRNA.2.PCG(s)
  lncRNA.2.PCG.m <- x2yMatrix(lncRNA.2.PCG)
  tf.2.PCG.m <- x2yMatrix(tf.2.PCG)
  fix <- x2yMatrixAdjust(lncRNA.2.PCG.m, tf.2.PCG.m,y.names = pcg$GeneID)
  fix
}


lncTF.all <- function(s=0) {
  fix <- getFix(s)
  message(paste('Calculate Phi of',ncol(fix$a), 'lncRNA and', ncol(fix$b),'TF'))
  lncRNA.tf.list <- xyCor(fix$a, fix$b)
  write.csv(lncRNA.tf.list, file = paste0('./data/lncRNA.tf.list/lncRNA.tf.',tf,'.',s,'.csv'))
  lncRNA.tf.list
}

# system.time(lncRNA.tf.fimo.0.1 <- lncTF.all(0.1))
# system.time(lncRNA.tf.fimo.0.2 <- lncTF.all(0.2))
# system.time(lncRNA.tf.fimo.0.3 <- lncTF.all(0.3))
# system.time(lncRNA.tf.fimo.0.4 <- lncTF.all(0.4))
# system.time(lncRNA.tf.fimo.0.5 <- lncTF.all(0.5))
# system.time(lncRNA.tf.fimo.0.6 <- lncTF.all(0.6))
# system.time(lncRNA.tf.fimo.0.7 <- lncTF.all(0.7))
# system.time(lncRNA.tf.fimo.0.8 <- lncTF.all(0.8))
# system.time(lncRNA.tf.fimo.0.9 <- lncTF.all(0.9))

# save(lncRNA.tf.fimo.0.1,
#      lncRNA.tf.fimo.0.2,
#      lncRNA.tf.fimo.0.3,
#      lncRNA.tf.fimo.0.4,
#      lncRNA.tf.fimo.0.5,
#      lncRNA.tf.fimo.0.6,
#      lncRNA.tf.fimo.0.7,
#      lncRNA.tf.fimo.0.8,
#      lncRNA.tf.fimo.0.9,
#      file = './cache/lncRNA.tf.fimo.x.rda'
#      )

############################################################
load(file = './cache/lncRNA.tf.fimo.x.rda')

fix1 <- getFix(0.1)
fix2 <- getFix(0.2)
fix3 <- getFix(0.3)
fix4 <- getFix(0.4)
fix5 <- getFix(0.5)
fix6 <- getFix(0.6)
fix7 <- getFix(0.7)
fix8 <- getFix(0.8)
fix9 <- getFix(0.9)


lncTP.0.1 <- new('XY2Z', raw=lncRNA.tf.fimo.0.1, x=fix1$a, y=fix1$b)
lncTP.0.2 <- new('XY2Z', raw=lncRNA.tf.fimo.0.2, x=fix2$a, y=fix2$b)
lncTP.0.3 <- new('XY2Z', raw=lncRNA.tf.fimo.0.3, x=fix3$a, y=fix3$b)
lncTP.0.4 <- new('XY2Z', raw=lncRNA.tf.fimo.0.4, x=fix4$a, y=fix4$b)
lncTP.0.5 <- new('XY2Z', raw=lncRNA.tf.fimo.0.5, x=fix5$a, y=fix5$b)
lncTP.0.6 <- new('XY2Z', raw=lncRNA.tf.fimo.0.6, x=fix6$a, y=fix6$b)
lncTP.0.7 <- new('XY2Z', raw=lncRNA.tf.fimo.0.7, x=fix7$a, y=fix7$b)
lncTP.0.8 <- new('XY2Z', raw=lncRNA.tf.fimo.0.8, x=fix8$a, y=fix8$b)
lncTP.0.9 <- new('XY2Z', raw=lncRNA.tf.fimo.0.9, x=fix9$a, y=fix9$b)


save(lncTP.0.1, lncTP.0.2, lncTP.0.3, lncTP.0.4, lncTP.0.5, lncTP.0.6, lncTP.0.7, lncTP.0.8, lncTP.0.9, file = './cache/lncTP.0.x.rda')
load('./cache/lncTP.0.x.rda')


lncTP.0.1
lncTP.0.2
lncTP.0.3
lncTP.0.4
lncTP.0.5
lncTP.0.6
lncTP.0.7
lncTP.0.8
lncTP.0.9

par(mfrow=c(3,3))
for(s in seq(0.1,0.9,0.1)) {
  lnctp <- get(paste0('lncTP.',s))
  cat('---------------------------',s,'----------------------------------------');cat('\n')
  print(lnctp)
  print(sort(lnctp@nodes$y))
  print(ncol(lnctp@detail))
  if (nrow(lnctp@detail) > 0) {
    barplot(sort(lnctp@detail$c11), main=s)
  }
}


write.csv(lncTP.0.3@nodes$y, './data/lnctp.tf.csv')
write.csv(lncTP.0.3@nodes$x, './data/lnctp.lncRNA.csv')
write.csv(lncTP.0.3@detail, './data/lnctp.detail.csv')
