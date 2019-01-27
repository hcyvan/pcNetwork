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

fimo.gss.pairs <- readRDS('./cache/fimo.gss.set.rds')
fimo.tss.pairs <- readRDS('./cache/fimo.tss.set.rds')
trrust.pairs <- readRDS('./cache/trrust.set.rds')
gtrd.pairs <- readRDS('./cache/gtrd.set.rds')

cor.pairs.format <- function(s=0.3) {
    biomart <- pf.get.biomart()
    cor.pairs <- readRDS('./cache/cor.pearson.pairs.rds')
    t1 <- biomart[match(cor.pairs$v1, biomart$ensembl_gene_id), 'gene_biotype']
    t2 <- biomart[match(cor.pairs$v2, biomart$ensembl_gene_id), 'gene_biotype']
    data.frame(cor.pairs,t1=t1,t2=t2) %>%
        filter((t1%in%pv.lncRNA & t2%in%pv.pcg), FDR<0.05, abs(r) >= s) %>%
        select(lncRNA=v1, gene=v2)
}
cp.3 <- cor.pairs.format(s=0.3)
cp.4 <- cor.pairs.format(s=0.4)
cp.5 <- cor.pairs.format(s=0.5)
cp.8 <- cor.pairs.format(s=0.8)

dEnricher <- function(pair1, pair2, background=NULL, rds=NULL, cores=6, refresh=FALSE) {
    if (refresh || is.null(rds) || !file.exists(rds)) {
        message('Calulate XY2Z ...\n')
        m.adj <-getX2yMatrixAdjust(pair1, pair2, background)
        x2y <- xyCor(m.adj$a, m.adj$b, cores=6)
        xy2z <- new('XY2Z', raw=x2y, x=m.adj$a, y=m.adj$b)
        if (!is.null(rds)) {
            saveRDS(xy2z, file = rds)
        }
        xy2z
    } else {
        message(paste0('Load XY2Z from ', rds, '...\n'))
        readRDS(rds)
    }
}

fimo.gss.3 <- dEnricher(cp.3, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.3.xy2z.rds')
fimo.tss.3 <- dEnricher(cp.3, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.3.xy2z.rds')
trrust.3 <- dEnricher(cp.3, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.3.xy2z.rds')
gtrd.3 <- dEnricher(cp.3, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.3.xy2z.rds')

fimo.gss.4 <- dEnricher(cp.4, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.4.xy2z.rds')
fimo.tss.4 <- dEnricher(cp.4, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.4.xy2z.rds')
trrust.4 <- dEnricher(cp.4, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.4.xy2z.rds')
gtrd.4 <- dEnricher(cp.4, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.4.xy2z.rds')

fimo.gss.5 <- dEnricher(cp.5, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.5.xy2z.rds')
fimo.tss.5 <- dEnricher(cp.5, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.5.xy2z.rds')
trrust.5 <- dEnricher(cp.5, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.5.xy2z.rds')
gtrd.5 <- dEnricher(cp.5, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.5.xy2z.rds')

fimo.gss.8 <- dEnricher(cp.8, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.8.xy2z.rds')
fimo.tss.8 <- dEnricher(cp.8, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.8.xy2z.rds')
trrust.8 <- dEnricher(cp.8, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.8.xy2z.rds')
gtrd.8 <- dEnricher(cp.8, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.8.xy2z.rds')


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
