source('./lib/globals.R')
source('./lib/helpers.R')

diff.gene <- helper.get.lncRNA.PCG()
data.fpkm <- helper.get.fpkm.count()
biomart <- helper.get.biomart()

genes.fpkm.f <- data.fpkm[rownames(data.fpkm)%in%diff.gene$GeneID,]
genes.fpkm <- data.fpkm[match(diff.gene$GeneID,rownames(data.fpkm)),]


genes.anno <- left_join(diff.gene, biomart, by=c('GeneID'='Gene.stable.ID')) %>%
  select(colnames(diff.gene), 
         tss = Transcription.start.site..TSS.,
         geneStart = Transcript.start..bp., 
         geneEnd = Transcript.end..bp., 
         chromosome = Chromosome.scaffold.name)

## Cor pvalue
### pvalue
cor.test.multi <- function(data) {
  data <- as.matrix(data)
  i <- 0
  diff.pvalue <- apply(data,1,function(x){
    print(i<<-i+1)
    apply(data,1,function(y){
      cor.test(x,y)$p.value
    })
  })
}
load('./cache/diff.qlf.2877.cor.pvalue.rda')
# system.time(cor.pvalue <- cor.test.multi(genes.fpkm))
# save(cor.pvalue, file='cache/diff.qlf.2877.cor.pvalue.rda')

getCorDistTable <- function(genes, anno, pvalue) {
  dimnames(pvalue) <- list(anno$GeneID, anno$GeneID)
  message('1. Reshape p-value...')
  pvalue.pairs <- filter(reshape2::melt(pvalue), match(Var1, anno$GeneID) < match(Var2, anno$GeneID))
  message('2. Calculate and add cor...')
  corr <- cor(t(as.matrix(genes)))
  rownames(corr) <- dimnames(pvalue)[[1]]
  colnames(corr) <- dimnames(pvalue)[[2]]
  corr.pairs <- filter(reshape2::melt(corr), match(Var1, anno$GeneID) < match(Var2, anno$GeneID))
  corr.pairs <- left_join(corr.pairs, pvalue.pairs, by=c('Var1'='Var1', 'Var2'='Var2')) %>%
                select(Var1, Var2, r=value.x, p.value=value.y)
  message('3. Add gene type, distance and significant...')
  g1 <- anno[match(pvalue.pairs$Var1, anno$GeneID),]
  g2 <- anno[match(pvalue.pairs$Var2, anno$GeneID),]
  corr.pairs <- data.frame(corr.pairs,
                           type1=g1$GeneType,
                           type2=g2$GeneType,
                           dist=ifelse(g1$chromosome==g2$chromosome, g1$tss-g2$tss,Inf),
                           significant=pvalue.pairs$value<(0.05/nrow(pvalue.pairs)))
  corr.pairs
}

diff.cor.pairs <- getCorDistTable(genes.fpkm, genes.anno, cor.pvalue)
save(diff.cor.pairs, file = './cache/diff.qlf.2877.pairs.rda')


################################## lncRNA vs mRNA
corr <- cor(t(as.matrix(genes.fpkm)))
lncRNA_pcg_corr <- corr[genes.anno$GeneType%in%config$PCGs,genes.anno$GeneType%in%config$lncRNA]
heatmap(lncRNA_pcg_corr[sample(nrow(lncRNA_pcg_corr), 200),sample(ncol(lncRNA_pcg_corr), 50)])

