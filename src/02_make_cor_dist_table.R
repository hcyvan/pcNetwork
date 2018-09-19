library('ProjectTemplate')
library(data.table)
load.project()

diff.lncrna.pcg <- diff.lncrna.pcg.551
data <- helper.getFpkm()$data

genes <- filter(data, GeneID%in%diff.lncrna.pcg$GeneID)
## Distance
genes.with.anno <- left_join(diff.lncrna.pcg, bioMart, by=c('GeneID'='Gene.stable.ID')) %>%
  select(colnames(diff.lncrna.pcg), 
         tss = Transcription.start.site..TSS.,
         geneStart = Transcript.start..bp., 
         geneEnd = Transcript.end..bp., 
         chromosome = Chromosome.scaffold.name)

## Cor pvalue
### pvalue
genes.fpkm <- as.matrix(genes[,-1])

cor.test.multi <- function(data) {
  i <- 0
  diff.pvalue <- apply(data,1,function(x){
    print(i<<-i+1)
    apply(data,1,function(y){
      cor.test(x,y)$p.value
    })
  })
}
#system.time(diff.pvalue.106 <- cor.test.multi(genes.fpkm))
#save(diff.pvalue.106, file='cache/diff.pvalue.106.rda')

getCorDistTable <- function(pvalue, genes) {
  message('1. Reshape p-value...')
  pvalue.pairs <- filter(reshape2::melt(pvalue), match(Var1, genes$GeneID) < match(Var2, genes$GeneID))
  message('2. Calculate and add cor...')
  corr <- cor(t(as.matrix(genes[,-1])))
  rownames(corr) <- dimnames(pvalue)[[1]]
  colnames(corr) <- dimnames(pvalue)[[2]]
  corr.pairs <- filter(reshape2::melt(corr), match(Var1, genes$GeneID) < match(Var2, genes$GeneID))
  corr.pairs <- left_join(corr.pairs, pvalue.pairs, by=c('Var1'='Var1', 'Var2'='Var2')) %>%
                select(Var1, Var2, r=value.x, p.value=value.y)
  message('3. Add gene type, distance and significant...')
  g1 <- genes.with.anno[match(pvalue.pairs$Var1,genes.with.anno$GeneID),]
  g2 <- genes.with.anno[match(pvalue.pairs$Var2,genes.with.anno$GeneID),]
  corr.pairs <- data.frame(corr.pairs,
                           type1=g1$GeneType,
                           type2=g2$GeneType,
                           dist=ifelse(g1$chromosome==g2$chromosome, g1$tss-g2$tss,Inf),
                           significant=pvalue.pairs$value<(0.05/nrow(pvalue.pairs)))
  corr.pairs
}

dimnames(diff.pvalue) <- list(genes$GeneID, genes$GeneID)
diff.cor.pairs <- getCorDistTable(diff.pvalue, genes)
save(diff.cor.pairs, file = './cache/diff.cor.pairs.rda')
