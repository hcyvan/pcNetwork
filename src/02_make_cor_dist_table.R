library('ProjectTemplate')
library(data.table)
load.project()

genes <- filter(fpkm.data, GeneID %in% diffLncRnaAndPcg[,'GeneID'])

## Distance
genes.with.anno <- left_join(diffLncRnaAndPcg, bioMart, by=c('GeneID'='Gene.stable.ID')) %>%
  select(colnames(diffRna), 
         tss = Transcription.start.site..TSS.,
         geneStart = Transcript.start..bp., 
         geneEnd = Transcript.end..bp., 
         chromosome = Chromosome.scaffold.name)

## Cor pvalue
### pvalue
genes.fpkm <- apply(as.matrix(genes)[,-1],2,as.numeric)
# i <- 0
# diff.pvalue <- apply(genes.fpkm,1,function(x){
#   print(i<<-i+1)
#   apply(genes.fpkm,1,function(y){
#     cor.diff.pvalue.pairs(x,y)$p.value
#   })
# })
# save(diff.pvalue, file='cache/diff.pvalue.rda')
dimnames(diff.pvalue) <- list(genes$GeneID, genes$GeneID)
diff.pvalue.pairs <- reshape2::melt(diff.pvalue)
diff.pvalue.pairs <- dplyr::filter(diff.pvalue.pairs,match(Var1,genes$GeneID)<match(Var2,genes$GeneID))
### Cor
diff.cor <- cor(t(genes.fpkm))
dimnames(diff.cor) <- list(genes$GeneID, genes$GeneID)
diff.cor.pairs <- reshape2::melt(diff.cor)
diff.cor.pairs <- filter(diff.cor.pairs,match(Var1,genes$GeneID)<match(Var2,genes$GeneID))
diff.cor.pairs <- left_join(diff.cor.pairs, diff.pvalue.pairs, by=c('Var1'='Var1', 'Var2'='Var2')) %>%
                  select(Var1, Var2, r=value.x, p.value=value.y)
### Add dist and significant
g1 <- genes.with.anno[match(diff.pvalue.pairs$Var1,genes.with.anno$GeneID),]
g2 <- genes.with.anno[match(diff.pvalue.pairs$Var2,genes.with.anno$GeneID),]
diff.cor.pairs <- data.table(diff.cor.pairs,
                             type1=genes.with.anno[match(diff.cor.pairs$Var1, genes.with.anno$GeneID), 'GeneType'],
                             type2=genes.with.anno[match(diff.cor.pairs$Var2, genes.with.anno$GeneID), 'GeneType'],
                             dist=ifelse(g1$chromosome==g2$chromosome, g1$tss-g2$tss,Inf),
                             significant=diff.pvalue.pairs$value<(0.05/nrow(diff.pvalue.pairs)))

save(diff.cor.pairs, file = './cache/diff.cor.pairs.rda')


