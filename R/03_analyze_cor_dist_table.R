load('cache/diff.cor.pairs.rda')
# load('cache/diff.cor.pairs.106.rda')

getCisGene <- function(pairs, b1, b2=0, significant=TRUE) {
  genes <- pairs %>%
          filter(abs(dist) < b1 & abs(dist) >= b2) %>%
          filter((type1%in%config$lncRNA & type2%in%config$PCGs) | (type1%in%config$PCGs & type2%in%config$lncRNA))
  if (significant) {
    genes <- filter(genes, significant)
  }
  genes <- genes %>%
          mutate(symbol1=helper.getGeneSymbol(Var1), symbol2=helper.getGeneSymbol(Var2)) %>%
          mutate_if(is.factor, as.character)
  genes
}


cis.2k.551 <- getCisGene(diff.cor.pairs, 2000) 
write.csv(cis.2k.551, 'reports/cis.2k.551.csv')
cis.63k.127k <- getCisGene(diff.cor.pairs, 127000, 63000)
write.csv(cis.63k.127k, file = 'reports/cis.63k.127k.csv')
cis.65535k.131071k <- getCisGene(diff.cor.pairs, 131071000, 65535000)
write.csv(cis.65535k.131071k, file = 'reports/cis.65535k.131071k.csv')


## --------------------------------------------- True Ratio -------------------------------------------
calculateSigRatio <- function(bound) {
  b1 <- bound[2]
  b2 <- bound[1]
  message('bound ', b2, ' to ', b1)
  genes <- getCisGene(diff.cor.pairs, b1, b2, FALSE)
  total <- table(genes$significant)
  f <- as.numeric(total['FALSE'])
  t <- as.numeric(total['TRUE'])
  t <- ifelse(is.na(t), 0, t)
  f <- ifelse(is.na(f), 0, f)
  c(b2, b1, t, f, t/(t+f))
}

a <- 2^seq(0,18)*1000 -1000
b <- 2^seq(1,19)*1000 -1000
cor.ratio.by.dist <- t(apply(cbind(a,b), 1, calculateSigRatio))
colnames(cor.ratio.by.dist) <- c('from', 'to', 'true', 'false', 'true.ratio')

# ============================
write.csv(getCisGene(diff.cor.pairs, Inf,0), file = './reports/diff.cor.pairs.csv')
