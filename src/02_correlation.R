source('./lib/globals.R')
source('./lib/helpers.R')

diff.gene <- helper.get.lncRNA.PCG()
data.fpkm <- helper.get.fpkm.count()
biomart <- helper.get.biomart()

genes.fpkm <- data.fpkm[match(diff.gene$GeneID,rownames(data.fpkm)),]


genes.anno <- left_join(diff.gene, biomart, by=c('GeneID'='Gene.stable.ID')) %>%
  select(colnames(diff.gene), 
         tss = Transcription.start.site..TSS.,
         geneStart = Transcript.start..bp., 
         geneEnd = Transcript.end..bp., 
         chromosome = Chromosome.scaffold.name)


multicor <- function(data, method= c('pearson', 'kendall', 'spearman'), rds=NA, rewrite=FALSE, verbose=TRUE) {
  # Data structure like this:
  #
  #       sample1  sample2  sample3
  # gene1   xx      xx      xx
  # gene2   xx      xx      xx
  # gene3   xx      xx      xx
  #
  # Return:
  #         v1, v2, r, p.value, FDR
  
  if(!is.na(rds)) {
    if(file.exists(rds) && !rewrite) {
      message(paste('Load data from', rds))
      return(readRDS(rds))
    }
  }
  
  method <- match.arg(method)
  data <- as.matrix(data)
  r <- matrix(nrow=nrow(data),ncol=nrow(data))
  if (verbose) {
    message(paste('Calculating:', nrow(data),'rounds needed!'))
  }
  i <- 0
  p <- apply(data,1,function(x){
    i<<-i+1
    t0 <- Sys.time()
    if (verbose) {
      cat(paste0(i, ' '))
    }
    j<-0
    apply(data,1,function(y){
      j<<-j+1
      if (i > j) {
        test <- cor.test(x,y,method = method)
        r[j,i] <<- test$estimate
        test$p.value
      } else {
        NaN
      }
    })
  })
  if (verbose) {
    cat('\n')
  }
  dimnames(r) <- dimnames(p)
  list(r,p)
  r.melt <- reshape2::melt(r) %>% filter(!is.na(value))
  p.melt <- reshape2::melt(p) %>% filter(!is.na(value))
  ret <- data.frame(v1=r.melt$Var1, v2=r.melt$Var2, r=r.melt$value, p.value=p.melt$value, FDR=p.adjust(p.melt$value, method = 'BH'))
  if (!is.na(rds)) {
    saveRDS(ret, rds)
  }
  ret
}

cor.pairs <- multicor(genes.fpkm, rds = './cache/test.rds')



annoPairs <- function(cor.pairs, annotation) {
  g1 <- annotation[match(cor.pairs$v1, annotation$GeneID),]
  g2 <- annotation[match(cor.pairs$v2, annotation$GeneID),]
  data.frame(v1=cor.pairs$v1,
             v2=cor.pairs$v2,
             type1=g1$GeneType,
             type2=g2$GeneType,
             r=cor.pairs$r,
             p.value=cor.pairs$p.value,
             FDR=cor.pairs$FDR,
             scaffold1=g1$chromosome,
             scaffold2=g2$chromosome,
             distance=ifelse(g1$chromosome==g2$chromosome, g1$tss-g2$tss,Inf))
}

system.time(cor.pairs.info <- annoPairs(cor.pairs, genes.anno))
saveRDS(cor.pairs.info, file = './cache/cor.pairs.info.rds')



############################################################### Old Version
## Cor pvalue
### pvalue
cor.test.multi <- function(data) {
  data <- as.matrix(data)
  i <- 0
  diff.pvalue <- apply(data,1,function(x){
    i<<-i+1
    cat(paste0(i, ' '))
    apply(data,1,function(y){
      cor.test(x,y)$p.value
    })
  })
}

system.time(cor.pvalue <- cor.test.multi(genes.fpkm))
# saveRDS(cor.pvalue, './cache/diff.qlf.2877.cor.pvalue.rds')
cor.pvalue <- readRDS('./cache/diff.qlf.2877.cor.pvalue.rds')

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

# diff.cor.pairs <- getCorDistTable(genes.fpkm, genes.anno, cor.pvalue)
# saveRDS(diff.cor.pairs, file = './cache/diff.qlf.2877.pairs.rds')
diff.cor.pairs<-readRDS('./cache/diff.qlf.2877.pairs.rds')

################################## lncRNA vs mRNA
corr <- cor(t(as.matrix(genes.fpkm)))
lncRNA_pcg_corr <- corr[genes.anno$GeneType%in%config$PCGs,genes.anno$GeneType%in%config$lncRNA]
heatmap(lncRNA_pcg_corr[sample(nrow(lncRNA_pcg_corr), 200),sample(ncol(lncRNA_pcg_corr), 50)])

