source('./lib/globals.R')
source('./lib/helpers.R')

diff.gene <- helper.get.lncRNA.PCG()
genes.fpkm <- helper.get.diff.fpkm()
genes.anno <- helper.get.diff.anno()

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

cor.pairs <- multicor(genes.fpkm, rds = './cache/cor.pearson.pairs.rds')

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
