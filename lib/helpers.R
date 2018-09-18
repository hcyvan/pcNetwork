source('./lib/lib.R')

helper.getGeneSymbol <- function(gene.id) {
  as.data.frame(bioMart)[match(gene.id, bioMart$Gene.stable.ID), 'HGNC.symbol']
}

helper.fixSample <- function() {
  samples <- colnames(fpkm.data)[-1]
  samples.t <- samples[1:499]
  samples.n <- samples[500:551]
  samples.t.id <- str_split_fixed(samples.t, '\\.', 4)[,3]
  samples.n.id <- str_split_fixed(samples.n, '\\.', 4)[,3]
  samples.t.paired <- samples.t[samples.t.id%in%samples.n.id]
  samples.paired <- c(samples.t.paired, samples.n)
  samples.paired
}

helper.getFpkm <- function(sample.fix=FALSE) {
  sample <- list()
  if (sample.fix) {
    sample$data <- data.frame(fpkm.data[,1], fpkm.data[,helper.fixSample()])
    sample$group <- c(rep('T', 54),rep('N', 52))
    sample
  } else {
    sample$data <- as.data.frame(fpkm.data)
    sample$group <- c(rep('T', 499),rep('N', 52))
    sample
  }
}

helper.plotGeneCor <- function(a, b, sample.fix=FALSE) {
  sample <- helper.getFpkm(sample.fix = sample.fix)
  a.fpkm <- as.numeric(as.matrix(filter(sample$data, GeneID==a)[,-1]))
  b.fpkm <- as.numeric(as.matrix(filter(sample$data, GeneID==b)[,-1]))
  dat <- data.frame(a=a.fpkm, b=b.fpkm, g=sample$group)
  dat.lm <- lm(b~a, data=dat)
  p <- ggplot(dat, aes(x = a, y = b, colour = g)) + geom_point()
  p + geom_abline(intercept = coef(dat.lm)[1],slope = coef(dat.lm)[2])
}

helper.plotGeneBox <- function(a, b, sample.fix=FALSE) {
  sample <- helper.getFpkm(sample.fix = sample.fix)
  a.fpkm <- as.numeric(as.matrix(filter(sample$data, GeneID==a)[,-1]))
  b.fpkm <- as.numeric(as.matrix(filter(sample$data, GeneID==b)[,-1]))
  dat.a <- data.frame(gene=a, group=sample$group, fpkm=a.fpkm)
  dat.b <- data.frame(gene=b, group=sample$group, fpkm=b.fpkm)
  dat <- rbind(dat.a, dat.b)
  ggplot(dat, aes(x=gene, y=fpkm, fill=group)) + geom_boxplot()
}

helper.plotGenePointAndBox <- function(a, b, sample.fix=FALSE) {
  sample <- helper.getFpkm(sample.fix = sample.fix)
  a.fpkm <- as.numeric(as.matrix(filter(sample$data, GeneID==a)[,-1]))
  b.fpkm <- as.numeric(as.matrix(filter(sample$data, GeneID==b)[,-1]))
  
  dat1 <- data.frame(a=a.fpkm, b=b.fpkm, g=sample$group)
  dat1.lm <- lm(b~a, data=dat1)
  p1 <- ggplot(dat1, aes(x = a, y = b, colour = g)) + geom_point() + labs(x=a,y=b)
  p1 <- p1 + theme(axis.title.x =element_text(size=10), axis.title.y=element_text(size=10))
  p1 <- p1 + geom_abline(intercept = coef(dat1.lm)[1],slope = coef(dat1.lm)[2])

  dat.a <- data.frame(gene=a, group=sample$group, fpkm=a.fpkm)
  dat.b <- data.frame(gene=b, group=sample$group, fpkm=b.fpkm)
  dat <- rbind(dat.a, dat.b)
  p2 <- ggplot(dat, aes(x=gene, y=fpkm, fill=group)) + geom_boxplot()
  .multiplot(p1,p2, cols = 2)
}
