library(stringr)
library(dplyr)
source('./lib/lib.R')

helper.get.data.count <- function() {
  data.count <- read.delim('./data/data.count.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE, row.names = 'GeneID')
  rownames(data.count) <- str_split_fixed(rownames(data.count), '\\.',2)[,1]
  data.count
}

helper.get.fpkm <- function(refresh=FALSE) {
  rds.path <- './cache/data.fpkm.rds'
  csv.path <- './data/data.fpkm.csv'
  if (file.exists(rds.path)&& !refresh) {
    data.fpkm <- readRDS(rds.path)
  } else {
    data.fpkm <- read.delim(csv.path, sep = ',', header = TRUE, stringsAsFactors = FALSE, row.names = 'GeneID')
    rownames(data.fpkm) <- str_split_fixed(rownames(data.fpkm), '\\.',2)[,1]
    saveRDS(data.fpkm, file = rds.path)
  }
  data.fpkm
}

helper.get.biomart <- function(refresh=FALSE) {
  rds.path <- './cache/mart_export.rds'
  csv.path <- './data/mart_export.csv'
  if (file.exists(rds.path)&& !refresh) {
    biomart <- readRDS(rds.path)
  } else {
    biomart <- read.delim(csv.path, sep = ',', header = TRUE, stringsAsFactors = FALSE)
    biomart <- distinct(biomart, Gene.stable.ID, .keep_all = TRUE)
    saveRDS(biomart, file = rds.path)
  }
  biomart
}

helper.get.diff.fpkm <- function(refresh=FALSE) {
  data <- helper.get.fpkm(refresh=refresh)
  diff.gene <- helper.get.lncRNA.PCG()
  genes.fpkm <- data[match(diff.gene$GeneID,rownames(data)),]
}

helper.get.diff.anno <- function(refresh=FALSE) {
  biomart <- helper.get.biomart(refresh=refresh)
  diff.gene <- helper.get.lncRNA.PCG()
  genes.anno <- left_join(diff.gene, biomart, by=c('GeneID'='Gene.stable.ID')) %>%
    select(colnames(diff.gene),
           tss = Transcription.start.site..TSS.,
           geneStart = Transcript.start..bp., 
           geneEnd = Transcript.end..bp., 
           chromosome = Chromosome.scaffold.name)
  print(head(genes.anno))
  genes.anno
}

helper.get.diff.gene <- function() {
  read.delim('./data/diff.qlf.csv', sep = ',', stringsAsFactors = FALSE)
}
helper.get.lncRNA.PCG <- function() {
  lncRNA <- read.delim('./data/diff.qlf.lncRNA.csv', sep = ',', stringsAsFactors = FALSE)
  pcg <- read.delim('./data/diff.qlf.pcg.csv', sep = ',', stringsAsFactors = FALSE)
  rbind(lncRNA, pcg)
}
helper.get.PCG <- function() {
  read.delim('./data/diff.qlf.pcg.csv', sep = ',', stringsAsFactors = FALSE)
}
helper.get.lncRNA <- function() {
  read.delim('./data/diff.qlf.lncRNA.csv', sep = ',', stringsAsFactors = FALSE)
}

helper.get.cell.localization <- function(id) {
  env = globalenv()
  key='.cell.location'
  if (exists(key, envir = env)) {
    get(key,envir = env)
  } else {
    cell.location <- readRDS('./cache/cell.location.rds')
    assign(key, cell.location, envir = env)
    cell.location
  }
}

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

.getRawData <- function(rawdata, sample.fix=FALSE) {
  sample <- list()
  if (sample.fix) {
    sample$data <- data.frame(rawdata[,1], rawdata[,helper.fixSample()])
    sample$group <- c(rep('T', 54),rep('N', 52))
    sample
  } else {
    sample$data <- as.data.frame(rawdata)
    sample$group <- c(rep('T', 499),rep('N', 52))
    sample
  }
}

helper.getCount <- function(sample.fix=FALSE) {
  .getRawData(count.data, sample.fix)  
}

helper.getFpkm <- function(sample.fix=FALSE) {
  .getRawData(fpkm.data, sample.fix)
}

.getAB <- function(a, b, sample) {
  .findGene <- function(id, sample.data) {
    gene <- filter(sample.data, GeneID==id)
    if (is.na(gene$GeneID[1])) {
      mart.export %>% filter(HGNC.symbol==id) -> tmp
      id.new <- tmp$Gene.stable.ID[1]
      gene <- filter(sample.data, GeneID==id.new)
    }
    if (is.na(gene$GeneID[1])) {
      stop('Can parse Gene ID: ', id)
    }
    as.numeric(as.matrix(gene[,-1]))
  }
  a.gene <- .findGene(a, sample)
  b.gene <- .findGene(b, sample)
  data.frame(a=a.gene, b=b.gene)
}

helper.plotGeneCor <- function(a, b, sample.fix=FALSE) {
  sample <- helper.getFpkm(sample.fix = sample.fix)
  ab <- .getAB(a, b, sample$data)
  
  dat <- data.frame(a=ab$a, b=ab$b, g=sample$group)
  dat.lm <- lm(b~a, data=dat)
  p <- ggplot(dat, aes(x = a, y = b, colour = g)) + geom_point()
  p + geom_abline(intercept = coef(dat.lm)[1],slope = coef(dat.lm)[2])
}

helper.plotGeneBox <- function(a, b, sample.fix=FALSE) {
  sample <- helper.getFpkm(sample.fix = sample.fix)
  ab <- .getAB(a, b, sample$data)
  
  dat.a <- data.frame(gene=a, group=sample$group, fpkm=ab$a)
  dat.b <- data.frame(gene=b, group=sample$group, fpkm=ab$b)
  dat <- rbind(dat.a, dat.b)
  ggplot(dat, aes(x=gene, y=fpkm, fill=group)) + geom_boxplot()
}

helper.plotGenePointAndBox <- function(a, b, sample.fix=FALSE) {
  sample <- helper.getFpkm(sample.fix = sample.fix)
  ab <- .getAB(a, b, sample$data)
  
  dat1 <- data.frame(a=ab$a, b=ab$b, g=sample$group)
  dat1.lm <- lm(b~a, data=dat1)
  p1 <- ggplot(dat1, aes(x = a, y = b, colour = g)) + geom_point() + labs(x=a,y=b)
  p1 <- p1 + theme(axis.title.x =element_text(size=10), axis.title.y=element_text(size=10))
  p1 <- p1 + geom_abline(intercept = coef(dat1.lm)[1],slope = coef(dat1.lm)[2])

  dat.a <- data.frame(gene=a, group=sample$group, fpkm=ab$a)
  dat.b <- data.frame(gene=b, group=sample$group, fpkm=ab$b)
  dat <- rbind(dat.a, dat.b)
  p2 <- ggplot(dat, aes(x=gene, y=fpkm, fill=group)) + geom_boxplot()
  .multiplot(p1,p2, cols = 2)
}

