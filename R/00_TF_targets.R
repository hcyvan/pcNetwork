library(pcProfile)
library(data.table)
library(ggplot2)
library(cowplot)

source('./R/lib.R')

biomart <- data.table(pf.get.biomart())

tf2gene.jasper <- data.table(tf2gene.jasper)
tf2gene.gtrd <- data.table(tf2gene.gtrd)
tf2gene.trrust <- data.table(tf2gene.trrust)

statistic <- function(tf2gene) {
  genes <- length(unique(tf2gene$gene))
  tfs <- length(unique(tf2gene$tf))
  tf.gene <- dim(tf2gene)[1]
  print(paste('genes:',genes,';tfs:', tfs, ';tf.gene:', tf.gene))
  # print(table(pf.ensembl2biotype(unique(tf2gene$gene))))
  
}

statistic(tf2gene.jasper)
statistic(tf2gene.gtrd)
statistic(tf2gene.trrust)


fimo.tss.hist <- ggplot(tf2gene.jasper[,.(.N), by=(tf)], aes(x=N))+geom_histogram(binwidth = 200)
gtrd.hist <- ggplot(tf2gene.gtrd[,.(.N),by=(tf)], aes(x=N))+geom_histogram(binwidth = 100)
trrust.hist <- ggplot(tf2gene.trrust[,.(.N),by=(tf)], aes(x=N))+geom_histogram(binwidth = 1)

plot_grid(fimo.tss.hist, gtrd.hist, trrust.hist, labels = c('A','B','C','D'))
