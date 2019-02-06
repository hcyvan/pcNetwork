library(pcProfile)
library(data.table)
library(ggplot2)
library(cowplot)
library(VennDiagram)

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
  #' print(table(pf.ensembl2biotype(unique(tf2gene$gene))))
}

statistic(tf2gene.jasper)
statistic(tf2gene.gtrd)
statistic(tf2gene.trrust)

p1 <- ggplot(tf2gene.jasper[,.(.N), by=(tf)], aes(x=N)) +
                geom_histogram(binwidth = 100) +
                labs(x='Target Gene numbers', y='TF numbers') +
                theme(axis.title = element_text(size = rel(.8)))
info1 <- ggplot_build(p1)
p1 <- p1 + annotate("text",
                    x = sum(info1$layout$panel_params[[1]]$x.range)/2,
                    y = sum(info1$layout$panel_params[[1]]$y.range)/2,
                    label = "tf2gene.jasper",size=rel(5))
              
p2 <- ggplot(tf2gene.gtrd[,.(.N),by=(tf)], aes(x=N)) +
                geom_histogram(binwidth = 100) +
                labs(x='Target Gene numbers', y='TF numbers') +
                theme(axis.title = element_text(size = rel(.8)))
info2 <- ggplot_build(p2)
p2 <- p2 + annotate("text",
                    x = sum(info2$layout$panel_params[[1]]$x.range)/2,
                    y = sum(info2$layout$panel_params[[1]]$y.range)/2,
                    label = "tf2gene.gtrd",size=rel(5))

p3 <- ggplot(tf2gene.trrust[,.(.N),by=(tf)], aes(x=N)) +
                geom_histogram(binwidth = 1) +
                labs(x='Target Gene numbers', y='TF numbers') +
                theme(axis.title = element_text(size = rel(.8)))
info3 <- ggplot_build(p3)
p3 <- p3 + annotate("text",
                    x = sum(info3$layout$panel_params[[1]]$x.range)/2,
                    y = sum(info3$layout$panel_params[[1]]$y.range)/2,
                    label = "tf2gene.trrust",size=rel(5))


T<-venn.diagram(list(JASPER=tf2gene.jasper$tf,
                     GTRD=tf2gene.gtrd$tf,
                     TRRUST=tf2gene.trrust$tf),
                filename = NULL,
                col=c('red','green','blue'),
                lty=1,
                lwd=3,
                cat.dist=0.1,
                margin=0.1,
                reverse=TRUE)

# png("./reports/thesis//tf2gene.png", height = 800, width = 800)
win.metafile(filename="./reports/thesis//tf2gene.emf",width=7,height=7)
plot_grid(p1, p2, p3, grobTree(T), labels = c('A','B','C','D'), label_size = 20)
dev.off()



