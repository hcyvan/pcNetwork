library(dplyr)
library(gplots)
library(ggplot2)
library(illuminaHumanv3.db)
library(VennDiagram)

source('./R/lib.R')

ls("package:illuminaHumanv3.db")
xx <- as.list(illuminaHumanv3SYMBOL[mappedkeys(illuminaHumanv3SYMBOL)])
prob.symbol <- as.data.frame(cbind(prob=names(xx),symbol=unlist(xx)),stringsAsFactors=FALSE)


c0 <- read.csv('./data/AR/GSE21245/beadchip_0.txt', comment.char = '#', sep = '\t')
c20min <- read.csv('./data/AR/GSE21245/beadchip_20min.txt', comment.char = '#', sep = '\t')
c40min <- read.csv('./data/AR/GSE21245/beadchip_40min.txt', comment.char = '#', sep = '\t')
c1h <- read.csv('./data/AR/GSE21245/beadchip_1h.txt', comment.char = '#', sep = '\t')
c2h <- read.csv('./data/AR/GSE21245/beadchip_2h.txt', comment.char = '#', sep = '\t')
c4h <- read.csv('./data/AR/GSE21245/beadchip_4h.txt', comment.char = '#', sep = '\t')
c8h <- read.csv('./data/AR/GSE21245/beadchip_8h.txt', comment.char = '#', sep = '\t')
c16h <- read.csv('./data/AR/GSE21245/beadchip_16h.txt', comment.char = '#', sep = '\t')
c24h <- read.csv('./data/AR/GSE21245/beadchip_24h.txt', comment.char = '#', sep = '\t')
c48h <- read.csv('./data/AR/GSE21245/beadchip_48h.txt', comment.char = '#', sep = '\t')


prob<-as.vector(c0$ID_REF)
value.raw <- data.frame(prob=prob,
                    t0=c0$VALUE,
                    t20min=c20min$VALUE,
                    t40min=c40min$VALUE,
                    t1h=c1h$VALUE,
                    t2h=c2h$VALUE,
                    t4h=c4h$VALUE,
                    t8h=c8h$VALUE,
                    t16h=c16h$VALUE,
                    t24h=c24h$VALUE,
                    t48h=c48h$VALUE,
                    stringsAsFactors = FALSE)

detection.pval.raw <- data.frame(prob=prob,
                            t0=c0$Detection.Pval,
                            t20min=c20min$Detection.Pval,
                            t40min=c40min$Detection.Pval,
                            t1h=c1h$Detection.Pval,
                            t2h=c2h$Detection.Pval,
                            t4h=c4h$Detection.Pval,
                            t8h=c8h$Detection.Pval,
                            t16h=c16h$Detection.Pval,
                            t24h=c24h$Detection.Pval,
                            t48h=c48h$Detection.Pval,
                            stringsAsFactors = FALSE)

diffscore.raw <- data.frame(prob=prob,
                        t0=c0$DiffScore,
                        t20min=c20min$DiffScore,
                        t40min=c40min$DiffScore,
                        t1h=c1h$DiffScore,
                        t2h=c2h$DiffScore,
                        t4h=c4h$DiffScore,
                        t8h=c8h$DiffScore,
                        t16h=c16h$DiffScore,
                        t24h=c24h$DiffScore,
                        t48h=c48h$DiffScore,
                        stringsAsFactors = FALSE)


beadchip <- list(
  value=distinct(inner_join(prob.symbol, value.raw, by='prob'), symbol, .keep_all = TRUE),
  detection.pval=distinct(inner_join(prob.symbol, detection.pval.raw, by='prob'), symbol, .keep_all = TRUE),
  diffscore=distinct(inner_join(prob.symbol, diffscore.raw, by='prob'), symbol, .keep_all = TRUE)
)


big.diffscore <- abs(beadchip$diffscore[,c(-1,-2)])>13
valid.detection.pval <- beadchip$detection.pval[,c(-1,-2)]<0.05
diff.gene <- (big.diffscore + valid.detection.pval)==2

################## Diff Gene Num
# Total number of different gene
sum(apply(diff.gene, 1, sum) >=1) # 5317 vs 5336

# diff.gene each time point
dge <- lapply(as.data.frame(diff.gene), function(x){
  beadchip$value$symbol[x]
})
dge.all <- Reduce(union, dge)
dge.early <- Reduce(union, dge[c('t20min','t40min','t1h','t2h','t4h')])
dge.late <- Reduce(union, dge[c('t8h','t16h','t24h','t48h')])
dge.early.only <- setdiff(dge.early, dge.both)
dge.late.only <- setdiff(dge.late, dge.both)
dge.both <- intersect(dge.early,dge.late)

length(dge.early) # 4009
length(dge.late) # 3501
length(dge.both) # 2193

diff.gene.index <- apply(diff.gene,1, sum)>=1
beadchip.diff <- list(
  value=beadchip$value[diff.gene.index,],
  detection.pval=beadchip$detection.pval[diff.gene.index,],
  diffscore=beadchip$diffscore[diff.gene.index,]
)

##################### DiffScore Z-Score Heatmap
labels <- c('20min', '40min', '1h', '2h', '4h', '8h', '16h', '24h', '48h')
png(filename=paste0('./reports/thesis/androgen_response_gene_A.png'),width=600,height=800)
# win.metafile(paste0('./reports/thesis/androgen_response_gene_A.emf'))
# pdf(filename=paste0('./reports/thesis/androgen_response_gene_A.pdf'),width=1024,height=728)

m <- as.matrix(beadchip.diff$diffscore[,4:12])
colnames(m) <- labels
par(cex.main=50)
heatmap.2(m,
          scale = 'row',
          Colv=FALSE,
          col=bluered,
          trace = 'none',
          cexCol = 2,
          srtCol = 45,
          key.title = NA,
          key.xlab = 'Z-Score of DiffScore',
          key.ylab = NA,
          key.par = list(),
          lhei = c(1,4),
          colsep=5)
dev.off()

##################### Time course distribution
diff.count <- sapply(dge, length)
df <- data.frame(point=factor(labels, levels=labels), num=diff.count[-1])
rownames(df) <- NULL
p1 <- ggplot(df, aes(x=point, y=num, group=1)) +
  geom_line(size=1) +
  geom_point() +
  labs(x=NULL,y='Differential gene number') +
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
p2 <- venn.diagram(list(EARLY=dge.early,
                        LATE=dge.late),
                   filename = NULL,
                   lty=1,
                   lwd=5,
                   col=c('red','blue'),
                   cex=2,
                   
                   fill=NA,
                   cat.cex=2,
                   cat.dist=0.05,
                   margin=0.1,
                   rotation.degree=45,
                   reverse=TRUE)


dev.off()
# png(filename=paste0('./reports/thesis/androgen_response_gene_B.png'),width=1200,height=600)
win.metafile(filename=paste0('./reports/thesis/androgen_response_gene_B.emf'),width=14,height=7)
plot_grid(p1,  grobTree(p2), labels = c('A','B'), label_size = 30)
dev.off()

#####################################################################################
library(pcProfile)

get.ar.target <- function(tf2gene) {
  ar.target <- filter(tf2gene,tf=='AR')$gene
  ar.target <- data.frame(gene=pf.ensembl2symbol(ar.target), type=pf.ensembl2biotype(ar.target), stringsAsFactors = FALSE)
  genes <- unique(ar.target$gene)
  genes[-which(genes=='')]
}


ar.target.gtrd <- get.ar.target(tf2gene.gtrd)
ar.target.jasper <- get.ar.target(tf2gene.jasper)
ar.target.trrust <- get.ar.target(tf2gene.trrust)

length(ar.target.gtrd)
length(ar.target.jasper)
length(ar.target.trrust)




length(dge.early.only)
length(dge.late.only)

length(intersect(ar.target.gtrd, dge.early.only))
length(intersect(ar.target.gtrd, dge.late.only))
length(intersect(ar.target.gtrd, dge.all))

length(intersect(ar.target.jasper, dge.early.only))
length(intersect(ar.target.jasper, dge.late.only))
length(intersect(ar.target.jasper, dge.all))

