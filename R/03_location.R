suppressPackageStartupMessages(source('./lib/globals.R'))
suppressPackageStartupMessages(source('./lib/helpers.R'))
suppressPackageStartupMessages(library(karyoploteR))
load_all('./package/x2y/')

cor.pairs.info <- readRDS('./cache/cor.pairs.info.rds')
biomart <- helper.get.biomart()
cor.lnc2all <- filter(cor.pairs.info, type1%in%config$lncRNA, type2%in%config$PCGs,FDR<0.05,abs(r)>=0.3)

###################################### plot
get.granges <- function(lnc) {
    loc <- biomart[match(unique(lnc$v2),biomart$Gene.stable.ID),]%>%
        select(chr=Chromosome.scaffold.name, start=Gene.start..bp., end=Gene.end..bp.)%>%
        mutate(chr=paste0('chr',chr))
    toGRanges(loc)
}


plot.location <- function(lnc, main='',file.name='',out.dir='./', save=TRUE) {
    colors <- ifelse(lnc$r>0,'red','green')
    grange <- get.granges(lnc)

    plot.type <- 1
    pp <- getDefaultPlotParams(plot.type = plot.type)
    pp$data1height=500
    pp$leftmargin <- 0.04

    if (save) {
      tiff(file.path(out.dir,paste0(ifelse(file.name=='',main, file.name),'.tif')),width=1200, height=800)
    }
    kp <- plotKaryotype(genome = 'hg38', plot.type=plot.type, plot.params = pp,cex=1, main=main)
    kpDataBackground(kp, color = "#FFFFFFAA")
    kpPlotDensity(kp, grange, col="#3388FF55", border="#3388FF")
    kpPlotRegions(kp, data=grange, col=colors)
    if (save) {
      dev.off()
    }
}

################################### plot.lnc.split
plot.lnc.split <- function(lnc.split, outdir='./reports/grange/all') {
  n.na <- n.nuclear <- n.cytoplasmic <- 0
  lapply(seq_along(lnc.split),function(i){
    lnc.name <- names(lnc.split)[[i]]
    rci <- helper.get.cell.localization()[lnc.name,]$rci
    if (is.na(rci)) {
      n.na <<- n.na + 1
      outdir <- file.path(outdir,'na')
    }else if (rci <0) {
      n.nuclear <<- n.nuclear + 1
      outdir <- file.path(outdir,'nuclear')
    } else {
      n.cytoplasmic <<- n.cytoplasmic + 1
      outdir <- file.path(outdir,'cytoplasmic')
    }
    if (!file.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }
    main <- paste(lnc.name, rci)
    plot.location(lnc.split[[i]], main=main, file.name = lnc.name,out.dir=outdir)
  })->tmp
  message(paste0('Nuclear: ', n.nuclear, "; Cytoplasmic: ", n.cytoplasmic, "; NA: ", n.na))
}



######################################### lnc.split & Filtered
lnc.split <- split(cor.lnc2all,as.vector(cor.lnc2all$v1))

load('./cache/lncTP.0.x.rda')
lncTP.pcgs <- getZByX(lncTP.0.3)
sapply(seq_along(lnc.split),function(i){
  lnc.name <- names(lnc.split)[[i]]
  rci <- helper.get.cell.localization()[lnc.name,]$rci
  n <- nrow(lnc.split[[i]])
  r <- lnc.split[[i]]$r
  ratio <- sum(r>0)/n
  tf.ratio <-length(intersect(lncTP.pcgs[[lnc.name]], lnc.split[[i]]$v2))/n

  c(lnc.name, rci, n, ratio,mean(r),sd(r),tf.ratio)
})->lnc.m
lnc.summary <- data.frame(name=lnc.m[1,],
                          rci=as.numeric(lnc.m[2,]),
                          n=as.numeric(lnc.m[3,]),
                          tf.ratio=as.numeric(lnc.m[7,]),
                          p.ratio=as.numeric(lnc.m[4,]),
                          mean=as.numeric(lnc.m[5,]),
                          sd=as.numeric(lnc.m[6,]),
                          stringsAsFactors = F)

par(mfrow=c(2,3))
plot(sort(lnc.summary$rci), ylab='NC RCI')
plot(sort(lnc.summary$n), ylab='PCG Number')
plot(sort(lnc.summary$p.ratio), ylab='Positive Corralation Ratio')
plot(sort(lnc.summary$mean), ylab='r Mean')
plot(sort(lnc.summary$sd), ylab='r SD')
par(mfrow=c(1,1))

n.cutoff <- quantile(lnc.summary$n, 0.75)
sd.cutoff <- quantile(lnc.summary$sd, 0.5,na.rm=TRUE)
pr.cutoff1 <- quantile(lnc.summary$p.ratio, 0.1)
pr.cutoff2 <- quantile(lnc.summary$p.ratio, 0.5)

lnc.filtered <- filter(lnc.summary, n > n.cutoff, p.ratio < pr.cutoff1|p.ratio >= pr.cutoff2)
lnc.filtered$name
lncTP.0.3@nodes$x
intersect(lnc.filtered$name, lncTP.0.3@nodes$x)

###########################################################################
p.ratio.tag <- ifelse(lnc.filtered$p.ratio > 0.5, 'positive', 'negative')
rci.tag<-c()
for (x in lnc.filtered$rci) {
  if (is.na(x)) {
    tag <- 'no'
  }
  else if (x < 0) {
    tag <- 'nuclear'
  }
  else {
    tag <- 'cytoplasmic'
  }
  rci.tag <- c(rci.tag, tag)
}
lnc.filtered.table <- table(rci.tag, p.ratio.tag)
lnc.filtered.table

###
par(mfrow=c(2,3))
barplot(sort(lnc.filtered$tf.ratio), main = 'Total')
barplot(sort(filter(lnc.filtered,rci <0)$tf.ratio), main = 'nuclear')
barplot(sort(filter(lnc.filtered,rci >0)$tf.ratio), main = 'cytoplasmic')
barplot(sort(filter(lnc.filtered,p.ratio > 0.5)$tf.ratio), main = 'positive')
barplot(sort(filter(lnc.filtered,p.ratio < 0.5)$tf.ratio), main = 'negative')
par(mfrow=c(1,1))

##
lnc.candidate <- filter(lnc.filtered,rci<0, tf.ratio>0)
saveRDS(lnc.candidate, file='./cache/lnc.candidate.location.rds') # <== This is Candidate
surv <- readRDS(file='./cache/candidate.surv.rds')
lnc.tfpcg <- getYZByX(lncTP.0.3)

###############################################################################
lnc.split.filtered <- lnc.split[lnc.filtered$name]
plot.lnc.split(lnc.split.filtered, outdir = "./reports/grange/filtered")
plot.lnc.split(helper.get.final(),.outdir='./reports/final/chr-location')
plot.location <- function(lnc, main='',file.name='',out.dir='./', save=TRUE) {
  colors <- ifelse(lnc$r>0,'red','green')
  grange <- get.granges(lnc)

  plot.type <- 1
  pp <- getDefaultPlotParams(plot.type = plot.type)
  pp$data1height=500
  pp$leftmargin <- 0.04

  if (save) {
    tiff(file.path(out.dir,paste0(ifelse(file.name=='',main, file.name),'.tif')),width=1200, height=800)
  }
  kp <- plotKaryotype(genome = 'hg38', plot.type=plot.type, plot.params = pp,cex=1, main=main)
  kpDataBackground(kp, color = "#FFFFFFAA")
  kpPlotDensity(kp, grange, col="#3388FF55", border="#3388FF")
  kpPlotRegions(kp, data=grange, col=colors)
  if (save) {
    dev.off()
  }
}

plot.location(lnc.split.filtered[[1]], save = FALSE)


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
####################################
gr1 <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(103,106),
  strand = "+",
  score = 5L, GC = 0.45)

gr2 <- GRanges(
  seqnames = c("chr1", "chr1"),
  ranges = IRanges(c(107,113), width = 3),
  strand = c("+", "-"),
  score = 3:4, GC = c(0.3,0.5)
)

grl = GRangesList("txA" = gr1, "txB" = gr2)




##################################3
?table(cor.lnc2all$scaffold2)

head(cor.lnc2all)
dim(cor.lnc2all)

as.vector(unique(cor.lnc2all$scaffold2))

lapply(lnc2all.splitsplit(cor.lnc2all,as.vector(cor.lnc2all$v1)),function(lnc){
    n <- nrow(lnc)
    n1 <- length(lnc$r[lnc$r>0])
    n2 <- length(lnc$r[lnc$r<0])
    n1_ratio <- length(lnc$r[lnc$r>0])/nrow(lnc)
    chr=c()
    for(i in c(1:23,'x')){
        gene.num <- nrow(filter(lnc, scaffold2==i))
        chr=c(chr, gene.num)
    }
    c(n, n1, n2, n1_ratio, chr)
})


head(cor.lnc2all)


#############################################


dim(biomart%>%filter(Gene.type%in%config$lncRNA)->lncRNAs)
lncRNAs$Gene.stable.ID%>%unique()->lncRNAs.unique

dim(lncRNAs)
length(lncRNAs.unique)
head(lncRNAs)

####
diff.lncRNA <- helper.get.lncRNA()
head(diff.lncRNA)
write.table(diff.lncRNA$GeneID, file='./reports/diff.lncRNA.646.list.csv', row.names=FALSE,col.names=FALSE,quote=FALSE)

cor.pairs$

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("karyoploteR")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
  ########################################################
   library(karyoploteR)

  kp <- plotKaryotype(genome = 'hg38', chromosomes=c("chr10", "chr12", "chr2"))
  kpAddBaseNumbers(kp)
  regions <- createRandomRegions(nregions=400, length.mean = 3e6, mask=NA,genome='hg38')
