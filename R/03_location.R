source('./lib/globals.R')
source('./lib/helpers.R')
library(karyoploteR)

cor.pairs.info <- readRDS('./cache/cor.pairs.info.rds')
biomart <- helper.get.biomart()



cor.lnc2all <- cor.pairs.info %>% filter(abs(r)>0.3,FDR<0.05) %>% filter(type1%in%config$lncRNA)

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

plot.location(cor.lnc2all, save = FALSE)

lnc.split <- split(cor.lnc2all,as.vector(cor.lnc2all$v1))
tmp<-lapply(seq_along(lnc.split),function(i){
    lnc.name <- names(lnc.split)[[i]]
    rci <- helper.get.cell.localization()[lnc.name,]$rci
    outdir <- './reports/grange/'
    if (is.na(rci)) {
      outdir <- paste0(outdir,'na')
    }else if (rci <0) {
      outdir <- paste0(outdir,'nuclear')
    } else {
      outdir <- paste0(outdir,'cytoplasmic')
    }
    main <- paste(lnc.name, rci)
    print(i)
    plot.location(lnc.split[[i]], main=main, file.name = lnc.name,out.dir=outdir)
})


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
  
  
