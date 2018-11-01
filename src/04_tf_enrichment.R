library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

source('./lib/globals.R')
source('./lib/helpers.R')

modules <- read.delim('./data/diff.qlf.2877.wgcna.color.csv', sep = ',',header = TRUE, stringsAsFactors = FALSE)

bg1 <- read.delim('./data/ENCODE_TF_ChIP-seq_2015.handle.txt', sep = ',',header = TRUE, stringsAsFactors = FALSE)
bg2 <- read.delim('./data/ChEA_2016.handle.txt', sep = '\t',header = TRUE, stringsAsFactors = FALSE)

bg1 <- bg1[,c('tf','symbol')]
bg3 <- rbind(bg1,bg2)

bg1 <- distinct(bg1)
bg2 <- distinct(bg2)
bg3 <- distinct(bg3)

barplot(sort(table(bg$tf)))
barplot(sort(table(bg$symbol)))


## ------------------------------- over-representation analysis ----------------------------------
### ---------------- for each color
colorEnricher <- function(c, bg, pvalueCutoff=0.05) {
  tf2gene <- bg[,c('tf','symbol')]
  m.c <- dplyr::filter(modules, color==c, type%in%config$PCGs)
  geneList <- bitr(m.c$id, fromType='ENSEMBL',toType=c('ENTREZID', 'SYMBOL'), OrgDb = org.Hs.eg.db)
  print(dim(geneList))
  x <- enricher(geneList$SYMBOL, TERM2GENE = tf2gene, pvalueCutoff = pvalueCutoff, minGSSize=NULL, maxGSSize=NULL)
  x
}

colorEnricherEach <- function(bg) {
  enrichers <- list(detail=list(), tfs=list())
  for(c in unique(modules$color)){
    print(c)
    result <- colorEnricher(c,bg=bg)
    enrichers$tfs[[c]] <- as.data.frame(result)$ID
    enrichers$detail[[c]] <- result
  }
  enrichers
}

enricher.bg1 <- colorEnricherEach(bg1)
enricher.bg2 <- colorEnricherEach(bg2)
enricher.bg3 <- colorEnricherEach(bg3)

### ---------------------for all
colorEnricherAll <- function(bg) {
  tf2gene <- bg[,c('tf','symbol')]
  m.c <- dplyr::filter(modules, type%in%config$PCGs)
  geneList <- bitr(m.c$id, fromType='ENSEMBL',toType=c('ENTREZID', 'SYMBOL'), OrgDb = org.Hs.eg.db)
  enricher.all <- enricher(geneList$SYMBOL, TERM2GENE = tf2gene, pvalueCutoff = 0.05, minGSSize=NULL, maxGSSize=NULL)
  list(detail=enricher.all, tfs=as.data.frame(enricher.all)$ID)
}
enricher.all.bg1 <- colorEnricherAll(bg1)
enricher.all.bg2 <- colorEnricherAll(bg2)
enricher.all.bg3 <- colorEnricherAll(bg3)


## ------------------------------- GSEA -----------------------------------------
### ---------------- for each color
colorGSEA <- function(c, bg, pvalueCutoff=0.05) {
  tf2gene <- bg[,c('tf','symbol')]
  m.c <- dplyr::filter(modules, color==c, type%in%config$PCGs) %>% arrange(desc(mg.cor))
  id.map <- bitr(m.c$id, fromType='ENSEMBL',toType=c('ENTREZID', 'SYMBOL'), OrgDb = org.Hs.eg.db)
  m.c <- m.c %>% inner_join(id.map, by = c("id"="ENSEMBL"))
  geneList<-m.c$mg.cor
  names(geneList)<-m.c$SYMBOL

  x <- GSEA(geneList, TERM2GENE=tf2gene, minGSSize=1, maxGSSize=6000, pvalueCutoff = pvalueCutoff)
  x
}

colorGSEAEach <- function(bg) {
  gseaers <- list(detail=list(), tfs=list())
  for(c in unique(modules$color)){
    print(c)
    result <- colorGSEA(c,bg)
    gseaers$detail[[c]] <- result
    gseaers$tfs[[c]] <- as.data.frame(result)$ID
  }
  gseaers
}

gseaer.bg1 <- colorGSEAEach(bg1)
gseaer.bg2 <- colorGSEAEach(bg2)
gseaer.bg3 <- colorGSEAEach(bg3)


### -------------------- for all
colorGSEAAll <- function(bg) {
  m.c <- dplyr::filter(modules, type%in%config$PCGs) %>% arrange(desc(mg.cor))
  id.map <- bitr(m.c$id, fromType='ENSEMBL',toType=c('ENTREZID', 'SYMBOL'), OrgDb = org.Hs.eg.db)
  m.c <- m.c %>% inner_join(id.map, by = c("id"="ENSEMBL"))
  geneList<-m.c$mg.cor
  names(geneList)<-m.c$SYMBOL
  
  gseaer.all <- GSEA(geneList, TERM2GENE=tf2gene, minGSSize=1, maxGSSize=6000, pvalueCutoff = 0.05)
  list(detail=gseaer.all, tfs=as.data.frame(gseaer.all)$ID)
}
gseaer.all.bg1 <- colorGSEAAll(bg1)
gseaer.all.bg2 <- colorGSEAAll(bg2)
gseaer.all.bg3 <- colorGSEAAll(bg3)














