library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

source('./lib/globals.R')
source('./lib/helpers.R')
diff.gene <- helper.get.lncRNA.PCG()
pcg <- diff.gene %>% filter(GeneType%in%config$PCGs) %>% dplyr::select(id=GeneID, fc=logFC) %>% arrange(desc(fc)) 


# modules <- read.delim('./data/diff.qlf.2877.wgcna.color.csv', sep = ',',header = TRUE, stringsAsFactors = FALSE)

bg1 <- read.delim('./data/enricher/ENCODE_TF_ChIP-seq_2015.handle.txt', sep = ',',header = TRUE, stringsAsFactors = FALSE)
bg2 <- read.delim('./data/enricher/ChEA_2016.handle.txt', sep = '\t',header = TRUE, stringsAsFactors = FALSE)
bg3 <- read.delim('./data/enricher/TRANSFAC_and_JASPAR_PWMs.handle.txt', sep = '\t',header = TRUE, stringsAsFactors = FALSE)

bg1 <- bg1[,c('tf','symbol')]
bg0 <- rbind(bg1,bg2,bg3)

bg1 <- distinct(bg1)
bg2 <- distinct(bg2)
bg3 <- distinct(bg3)
bg0 <- distinct(bg0)

dim(bg1)
dim(bg2)
dim(bg3)
dim(bg0)

length(unique(bg1$tf))
length(unique(bg2$tf))
length(unique(bg3$tf))
length(unique(bg0$tf))

barplot(sort(table(bg0$tf)))
barplot(sort(table(bg0$symbol)))

# ---------------------------------------------- filter bg
library(biomaRt)
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
searchAttributes(mart = ensembl, pattern = "ensembl")
biomart.symbol.biotype = getBM(attributes = c('ensembl_gene_id','hgnc_symbol', 'gene_biotype'), mart = ensembl)
save(biomart.symbol.biotype, file = './cache/biomart.symbol.biotype.rda')
filter(biomart.symbol.biotype, gene_biotype%in%config$PCGs) -> biomart.pcgs

bg0.new <- filter(bg0, symbol%in%biomart.pcgs$hgnc_symbol)
bg1.new <- filter(bg1, symbol%in%biomart.pcgs$hgnc_symbol)
bg2.new <- filter(bg2, symbol%in%biomart.pcgs$hgnc_symbol)
bg3.new <- filter(bg3, symbol%in%biomart.pcgs$hgnc_symbol)

# -------------------------------------------

###################################################### Diff Gene Enricher
colorEnricherAll <- function(bg, pvalueCutoff=0.05) {
  tf2gene <- bg[,c('tf','symbol')]
  geneList <- inner_join(pcg, biomart.pcgs, by=c('id'='ensembl_gene_id'))$hgnc_symbol
  enricher.all <- enricher(geneList, TERM2GENE = tf2gene, pvalueCutoff = pvalueCutoff, minGSSize=NULL, maxGSSize=NULL)
  list(detail=enricher.all, tfs=as.data.frame(enricher.all)$ID)
}
enricher.all.bg1.001 <- colorEnricherAll(bg1,0.01)
enricher.all.bg2.001 <- colorEnricherAll(bg2,0.01)
enricher.all.bg3.001 <- colorEnricherAll(bg3,0.01)
enricher.all.bg0.001 <- colorEnricherAll(bg0,0.01)
enricher.all.bg0.new.001 <- colorEnricherAll(bg0.new,0.01)


enricher.all.bg1 <- colorEnricherAll(bg1)
enricher.all.bg2 <- colorEnricherAll(bg2)
enricher.all.bg3 <- colorEnricherAll(bg3)
enricher.all.bg0 <- colorEnricherAll(bg0)
enricher.all.bg0.new <- colorEnricherAll(bg0.new) ### Use this

write.csv(enricher.all.bg0.new, file = './data/enricher.all.bg0.new.csv', row.names = FALSE)

####################################################### WGCNA Moudle
## ------------------------------- over-representation analysis ----------------------------------
### ---------------- for each color
colorEnricher <- function(c, bg, pvalueCutoff) {
  tf2gene <- bg[,c('tf','symbol')]
  m.c <- dplyr::filter(modules, color==c, type%in%config$PCGs)
  geneList <- bitr(m.c$id, fromType='ENSEMBL',toType=c('ENTREZID', 'SYMBOL'), OrgDb = org.Hs.eg.db)
  print(dim(geneList))
  x <- enricher(geneList$SYMBOL, TERM2GENE = tf2gene, pvalueCutoff = pvalueCutoff, minGSSize=NULL, maxGSSize=NULL)
  x
}

colorEnricherEach <- function(bg, pvalueCutoff=0.05) {
  enrichers <- list(detail=list(), tfs=list())
  for(c in unique(modules$color)){
    print(c)
    result <- colorEnricher(c,bg=bg, pvalueCutoff = pvalueCutoff)
    enrichers$tfs[[c]] <- as.data.frame(result)$ID
    enrichers$detail[[c]] <- result
  }
  enrichers
}

enricher.bg1.001 <- colorEnricherEach(bg1,0.01)
enricher.bg2.001 <- colorEnricherEach(bg2,0.01)
enricher.bg3.001 <- colorEnricherEach(bg3,0.01)

enricher.bg1 <- colorEnricherEach(bg1)
enricher.bg2 <- colorEnricherEach(bg2)
enricher.bg3 <- colorEnricherEach(bg3)
### ---------------------for all
colorEnricherAll <- function(bg, pvalueCutoff=0.05) {
  tf2gene <- bg[,c('tf','symbol')]
  m.c <- dplyr::filter(modules, type%in%config$PCGs)
  geneList <- bitr(m.c$id, fromType='ENSEMBL',toType=c('ENTREZID', 'SYMBOL'), OrgDb = org.Hs.eg.db)
  enricher.all <- enricher(geneList$SYMBOL, TERM2GENE = tf2gene, pvalueCutoff = pvalueCutoff, minGSSize=NULL, maxGSSize=NULL)
  list(detail=enricher.all, tfs=as.data.frame(enricher.all)$ID)
}
enricher.all.bg1.001 <- colorEnricherAll(bg1,0.01)
enricher.all.bg2.001 <- colorEnricherAll(bg2,0.01)
enricher.all.bg3.001 <- colorEnricherAll(bg3,0.01)
enricher.all.bg0.001 <- colorEnricherAll(bg3,0.01)

enricher.all.bg1 <- colorEnricherAll(bg1)
enricher.all.bg2 <- colorEnricherAll(bg2)
enricher.all.bg3 <- colorEnricherAll(bg3)
enricher.all.bg0 <- colorEnricherAll(bg0)

### -------- export
as.data.frame(enricher.bg3.001$detail$brown)

total <- c()
colors <- unique(modules$color)
for(c in colors){
  print(c)
  result <- as.data.frame(enricher.bg3.001$detail[[c]])
  total <- c(total, nrow(result))
  write.csv(result, row.names = FALSE, file = paste0('./reports/tf_enricher/',c,'.csv'))
}
write.csv(as.data.frame(enricher.all.bg3.001$detail), row.names = FALSE, file = paste0('./reports/tf_enricher/all.csv'))


tf.n <- data.frame(module=colors, n=total)
modules.table <- read.delim('./reports/wgcna/geneModuleTable.csv', sep = ',')
write.csv(data.frame(modules.table, tf=tf.n[match(modules.table$X, tf.n$module),2]), file = './reports/tf_enricher/module_tf_num.csv', row.names = FALSE)


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
gseaer.all.bg0 <- colorGSEAAll(bg0)














