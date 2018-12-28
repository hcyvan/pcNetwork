library(clusterProfiler)
library(org.Hs.eg.db)

source('./lib/globals.R')
source('./lib/helpers.R')
load_all('./package/x2y/')

cor.pairs.info <- readRDS('./cache/cor.pairs.info.rds')
biomart <- helper.get.biomart()
cor.lnc2all <- filter(cor.pairs.info, type1%in%config$lncRNA, type2%in%config$PCGs,FDR<0.05,abs(r)>=0.3)



load('./cache/lncTP.0.x.rda')
candidate <- helper.get.candidate()

cor.pairs.info <- readRDS('./cache/cor.pairs.info.rds')
biomart <- helper.get.biomart()
cor.lnc2all <- filter(cor.pairs.info,
                      type1%in%config$lncRNA, type2%in%config$PCGs,FDR<0.05,abs(r)>=0.3)


surv <- readRDS(file='./cache/candidate.surv.rds')
lnc.tfpcg <- getYZByX(lncTP.0.3)

options(digits = 3)
for (ca in candidate$name) {
    fd.km.lnc <- ifelse(ca%in%surv$fd.km$lncRNA,1,0)
    fd.km.pcg <- length(intersect(lnc.tfpcg[[ca]][[1]], surv$fd.km$pcg))/length(lnc.tfpcg[[ca]][[1]])
    o.km.lnc <- ifelse(ca%in%surv$o.km$lncRNA,1,0)
    o.km.pcg <- length(intersect(lnc.tfpcg[[ca]][[1]], surv$o.km$pcg))/length(lnc.tfpcg[[ca]][[1]])
    fd.cox.lnc <- ifelse(ca%in%surv$fd.cox$lncRNA,1,0)
    fd.cox.pcg <- length(intersect(lnc.tfpcg[[ca]][[1]], surv$fd.cox$pcg))/length(lnc.tfpcg[[ca]][[1]])
    o.cox.lnc <- ifelse(ca%in%surv$o.cox$lncRNA,1,0)
    o.cox.pcg <- length(intersect(lnc.tfpcg[[ca]][[1]], surv$o.cox$pcg))/length(lnc.tfpcg[[ca]][[1]])
    print(ca)
    print(c(fd.km.lnc, fd.km.pcg, o.km.lnc, o.km.pcg))
    print(c(fd.cox.lnc, fd.cox.pcg, o.cox.lnc, o.cox.pcg))
}

str(lncTP.0.3@raw)


lncTP.filterd<-filterByX(lncTP.0.3, candidate$GeneID)
lncTP.0.3
lncTP.filterd


lncTP.filterd@nodes$x
lncTP.filterd@nodes$y


enrich.all <- function(lnc) {
    set <- as.vector(filter(cor.lnc2all, v1==lnc)$v2)
    gene <- bitr(set, fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)


    gene.sorted <- filter(cor.lnc2all,v1==lnc)%>%select(v2,r)%>%arrange(desc(r))
    gene.names <-  bitr(gene.sorted$v2, fromType = "ENSEMBL",
                        toType = c("ENTREZID", "SYMBOL"),
                        OrgDb = org.Hs.eg.db)
    gene.names <- distinct(gene.names, ENSEMBL, .keep_all = TRUE)
    gene.joined <- inner_join(gene.sorted, gene.names, by=c('v2'='ENSEMBL'))
    geneList <- gene.joined$r
    names(geneList) <- gene.joined$ENTREZID

    ego.MF <- enrichGO(gene=gene$ENTREZID, OrgDb=org.Hs.eg.db, ont='MF')
    ego.BP <- enrichGO(gene=gene$ENTREZID, OrgDb=org.Hs.eg.db, ont='BP')
    ego.CC <- enrichGO(gene=gene$ENTREZID, OrgDb=org.Hs.eg.db, ont='CC')
    ekegg.path <- enrichKEGG(gene=gene$ENTREZID)

    ggo.MF <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont='MF')
    ggo.BP <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont='BP')
    ggo.CC <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont='CC')
    gkegg.path <- gseKEGG(geneList = geneList)

    list(
        ego.MF=ego.MF,
        ego.BP=ego.BP,
        ego.CC=ego.CC,
        ekegg.path=ekegg.path,
        ggo.MF=ggo.MF,
        ggo.BP=ggo.BP,
        ggo.CC=ggo.CC,
        gkegg.path=gkegg.path
    )
}

all.enrich.result <- list()
for (lnc in candidate$name) {
    all.enrich.result[[lnc]] <- enrich.all(lnc)
}

saveRDS(all.enrich.result, file = './cache/all.enrich.result.rds')


candidate
lnc <- 'ENSG00000196366'
lnctp <- filterByX(lncTP.0.3, lnc)
#set <- lnctp@nodes$z
set <- as.vector(filter(cor.lnc2all, v1==lnc)$v2)
########################################################### GO
gene <- bitr(set, fromType = "ENSEMBL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
ego.MF <- enrichGO(gene=gene$ENTREZID, OrgDb=org.Hs.eg.db, ont = 'MF')
ego.BP <- enrichGO(gene=gene$ENTREZID, OrgDb=org.Hs.eg.db, ont = 'BP')
ego.CC <- enrichGO(gene=gene$ENTREZID, OrgDb=org.Hs.eg.db, ont = 'CC')
ekegg.path <- enrichKEGG(gene=gene$ENTREZID)
########################################################### GSEA
gene.sorted <- filter(cor.lnc2all,v1==lnc)%>%select(v2,r)%>%arrange(desc(r))
gene.names <-  bitr(gene.sorted$v2, fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db)
gene.names <- distinct(gene.names, ENSEMBL, .keep_all = TRUE)
gene.joined <- inner_join(gene.sorted, gene.names, by=c('v2'='ENSEMBL'))
geneList <- gene.joined$r
names(geneList) <- gene.joined$ENTREZID
ggo <- gseGO(geneList = geneList, OrgDb=org.Hs.eg.db)
gkegg <- gseKEGG(geneList = geneList)

head(ggo)


              
bitr(candidate$name, OrgDb = org.Hs.eg.db, fromType = "ENSEMBL", toType=c("ENTREZID"))

###########################
 library(org.Hs.eg.db)
data(geneList, package = 'DOSE')
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
head(gene.df)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
?dropGO
?gofilter
getGOLevel


library(KEGGREST)
listDatabases()
org <- keggList("organism")
org[,1:2]
head(org)

queryables <- c(listDatabases(), org[,1], org[,2])
keggList('hsa')
hsa <- keggList('hsa')


hsa[1]
query <- keggGet(c("hsa:10458", "ece:Z5100"))
str(query[[1]])


?enrichKEGG
