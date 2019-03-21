source('./R/lib.R')
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)


diff <- pf.get.diff()
diff.pcg <- pf.get.diff('pcg')

pcg = bitr(diff.pcg$GeneID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego.CC <- enrichGO(gene          = pcg$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego.MF <- enrichGO(gene          = pcg$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego.BP <- enrichGO(gene          = pcg$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ek.pathway <- enrichKEGG(gene         = pcg$ENTREZID,
                         organism     = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

# saveRDS(ego.MF, file = './data/ego.MF.rds')
# saveRDS(ego.BP, file = './data/ego.BP.rds')
# saveRDS(ego.CC, file = './data/ego.CC.rds')
# saveRDS(ek.pathway, file = './data/ek.pathway.rds')

# ego.MF <- readRDS('./data/ego.MF.rds')
# ego.BP <- readRDS(file = './data/ego.BP.rds')
# ego.CC <- readRDS(file = './data/ego.CC.rds')
# ek.pathway <- readRDS(file = './data/ek.pathway.rds')

dim(head(ego.MF,n=100)[,-8])
dim(head(ego.BP,n=1000)[,-8])
dim(head(ego.CC,n=100)[,-8])
dim(head(ek.pathway,n=100)[,-8])
#######################################################
################ simplify
# https://github.com/GuangchuangYu/clusterProfiler/issues/28
######################################################
# ego.MF.sim<-simplify(ego.MF)
# ego.BP.sim<-simplify(ego.BP)
# ego.CC.sim<-simplify(ego.CC)
# 
# saveRDS(ego.MF.sim, file = './data/ego.MF.sim.rds')
# saveRDS(ego.BP.sim, file = './data/ego.BP.sim.rds')
# saveRDS(ego.CC.sim, file = './data/ego.CC.sim.rds')

ego.MF.sim <- readRDS(file = './data/ego.MF.sim.rds')
ego.BP.sim <- readRDS(file = './data/ego.BP.sim.rds')
ego.CC.sim<-readRDS(file = './data/ego.CC.sim.rds')

dim(head(ego.MF.sim,n=100)[,-8])
dim(head(ego.BP.sim,n=1000)[,-8])
dim(head(ego.CC.sim,n=100)[,-8])
dim(head(ek.pathway,n=100)[,-8])
##################################### 
############### plot
p1<-dotplot(ego.MF.sim)
p2<-dotplot(ego.BP.sim)
p3<-dotplot(ego.CC.sim)
p4<-dotplot(ek.pathway)
win.metafile(filename="./reports/thesis/diff.enrich.emf",width=9,height=14)
plot_grid(p1, p2, p3, p4, labels = c('A','B','C','D'), cols=1,label_size = 20)
dev.off()
# 
# plotGOgraph(ego.MF.sim,firstSigNodes=100,sigForAll=FALSE)
# plotGOgraph(ego.BP.sim,firstSigNodes=100,sigForAll=FALSE)
# plotGOgraph(ego.CC.sim,firstSigNodes=100,sigForAll=FALSE)



pcg = bitr(pf.symbol2emsembl(ac$ID), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ek.pathway <- enrichKEGG(gene         = pcg$ENTREZID,
                         organism     = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)

######################################################
########################### AR time course
#####################################################
ar.diff <- pf.get.diff.ar()
################################ early
early = bitr(pf.symbol2emsembl(ar.diff$early),
           fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
late = bitr(pf.symbol2emsembl(ar.diff$late),
             fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego.CC <- enrichGO(gene          = early$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
ego.MF <- enrichGO(gene          = early$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego.BP <- enrichGO(gene          = early$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ek.pathway <- enrichKEGG(gene         = early$ENTREZID,
                         organism     = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
ego.MF.sim<-simplify(ego.MF)
ego.BP.sim<-simplify(ego.BP)
ego.CC.sim<-simplify(ego.CC)
saveRDS(ego.MF.sim, file = './data/ar.early.ego.MF.sim.rds')
saveRDS(ego.BP.sim, file = './data/ar.early.ego.BP.sim.rds')
saveRDS(ego.CC.sim, file = './data/ar.early.ego.CC.sim.rds')


##########################################
diff.ar<-pf.get.diff.ar()
early<-diff.ar$early
late<-diff.ar$late
inter<-diff.ar$inter
early.only<-setdiff(early,inter)
late.only<-setdiff(late,inter)
early = bitr(pf.symbol2emsembl(early.only),
             fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
late = bitr(pf.symbol2emsembl(late.only),
            fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gcSample<-list(Early=early$ENTREZID,Late=late$ENTREZID)
ck.pathway <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG")
ck.CC <- compareCluster(geneCluster = gcSample, fun = "enrichGO",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
ck.MF <- compareCluster(geneCluster = gcSample, fun = "enrichGO",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        minGSSize = 0,
                        readable      = TRUE)
ck.BP <- compareCluster(geneCluster = gcSample, fun = "enrichGO",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)
############### plot
p1<-dotplot(ck.BP)
p2<-dotplot(ck.CC)
p3<-dotplot(ck.pathway)
win.metafile(filename="./reports/thesis/ar.early.late.enrich.emf",width=10,height=9)
plot_grid(p1, p2, p3, labels = c('A','B','C'), cols=1,label_size = 20)
dev.off()
