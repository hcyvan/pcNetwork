source('./lib/helpers.R')

library(WGCNA)

diff.gene <- pf.get.diff()
zfpkm <- pf.filter.zfpkm(diff.gene$GeneID)
biomart <- pf.filter.anno(diff.gene$GeneID)

# ----------------------------- WGCNA -------------------------------------
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

datExpr <- t(as.matrix(zfpkm))
# ----------------------------
corType <- "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers <- ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
# networkType <- 'unsigned'
networkType <- 'signed hybrid'

resultPath <- './reports/wgcna/'
TOMfile <- './cache/diff.qlf.2877.tom'

## -------------------------- Estimate Soft Power -------------------------
## Use this value: stf$powerEstimate
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr,
                        networkType = networkType,
                        corFnc = corFnc,
                        powerVector = powers,
                        verbose = 5)
# tiff(paste0(resultPath, 'sft.tiff'), width = 862, height = 571)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topoloy Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.9, col='red')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()
## ------------------------- Find Modules ----------------------------
system.time(
  net <- blockwiseModules(datExpr,
                          power = sft$powerEstimate,
                          corType = corType,
                          maxPOutliers = maxPOutliers,
                          networkType = networkType,
                          minModuleSize = 30,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = TOMfile,
                          verbose = 3)
)

load(net$TOMFiles)
load('./cache/diff.qlf.2877.tom-block.1.RData')
tom <- as.matrix(TOM)
saveRDS(tom, file = './cache/tom.rds')
rm(TOM)

## -------------------------------
sizeGrWindow(2, 9)
colors <- labels2colors(net$colors)
# tiff(paste0(resultPath, 'clusterDendrogram.tiff'), width = 862, height = 571)
plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]], 'Moudle colors',
                    dendroLabels = FALSE,
                    hang= 0.3,
                    addGuide = TRUE,
                    guideHang = 0.05)
# dev.off()

## ----    Module Eigengene relationship / Eigengene Networks ---------
MEs <- moduleEigengenes(datExpr, colors, softPower = sft$powerEstimate)$eigengenes
MEs <- orderMEs(MEs)
tiff(paste0(resultPath, 'ModuleRelationship.tiff'), width = 862, height = 571)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,16,2,16),
                      marHeatmap = c(3,16,2,16), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
## ---- gene VS. module -------
MGs <- list()
MGs.col <- data.frame()
for(c in unique(colors)){
  print(c)
  datExpr_c = datExpr[,colors==c]
  cor_c <- cor(datExpr_c, MEs[[paste0('ME',c)]])
  colnames(cor_c) <- 'cor'
  MGs[[c]] <- cor_c
  MGs.col <- rbind(MGs.col, data.frame(mg.cor=cor_c[,1]))
}
## --------------------------- Export result -------------------------
summaryTable <- function(data) {
  data <- as.data.frame.matrix(data)
  lncRNA <- rowSums(data[intersect(colnames(data), config$lncRNA)])
  PCG <- rowSums(data[intersect(colnames(data), config$PCGs)])
  data.frame(lncRNA=lncRNA, PCG=PCG)
}

gene.module <- data.frame(diff.gene, colors=colors)
gene.module <- data.frame(gene.module, mg.cor=MGs.col$mg.cor[match(gene.module$GeneID, rownames(MGs.col))])
gene.module.summary <- summaryTable(table(gene.module[c('colors', 'GeneType')]))
write.csv(gene.module, file = paste0(resultPath, 'geneModule.csv'))
write.csv(gene.module.summary , file=paste0(resultPath, 'geneModuleTable.csv'))
saveRDS(data.frame(id=gene.module$GeneID, type=gene.module$GeneType,colors=gene.module$colors, mg.cor=gene.module$mg.cor),
          file ='./cache/wgcna.colors.rds')



###################################################

tom.2 = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)
tom.2 = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate, networkType = 'signed')

adjacency.2 = adjacency(datExpr, power = sft$powerEstimate,type='signed')
adjacency = adjacency(datExpr, power = sft$powerEstimate)
tom.3 <- TOMsimilarity(adjacency.2)
tom.4 <- TOMsimilarity(adjacency, TOMType = 'signed')












