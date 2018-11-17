library(WGCNA)

source('./lib/globals.R')
source('./lib/helpers.R')

diff.gene <- helper.get.lncRNA.PCG()
data.fpkm <- helper.get.fpkm.count()
biomart <- helper.get.biomart()

genes.fpkm <- data.fpkm[match(diff.gene$GeneID,rownames(data.fpkm)),]
datExpr <- t(as.matrix(genes.fpkm))
# ----------------------------- WGCNA -------------------------------------
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ----------------------------
corType <- "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers <- ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
networkType <- 'signed'

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
tiff(paste0(resultPath, 'sft.tiff'), width = 862, height = 571)
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
dev.off()
## ------------------------- Find Modules ----------------------------
system.time(
  net <- blockwiseModules(datExpr,
                          power = sft$powerEstimate,
                          corType = corType,
                          maxPOutliers = maxPOutliers,
                          networkType = networkType,
                          TOMType = 'signed', 
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
rm(TOM)

## -------------------------------
sizeGrWindow(2, 9)
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
tiff(paste0(resultPath, 'clusterDendrogram.tiff'), width = 862, height = 571)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 'Moudle colors',
                    dendroLabels = FALSE,
                    hang= 0.3,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

## ----    Module Eigengene relationship / Eigengene Networks ---------
MEs <- moduleEigengenes(datExpr,
                        moduleColors,
                        softPower = sft$powerEstimate)$eigengenes
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
for(c in unique(moduleColors)){
  print(c)
  datExpr_c = datExpr[,moduleColors==c]
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

gene.module <- data.frame(diff.gene, moduleColor=moduleColors)
gene.module <- data.frame(gene.module, mg.cor=MGs.col$mg.cor[match(gene.module$GeneID, rownames(MGs.col))])
gene.module.summary <- summaryTable(table(gene.module[c('moduleColor', 'GeneType')]))
write.csv(gene.module, file = paste0(resultPath, 'geneModule.csv'))
write.csv(gene.module.summary , file=paste0(resultPath, 'geneModuleTable.csv'))
write.csv(data.frame(id=gene.module$GeneID, type=gene.module$GeneType,color=gene.module$moduleColor, mg.cor=gene.module$mg.cor),
          row.names = FALSE,
          file ='./data/diff.qlf.2877.wgcna.color.csv')


################################ step by step #####################################
softPower <- sft$powerEstimate
system.time(adjacency <- adjacency(datExpr, power = softPower))



