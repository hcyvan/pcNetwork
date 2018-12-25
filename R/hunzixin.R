library(WGCNA)

enableWGCNAThreads()


load('~/文档/huzixin/expr_for_WGCNA.rda')

datExpr <- expr_all

corType <- "pearson"
# corType <- "bicor"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers <- ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
networkType <- 'unsigned'

resultPath <- './reports/huzixin/'
TOMfile <- file.path(resultPath,'tom')

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
## --------------------------------------------------------------------------
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
#load(net$TOMFiles)
load("./reports/huzixin//tom-block.1.RData") # TOM
## -------------------------------------------------------------------------
sizeGrWindow(2, 9)
colors <- labels2colors(net$colors)
tiff(paste0(resultPath, 'clusterDendrogram.tiff'), width = 862, height = 571)
plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]], 'Moudle colors',
                    dendroLabels = FALSE,
                    hang= 0.3,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()
