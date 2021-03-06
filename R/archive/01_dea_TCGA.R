library(edgeR)
library(stringr)
library(dplyr)
library(pcProfile)

source('./R/lib.R')

# source('./lib/globals.R')
# source('./lib/helpers.R')

getDE <- function(def, lfc=1, p.value=0.05) {
  fdr <- p.adjust(def$table$PValue, method = 'BH')
  is.de <- decideTestsDGE(def, lfc=lfc, p.value = p.value)
  de <- def$table[is.de!=0,]
  de.genes <- def$genes[match(rownames(de), rownames(def$genes)),]
  cbind(de.genes, de,FDR=fdr[is.de!=0])
}

data("prad.rna.count") # put prad.rna.count into the global environment
biomart <- distinct(pf.get.biomart(), ensembl_gene_id, .keep_all = TRUE)

# 
# data.sample <- read.delim('./data/data.sample.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)

counts <- prad.rna.count
group <- c(rep('T', 499), rep('N',52)) # 1:499 tumor, 500:551: normal, ?prad.ran.count fro detail.
genes <- biomart[match(rownames(prad.rna.count), biomart$ensembl_gene_id),] %>%
  dplyr::select(GeneID=ensembl_gene_id,
         GeneType=gene_biotype,
         symbol=hgnc_symbol)


y <- DGEList(counts = counts, genes = genes)
## filter data
keep <- rowSums(cpm(counts) > 0.5) >= ncol(counts)/10
y <- y[keep, keep.lib.sizes=FALSE]
## normalization
system.time(y <- calcNormFactors(y))
## ------------- quasi-likelihood -----------------------------
## estimate Disp
tissue <- factor(data.sample$sample.type)
design <- model.matrix(~tissue)
system.time(y <- estimateDisp(y, design, robust=TRUE)) # 483s
plotBCV(y)
# saveRDS(y, './cache/y.glm.rds')
## get GLM model
system.time(fit <- glmQLFit(y, design, robust=TRUE))
## test
system.time(qlf <- glmQLFTest(fit))
lfc <- log2(2)
p.value <- 0.05

is.de <- decideTests(qlf, lfc=lfc, p.value = p.value)
summary(is.de)
plotMD(qlf, status = is.de, col=c("red","blue"))

## the NA line is because not Every gene is RNA-seq is annotated in biomart
diff.qlf <- getDE(qlf, lfc=lfc, p.value = p.value) # <--------------- USE THIS 
diff.qlf.lncRNA <- filter(diff.qlf , GeneType %in% config$lncRNA)
diff.qlf.pcg <- filter(diff.qlf , GeneType %in% config$PCGs)
write.csv(diff.qlf, file = './data/diff.qlf.csv', row.names = FALSE)
write.csv(diff.qlf.lncRNA, file = './data/diff.qlf.lncRNA.csv', row.names = FALSE)
write.csv(diff.qlf.pcg, file = './data/diff.qlf.pcg.csv', row.names = FALSE)

## ------ compare classical and quasi-likelihood
diff.classical <- read.delim('./data/diff.classical.csv', sep = ',', stringsAsFactors = FALSE)
diff.classical.lncRNA <- read.delim('./data/diff.classical.lncRNA.csv', sep = ',', stringsAsFactors = FALSE)
diff.classical.pcg <- read.delim('./data/diff.classical.pcg.csv', sep = ',', stringsAsFactors = FALSE)

diff <- function(d1, d2) {
  print(dim(d1))
  print(dim(d2))
  print(length(intersect(d1$GeneID, d2$GeneID)))
}
diff(diff.qlf, diff.classical)
diff(diff.qlf.lncRNA, diff.classical.lncRNA)
diff(diff.qlf.pcg, diff.classical.pcg)
## ------ classical method
y <- DGEList(counts = counts, group = group, genes = genes)
## filter data
keep <- rowSums(cpm(counts) > 0.5) >= ncol(counts)/10
y <- y[keep, keep.lib.sizes=FALSE]

system.time(y <- calcNormFactors(y))
system.time(y <- estimateDisp(y))
plotBCV(y)

system.time(et <- exactTest(y))
plotMD(y)

diff.all <- topTags(et, n=10000)$table
diff.all <- filter(diffGenes, (logFC > 0.7 | logFC < -0.7) & FDR < 0.01)
diff.lncRNA <- filter(diffGenes, GeneType %in% config$lncRNA)
diff.pcg <- filter(diffGenes, GeneType %in% config$PCGs)

write.csv(diff.all, file = './reports/diff.all.106.csv')
write.csv(diff.lncRNA, file = './reports/diff.lncRNA.106.csv')
write.csv(diff.pcg, file = './reports/diff.pcg.106.csv')
