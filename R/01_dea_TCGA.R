library(edgeR)
library(stringr)
library(dplyr)
library(pcProfile)

source('./R/lib.R')

data("prad.rna.count") # put prad.rna.count into the global environment
biomart <- distinct(pf.get.biomart(), ensembl_gene_id, .keep_all = TRUE)

counts <- prad.rna.count[,c(500:551,1:499)]# 1:499 tumor, 500:551: normal, ?prad.ran.count fro detail.
group <- c(rep('N',52),rep('T', 499))
genes <- biomart[match(rownames(prad.rna.count), biomart$ensembl_gene_id),] %>%
  dplyr::select(GeneID=ensembl_gene_id,
         GeneType=gene_biotype,
         symbol=hgnc_symbol)


y <- DGEList(counts = counts, genes = genes)
## filter data
keep <- rowSums(cpm(counts) > 0.5) >= ncol(counts)/10
y <- y[keep, keep.lib.sizes=FALSE]
## normalization: TMM
system.time(y <- calcNormFactors(y))
## ------------- quasi-likelihood -----------------------
design <- model.matrix(~factor(group))
## estimate Disp
system.time(y <- estimateDisp(y, design, robust=TRUE))
## get GLM model
system.time(fit <- glmQLFit(y, design, robust=TRUE))
## test
system.time(qlf <- glmQLFTest(fit, coef = 2))

qlf.top <- topTags(qlf, 10000)$table
diff <- filter(qlf.top,abs(logFC)>1 & FDR < 0.05)

saveRDS(diff, './support/diff.T_N.qlf.1.005.3069.rds')
a <- pf.get.diff('pcg')
