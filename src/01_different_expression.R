library(edgeR)
library(ProjectTemplate)

load.project()

data <- helper.getCount(TRUE)

counts <- data$data[,-1]
group <- data$group
genes <- bioMart[match(data$data$GeneID, bioMart$Gene.stable.ID),] %>%
        select(GeneID=Gene.stable.ID,
               GeneType=Gene.type,
               symbol=HGNC.symbol,
               GoTermName=GO.term.name)

y <- DGEList(counts = counts, group = data$group, genes = genes)
## filter data
keep <- rowSums(cpm(counts) > 0.5) >= ncol(counts)/2
y <- y[keep, keep.lib.sizes=FALSE]

system.time(y <- calcNormFactors(y))
system.time(y <- estimateCommonDisp(y))
system.time(y <- estimateTagwiseDisp(y))
system.time(et <- exactTest(y))

diff.all <- topTags(et, n=10000)$table
diff.all <- filter(diffGenes, (logFC > 0.7 | logFC < -0.7) & FDR < 0.01)
diff.lncRNA <- filter(diffGenes, GeneType %in% config$lncRNA)
diff.pcg <- filter(diffGenes, GeneType %in% config$PCGs)

write.csv(diff.all, file = './reports/diff.all.106.csv')
write.csv(diff.lncRNA, file = './reports/diff.lncRNA.106.csv')
write.csv(diff.pcg, file = './reports/diff.pcg.106.csv')
