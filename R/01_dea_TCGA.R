library(edgeR)
library(stringr)
library(dplyr)
library(pcProfile)
library(ggplot2)
library(gplots)
source('./R/lib.R')

biomart <- distinct(pf.get.biomart(), ensembl_gene_id, .keep_all = TRUE)
counts <- pf.get.count()[,c(500:551,1:499)]# 1:499 tumor, 500:551: normal
group <- c(rep('N',52),rep('T', 499))
genes <- biomart[match(rownames(counts), biomart$ensembl_gene_id),] %>%
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
diff <- filter(qlf.top, !is.na(GeneID), abs(logFC)>1 & FDR < 0.05)

saveRDS(diff, './support/diff.T_N.qlf.1.005.3069.rds')
# write.csv(diff, './data/diff.T_N.qlf.1.005.3069.csv', col.names = FALSE)

#############################################################################################
##################################### Analysis #############################################
###########################################################################################
# diff gene statistics
diff.pcg <- filter(diff, GeneType%in%pv.pcg)
diff.lncRNA <- filter(diff, GeneType%in%pv.lncRNA)
diff.other <- filter(diff, !(GeneType%in%c(pv.lncRNA,pv.pcg)))

dim(filter(diff, logFC>0))
dim(filter(diff, logFC<0))

n.all <- nrow(diff)
n.pcg <- nrow(diff.pcg)
n.lncRNA <- nrow(diff.lncRNA)
n.other <- nrow(diff.other)
n.all
n.pcg
n.lncRNA
n.other
# volcano
dd <- data.frame(qlf$table, FDR=p.adjust(qlf$table$PValue,method = 'BH'),symbol=qlf$genes$symbol)
png('./reports/thesis/tcga_diff_gene_valcano.png',width=600,height = 400)
ggplot(data=dd) + geom_point(aes(x=logFC,y=-log(FDR),color=logCPM)) +
  scale_colour_gradientn(colours=c("#FF0000" ,"#FFFF00" )) +
  geom_hline(aes(yintercept=-log(0.05)),linetype="dashed") +
  geom_vline(aes(xintercept=1),linetype="dashed") +
  geom_vline(aes(xintercept=-1),linetype="dashed") +
  annotate('text',x=-8,y=100,label='Downregulate\nlogFC < -1, FDR < 0.05',size=4)+
  annotate('text',x=5,y=200,label='Upregulate\nlogFC > 1, FDR < 0.05',size=4)+
  annotate("rect", xmin=1, xmax=Inf, ymin=-log(0.05), ymax=Inf,alpha=.1,fill="red") +
  annotate("rect", xmin=-Inf, xmax=-1, ymin=-log(0.05), ymax=Inf,alpha=.1,fill="green") +
  ggtitle('Tumor VS Normal')
dev.off()

#Diff Gene Heatmap
normal <- samples[500:551]
tumor <- samples[1:499]
tumor.2 <- c()
for (i in tumor) {
  for(j in normal) {
    if (str_sub(i,1,12)==str_sub(j,1,12)) {
      tumor.2 <- c(tumor.2, i)
      break
    }
  }
}
diff <- arrange(diff, desc(abs(logFC)))
diff.fpkm <- pf.filter.fpkm(diff$GeneID, rm.na = TRUE)
diff.fpkm.filter = diff.fpkm[apply(diff.fpkm==0,1,sum)==0,]
png(filename=paste0('./reports/thesis/tcga_diff_gene_heatmap.png'),width=1200,height=600)
m <- log2(as.matrix(diff.fpkm.filter)[,c(tumor.2, normal)])
a<-heatmap.2(m,
             key.title = NA,
             key.xlab = 'Z-Score of log2(FPKM)',
             key.ylab = NA,
             colsep=54,
             sepwidth = c(0.5,0.5),
            trace = 'none',
            scale = 'row',
            Colv=FALSE,
            lhei = c(1,4),
            col=greenred)
dev.off()
