source('./lib/helpers.R')
source('./R/lib.R')
devtools::load_all('./package/x2y/')

############################################### GET x.2.y
# plot.lncRNA.2.PCG()
# plot.tf.2.PCG()
################################################# tf-pcg
# get.tf.2.PCG.from.enricher <- function() {
#   tf.enricher <- read.delim('./data/enricher.all.bg0.new.csv', sep = ',', stringsAsFactors = F)
#   tf.enricher.list <- setNames(split(tf.enricher, seq(nrow(tf.enricher))), tf.enricher$detail.ID)
#
#   load('./cache/biomart.symbol.biotype.rda')
#
#   lapply(tf.enricher.list, function(tf){
#     symbol <- str_split(tf$detail.geneID, '/')[[1]]
#     filter(biomart.symbol.biotype, hgnc_symbol%in%symbol)$ensembl_gene_id
#   })
# }
#############################################################
pcg <- pf.get.diff('pcg')
lncRNA <- pf.get.diff('lncRNA')

fimo.gss.pairs <- readRDS('./cache/fimo.gss.set.rds')
fimo.tss.pairs <- readRDS('./cache/fimo.tss.set.rds')
trrust.pairs <- readRDS('./cache/trrust.set.rds')
gtrd.pairs <- readRDS('./cache/gtrd.set.rds')
gtrd.pc.pairs <- readRDS('./cache/gtrd.set.pc.rds')

cor.pairs.format <- function(s=0.3) {
    biomart <- pf.get.biomart()
    cor.pairs <- readRDS('./cache/cor.pearson.pairs.rds')
    t1 <- biomart[match(cor.pairs$v1, biomart$ensembl_gene_id), 'gene_biotype']
    t2 <- biomart[match(cor.pairs$v2, biomart$ensembl_gene_id), 'gene_biotype']
    data.frame(cor.pairs,t1=t1,t2=t2) %>%
        filter((t1%in%pv.lncRNA & t2%in%pv.pcg), FDR<0.05, abs(r) >= s) %>%
        select(lncRNA=v1, gene=v2)
}
cp.3 <- cor.pairs.format(s=0.3)
cp.4 <- cor.pairs.format(s=0.4)
cp.5 <- cor.pairs.format(s=0.5)
cp.8 <- cor.pairs.format(s=0.8)

fimo.gss.3 <- dEnricher(cp.3, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.3.xy2z.rds')
fimo.tss.3 <- dEnricher(cp.3, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.3.xy2z.rds')
trrust.3 <- dEnricher(cp.3, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.3.xy2z.rds')
gtrd.3 <- dEnricher(cp.3, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.3.xy2z.rds')
gtrd.pc.3 <- dEnricher(cp.3, gtrd.pc.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.pc.cp.3.xy2z.rds')

fimo.gss.4 <- dEnricher(cp.4, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.4.xy2z.rds')
fimo.tss.4 <- dEnricher(cp.4, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.4.xy2z.rds')
trrust.4 <- dEnricher(cp.4, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.4.xy2z.rds')
gtrd.4 <- dEnricher(cp.4, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.4.xy2z.rds')
gtrd.pc.4 <- dEnricher(cp.4, gtrd.pc.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.pc.cp.4.xy2z.rds')

fimo.gss.5 <- dEnricher(cp.5, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.5.xy2z.rds')
fimo.tss.5 <- dEnricher(cp.5, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.5.xy2z.rds')
trrust.5 <- dEnricher(cp.5, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.5.xy2z.rds')
gtrd.5 <- dEnricher(cp.5, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.5.xy2z.rds')
gtrd.pc.5 <- dEnricher(cp.5, gtrd.pc.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.pc.cp.5.xy2z.rds')

fimo.gss.8 <- dEnricher(cp.8, fimo.gss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.gss.cp.8.xy2z.rds')
fimo.tss.8 <- dEnricher(cp.8, fimo.tss.pairs, pcg$GeneID, rds='./cache/xy2z/fimo.tss.cp.8.xy2z.rds')
trrust.8 <- dEnricher(cp.8, trrust.pairs, pcg$GeneID, rds='./cache/xy2z/trrust.cp.8.xy2z.rds')
gtrd.8 <- dEnricher(cp.8, gtrd.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.cp.8.xy2z.rds')
gtrd.pc.8 <- dEnricher(cp.8, gtrd.pc.pairs, pcg$GeneID, rds='./cache/xy2z/gtrd.pc.cp.8.xy2z.rds')

#filter(lncTP.gtrd@detail, b=='MYC', a%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
