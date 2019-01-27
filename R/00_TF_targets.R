source('./R/lib.R')
library(data.table)
library(stringr)

tf2gene <- function(pairs) {
    lapply(split(pairs,as.vector(pairs$tf)), function(x){unique(as.vector(x$gene))})
}

### Fimo TSS
fimo.tss.1 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.1.tsv'), check.names = TRUE)
fimo.tss.2 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.2.tsv'), check.names = TRUE)
fimo.tss.3 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.3.tsv'), check.names = TRUE)
fimo.tss.4 <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.4.tsv'), check.names = TRUE)
fimo.tss <- rbind(fimo.tss.1, fimo.tss.2, fimo.tss.3, fimo.tss.4)
name <- str_split(fimo.tss$sequence_name, '\\|', simplify = TRUE)
fimo.tss[,c('gene','transcript','chromosome'):=list(name[,1], name[2], name[,3])]
fimo.tss.set <- unique(fimo.tss[p.value<=1e-5,.(tf=motif_alt_id, gene=gene)])
saveRDS(as.data.frame(fimo.tss.set), './cache/fimo.tss.set.rds')
### Fimo GSS
fimo.gss <- fread('./data/fimo.460.gss.tsv', stringsAsFactors = F, check.names = TRUE)
fimo.gss.set <- unique(fimo.gss[p.value<1e-5, .(tf=motif_alt_id, gene=sequence_name)])
saveRDS(as.data.frame(fimo.tss.set), './cache/fimo.gss.set.rds')
### TRUST
trrust<-fread('./data/trrust_rawdata.human.tsv')
trrust.set <- unique(trrust[,.(tf=V1,gene=pf.symbol2emsembl(V2))])
saveRDS(as.data.frame(trrust.set), './cache/trrust.set.rds')
### Gtrd
#gtrd <- fread('./data/Homo_sapiens_meta_clusters.interval', check.names = TRUE)
#gtrd.mid <- gtrd[,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]
#biomart.set <- as.data.table(pf.get.biomart())
#biomart.set <- biomart.set[chromosome_name%in%gtrd.mid$chrom,][,.(gene=ensembl_gene_id, transcript=ensembl_transcript_id, tss=transcription_start_site, up=transcription_start_site-1000,chrom=chromosome_name)]
#biomart.map<-split(biomart.set, by = 'chrom')
#i <- 0
#total <- nrow(gtrd.mid)
#genes<-mapply(function(chrom, start, end){
#    i <<- i + 1
#    print(paste(i, total))
#    tfbs <- biomart.map[[chrom]][up<=end&tss>=start,]
#    if (nrow(tfbs)==0){
#        gene <- NA
#    }else{
#        gene <- tfbs$gene[1]
#    }
#    gene
#}, gtrd.mid$chrom, gtrd.mid$start, gtrd.mid$end)
#gtrd.mid[,("gene"):=genes]
#saveRDS(gtrd.mid, './cache/gtrd.mid.rds')
gtrd.mid <- readRDS('./cache/gtrd.mid.rds')
gtrd.set<-unique(gtrd.mid[!is.na(gene), .(tf, gene)])
saveRDS(as.data.frame(gtrd.set), './cache/gtrd.set.rds')
gtrd.set <- as.data.table(readRDS('./cache/gtrd.set.rds'))
##################
library(ggplot2)
library(cowplot)

fimo.tss.hist <- ggplot(fimo.tss.set[,.(.N), by=(tf)], aes(x=N))+geom_histogram(binwidth = 100)
fimo.gss.hist <- ggplot(fimo.gss.set[,.(.N), by=(tf)], aes(x=N))+geom_histogram(binwidth = 100)
trrust.hist <- ggplot(trrust.set[,.(.N),by=(tf)], aes(x=N))+geom_histogram(binwidth = 1)
gtrd.hist <- ggplot(gtrd.set[,.(.N),by=(tf)], aes(x=N))+geom_histogram(binwidth = 100)
plot_grid(fimo.tss.hist, fimo.gss.hist, trrust.hist, gtrd.hist, labels = c(1,2,3,4))
