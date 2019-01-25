source('./R/lib.R')
library(data.table)
library(stringr)

tf2gene <- function(pairs) {
    lapply(split(pairs,as.vector(pairs$tf)), function(x){unique(as.vector(x$gene))})
}

### Fimo TSS
read.fimo.460.tss <- function() {
    ret <- data.table()
    for(i in 1:4) {
        fimo <- fread(paste0('./data/fimo.460.tss/fimo.460.tss.1000.', i, '.tsv'), check.names = TRUE)
        sequence_name <- str_split(fimo$sequence_name, '\\|', simplify = TRUE)
        gene <- sequence_name[,1]
        transcript <- sequence_name[,2]
        chromosome <- sequence_name[,3]
        fimo[,c('gene','transcript','chromosome'):=list(gene, transcript, chromosome)]
        fimo_uniq <- unique(fimo, by=c('motif_id', 'start', 'strand', 'chromosome'))
        ret <- rbind(ret, fimo_uniq)
    }
    ret
}
fimo <- read.fimo.460.tss()
fimo.tss.set <- fimo[p.value<1e-6, .(.N), by=.(motif_alt_id, gene)][,.(tf=motif_alt_id, gene)]
fimo.tss.set[,.(.N),by=(tf)]
saveRDS(as.data.frame(fimo.tss.set), './cache/fimo.tss.set.rds')
tf2gene.fimo.tss <- tf2gene(fimo.tss.set)

### Fimo GSS
get.tf.2.PCG.from.fimo <- function() {
  env = globalenv()
  key = '.tf.2.PCG'
  if (exists(key, envir = env)) {
    get(key, envir = env)
  } else {
    fimo <- read.delim('./data/enricher/fimo.460.2230.tsv', stringsAsFactors = F)
    fimo <- filter(fimo, motif_alt_id!='')
    tf.2.PCG <- lapply(split(fimo,as.vector(fimo$motif_alt_id)), function(x){unique(as.vector(x$sequence_name))})
    assign(key, tf.2.PCG, envir = env)
    tf.2.PCG
  }
}
tf2gene.fimo.gss <- get.tf.2.PCG.from.fimo()
### TRUST
trrust<-fread('./data/trrust_rawdata.human.tsv')
trrust.set <- trrust[,.(tf=V1,gene=pf.symbol2emsembl(V2))]
trrust.set[,.(.N),by=(tf)]
saveRDS(as.data.frame(trrust.set), './cache/trrust.set.rds')
tf2gene.trrust <- tf2gene(trrust.set)


### Gtrd
gtrd <- fread('./data/Homo_sapiens_meta_clusters.interval', check.names = TRUE)
gtrd.mid <- gtrd[,.(chrom=str_sub(X.CHROM, 4), start=START, end=END, len=END-START, tf=tfTitle)]

biomart.set <- as.data.table(pf.get.biomart())
biomart.set <- biomart.set[chromosome_name%in%gtrd.mid$chrom,][,.(gene=ensembl_gene_id, transcript=ensembl_transcript_id, tss=transcription_start_site, up=transcription_start_site-1000,chrom=chromosome_name)]
biomart.map<-split(biomart.set, by = 'chrom')

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
gtid.mid <- readRDS('./cache/gtrd.mid.rds')
gtrd.set<-unique(gtrd.mid[!is.na(gene), .(tf, gene)])
gtrd.set[,.(.N),by=(tf)]
saveRDS(as.data.frame(gtrd.set), './cache/gtrd.set.rds')

tf2gene.gtrd <- tf2gene(gtrd.set)
##################
par(mfrow=c(2,2))
plot(sort(sapply(tf2gene.fimo.tss, length)), ylab = 'pcg number', xlab = 'fimo1')
plot(sort(sapply(tf2gene.fimo.gss, length)), ylab = 'pcg number', xlab = 'fimo0')
plot(sort(sapply(tf2gene.trrust, length)), ylab = 'pcg number', xlab = 'trrust')
plot(sort(sapply(tf2gene.gtrd, length)), ylab = 'pcg number', xlab = 'gtrd')
