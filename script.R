load('./cache/00_TF_targets.Rdata')
i <- 0
total <- nrow(gtrd.mid)
genes<-mapply(function(chrom, start, end){
    i <<- i + 1
    print(paste(i, total))
    tfbs <- biomart.map[[chrom]][up<=end&tss>=start,]
    if (nrow(tfbs)==0){
        gene <- NA
    }else{
        gene <- tfbs$gene[1]
    }
    gene
}, gtrd.mid$chrom, gtrd.mid$start, gtrd.mid$end)

save.image('./cache/00_TF_targets.2.Rdata')
