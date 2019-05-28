source('R/lib.R')
############ change hg19 to GRCh38
# use this webtool to change hg19 to GRCh38 <http://uswest.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core>
ARA_lincRNA.38<-read.csv('./data/lncRNA2TF/ARA_lincRNAs.hg38.bed',sep = '\t',header = FALSE)
colnames(ARA_lincRNA.38)<-c('chrom','chromStart','chromEnd','name','score','strand')

############################
biomart<-pf.get.biomart()
biomart.split<-split(biomart,as.vector(biomart$chromosome_name))
ARA_lincRNA.38.split<-split(ARA_lincRNA.38,as.vector(ARA_lincRNA.38$chrom))
names(ARA_lincRNA.38.split)
ARAlincRNA38En<-do.call(rbind,lapply(names(ARA_lincRNA.38.split), function(x){
  print(x)
  bs<-biomart.split[[x]]
  aras<-ARA_lincRNA.38.split[[x]]
  do.call(rbind,lapply(split(aras,seq(nrow(aras))), function(y){
    m<-lapply(split(bs,seq(nrow(bs))), function(z){
      if(as.numeric(y$chromStart)<as.numeric(z$transcript_end) && as.numeric(z$transcript_star)<as.numeric(y$chromEnd)){
        c(z,as.vector(y$name))
      }else {
        NA
      }
    })
    m[is.na(m)]<-NULL
    do.call(rbind,m)
  }))
}))
ARAlincRNA38En<-as.data.frame(ARAlincRNA38En)
# saveRDS(ARAlincRNA38En,file = 'cache/ARAlincRNA38En.rds') 
# saveRDS(ARAlincRNA38En,file = 'cache/ARAlincRNA38En2.rds') # use this
ARAlincRNA38En<-readRDS('cache/ARAlincRNA38En2.rds')
ARAlincRNA38En.filter<-filter(ARAlincRNA38En,gene_biotype%in%pv.lncRNA)
ARAlincRNA<-unique(do.call(c,ARAlincRNA38En.filter$ensembl_gene_id))
a<-filter(lnctfv1,v2=='AR')
dim(a)
a<-a%>%arrange(FDR)
