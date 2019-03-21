source('./R/lib.R')
biomart <- pf.get.biomart()
diff<-pf.get.diff()
lncRNA<-pf.get.diff('lncRNA')
diff.anno <- filter(biomart,ensembl_gene_id%in%diff$GeneID)%>%
          mutate(start=transcript_start-5000,end=transcript_start)%>%
          select(ensembl_gene_id,chr=chromosome_name,start,end)
  
lnc.anno <- filter(biomart,ensembl_gene_id%in%lncRNA$GeneID)%>%
  select(ensembl_gene_id,chr=chromosome_name,tss=transcript_start)

i<-0
total<-nrow(lnc.anno)
lnc.gene.raw<-data.frame()
a<-lapply(split(lnc.anno, seq(nrow(lnc.anno))), function(lnc){
  i<<-i+1
  print(paste0(i,'/',total))
  diff.tmp<-filter(diff.anno,chr==lnc$chr)
  sapply(split(diff.tmp, seq(nrow(diff.tmp))), function(gene){
    if(lnc$tss>=gene$start&&lnc$tss<gene$end) {
      lnc.gene.raw<<-rbind(lnc.gene.raw,data.frame(lncRNA=lnc$ensembl_gene_id,
                                          gene=gene$ensembl_gene_id,
                                          dist=gene$end-lnc$tss))
    }
  })
})
saveRDS(lnc.gene.raw, 'support/lnc.cis.origin.rds')
lnc.cis<-pf.get.lnc.cis(up=1000)
#############################################
mutate(lnc.cis,
       lncRNA=pf.ensembl2symbol(lncRNA),
       gene=pf.ensembl2symbol(gene),
       chr=pf.ensembl2chr(lnc.cis$lncRNA))
