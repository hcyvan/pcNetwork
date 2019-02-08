fpkm <- pf.get.fpkm()

zfpkm <- t(apply(fpkm[apply(fpkm, 1, sum)!=0,],1,function(x){
  if(min(x)==0){x[x==0] <- min(x[x!=0])/1e10}
  scale(log2(x))
}))
colnames(zfpkm) <- colnames(fpkm)

saveRDS(zfpkm, './data/prad.rna.zfpkm.rds')
