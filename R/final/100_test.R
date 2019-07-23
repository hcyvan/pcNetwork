ac<-pf.get.ac()

hub<-read.csv('./reports/thesis/13059_2017_1266_MOESM3_ESM.csv')

sort(intersect(hub$TF.hubs,ac$ID))


##############################################
zfpkm<-pf.get.zfpkm()
samples<-unique(sapply(str_split(colnames(zfpkm),'\\.'),function(x){paste(x[1:3],collapse = '-')}))

clinical <- read.delim('./data/tcga-prad/nationwidechildrens.org_clinical_patient_prad.txt', stringsAsFactors=FALSE)[-c(1,2),]
clinical <- filter(clinical, bcr_patient_barcode%in%samples)
