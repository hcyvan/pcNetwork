setwd('../download-firehouse/gdac.broadinstitute.org_PRAD.Methylation_Preprocess.Level_3.2016012800.0.0/')

raw <- lapply(dir(),data.table::fread)
sapply(raw,dim)


############################
library(TCGAbiolinks)
query.met <- GDCquery(project = "TCGA-PRAD", 
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)
met <- GDCprepare(query = query.met,
                  save = TRUE, 
                  save.filename = "DNAmethylation_PRAD.rda",
                  summarizedExperiment = TRUE)