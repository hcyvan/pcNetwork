library(dplyr)
library(stringr)
library(ggpubr)
library(cowplot)
library(survival)
library(survminer)

source('./R/lib.R')

clinical1.raw <- read.delim('./data/tcga-prad/nationwidechildrens.org_clinical_patient_prad.txt', stringsAsFactors=FALSE)[-c(1,2),]
clin.df <- dplyr::select(clinical1.raw,
                         case=bcr_patient_barcode,
                         gleason=gleason_score,
                         new_tumor=days_to_biochemical_recurrence_first,
                         #new_tumor=new_tumor_event_dx_indicator,
                         last=last_contact_days_to)%>%
  mutate(new_tumor=suppressWarnings(as.numeric(new_tumor)),
         last=suppressWarnings(as.numeric(last)))
df<-data.frame(case=as.vector(clin.df$case),
               time=ifelse(is.na(clin.df$new_tumor), clin.df$last,clin.df$new_tumor),
               status=ifelse(is.na(clin.df$new_tumor), 1,2),
               stringsAsFactors = FALSE)

clinical2.raw <- read.delim('./data/tcga-prad/clinical/clinical.tsv', stringsAsFactors=FALSE)
clin.o <- dplyr::select(clinical2.raw,
                        case=submitter_id,
                        death=days_to_death,
                        last=days_to_last_follow_up)%>%
  mutate(death=suppressWarnings(as.numeric(death)),
         last=suppressWarnings(as.numeric(last)))

overall<-data.frame(case=as.vector(clin.o$case),
                    time=ifelse(is.na(clin.o$death), clin.o$last,clin.o$death),
                    status=ifelse(is.na(clin.o$death), 1,2),
                    stringsAsFactors = FALSE)

saveRDS(df,'./support/clinical.disesefree.rds')
saveRDS(overall,'./support/clinical.overall.rds')
