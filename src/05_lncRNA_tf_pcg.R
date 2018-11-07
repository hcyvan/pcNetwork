library(dplyr)

source('./lib/globals.R')
source('./lib/helpers.R')

bg1 <- read.delim('./data/ENCODE_TF_ChIP-seq_2015.handle.txt', sep = ',',header = TRUE, stringsAsFactors = FALSE)
bg2 <- read.delim('./data/ChEA_2016.handle.txt', sep = '\t',header = TRUE, stringsAsFactors = FALSE)

bg1 <- bg1[,c('tf','symbol')]
bg3 <- rbind(bg1,bg2)

bg1 <- distinct(bg1)
bg2 <- distinct(bg2)
bg3 <- distinct(bg3)

tfs <- read.delim('./reports/tf_enricher/all.csv',sep = ',')
tfs$ID

starbase <- read.delim('./data/lncRNA_rbp.txt')
unique(starbase$RBP)

intersect(tfs$ID, unique(starbase$RBP))
intersect(unique(bg3$tf), unique(starbase$RBP))

modules <- read.delim('./data/diff.qlf.2877.wgcna.color.csv', sep = ',',header = TRUE, stringsAsFactors = FALSE)

lncRNAs <- (modules %>% filter(type%in%config$lncRNA))$id
pcgs <- (modules %>% filter(type%in%config$PCGs))$id




for(c in unique(modules$color)){
  cat(paste0(c, ': '))
  cat(intersect(enricher.bg3$tfs[[c]],tf.base))
  cat('\n')
}


lnRNA2tf <- read.

tf.base <- c('ACIN1', 'ADAR', 'AIFM1', 'ALKBH5', 'ALYREF', 'AUH', 'BCCIP', 'BUD13', 'CAPRIN1', 'CBX7', 'CELF2', 'CNBP', 'CPSF6', 'CSTF2T', 'DDX3X', 'DDX42', 'DDX54', 'DGCR8', 'DHX9', 'DICER1', 'DIS3L2', 'DKC1', 'EIF3A', 'EIF3B', 'EIF3D', 'EIF3G', 'EIF4A1', 'EIF4A3', 'EIF4G1', 'EIF4G2', 'ELAVL1', 'ELAVL3', 'EWSR1', 'FAM120A', 'FBL', 'FKBP4', 'FMR1', 'FTO', 'FUS', 'FXR1', 'FXR2', 'GNL3', 'GTF2F1', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPD', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1', 'IGF2BP1', 'IGF2BP2', 'IGF2BP3', 'ILF3', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'KHSRP', 'LARP4B', 'LARP7', 'LIN28', 'LIN28A', 'LIN28B', 'LSM11', 'MBNL1', 'MBNL2', 'METTL14', 'METTL3', 'MOV10', 'MSI1', 'MSI2', 'NCBP3', 'NONO', 'NOP56', 'NOP58', 'NPM1', 'NUMA1', 'PAPD5', 'PCBP2', 'PRPF8', 'PTBP1', 'PUM1', 'PUM2', 'QKI', 'RANGAP1', 'RBFOX2', 'RBM10', 'RBM22', 'RBM27', 'RBM39', 'RBM47', 'RBM5', 'RBM6', 'RC3H1', 'RNF219', 'RTCB', 'SAFB2', 'SBDS', 'SF3A3', 'SF3B4', 'SLBP', 'SLTM', 'SMNDC1', 'SND1', 'SRSF1', 'SRSF10', 'SRSF3', 'SRSF7', 'SRSF9', 'TAF15', 'TARBP2', 'TARDBP', 'TIA1', 'TIAL1', 'TNRC6A', 'TRA2A', 'TROVE2', 'U2AF1', 'U2AF2', 'UPF1', 'VIM', 'WTAP', 'XRN2', 'YTHDC1', 'YTHDC2', 'YTHDF1', 'YTHDF2', 'YWHAG', 'ZC3H7B', 'ZFP36', 'ZNF184')
xrn2 <- read.delim('./data/starBaseV3_hg19_CLIP-seq_XRN2_all.csv', sep = ',')
fus <- read.delim('./data/starBaseV3_hg19_CLIP-seq_FUS_all.txt', comment.char = "#")

intersect(modules$id, xrn2$geneID) -> idx
modules %>% filter(id%in%idx)

intersect(modules$id, fus$geneID) -> idx
modules %>% filter(id%in%idx)