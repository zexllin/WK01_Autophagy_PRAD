library(here)
library(GDCRNATools)
# getdata
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-PRAD',data.type  = 'RNAseq', write.meta = FALSE)

#Filter duplicated samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
# Removed 3 samples

#Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
#Removed 1 samples

#Merge RNAseq data 
datapath <- here("data/TCGA-PRAD/")
dat <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                   path=datapath, # the folder in which the data stored
                   organized = FALSE, # if the data are in separate folders
                   data.type = 'RNAseq')
save(dat,metaMatrix.RNA,file=here("results/TCGA-PRAD_RawData.rda"))

metadata = read.table(here("data/Clinical_BCR_XML.merge.txt"), header = T,sep = "\t")

saveRDS(metaMatrix.RNA,file=here("data/TCGA-PRAD_metaMatrixRNA.Rds"))
saveRDS(dat,file=here("data/TCGA-PRAD_RawData.Rds"))
