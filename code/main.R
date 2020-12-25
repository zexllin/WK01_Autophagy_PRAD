rm(list = ls())
library(here)
source(here("code/function.R"))
source(here("code/packages.R"))
path = here("results/")
gtype = read.csv(here("material/mart_export.txt"),sep = ",",check.names = F)

# get data
load(paste0(path,"TCGA-PRAD_RawData.rda"))
metadata = read.table(here("data/Clinical_BCR_XML.merge.txt"), header = T,sep = "\t",stringsAsFactors = F)
autophagygenes = read.table(here("data/Autophagy genes.txt"),header = F,stringsAsFactors = F)$V1
autophagygenes = unique(autophagygenes)
# DE analysis
group_list = sapply(colnames(dat),function(x){
  if(grepl("-01$",x)){
    return("T")
  }else if(grepl("-11$",x)){
    return("N")
  }})

dir.create(paste0(path,"2.DeAnalysis"))
# DEGAll <- gdcDEAnalysis(counts = dat,group = group_list,comparison = 'T-N',filter = T)
# deall <- DEGAll[which(DEGAll$PValue<0.05&abs(DEGAll$logFC)>1),]

DEGAll <- DE.limma.CSDN(counts = dat,group_list = group_list,filter = F)
deall <- DEGAll[which(DEGAll$P.Value<0.05&abs(DEGAll$logFC)>1),]
save(DEGAll,file = paste0(path,"DEGAll.rda"))

rnaExpr <- gdcVoomNormalization(counts = dat, filter = FALSE)
volcanoplot(DEGAll,data.type = here::here("results/2.DeAnalysis/DeAnalysis"))

DE_BarPlot(deg = deall, angle = 45, data.type = 'RNAseq')

gT<-which(group_list=="T")
gN<-which(group_list=="N")
group_plot<-rep(c("Tumor","Normal"),c(length(gT),length(gN)))
HeatmapPlot_t(exprdata = rnaExpr[rownames(deall),c(gT,gN)],Group = group_plot,
              datatype = here("results/2.DeAnalysis/"),save = "pdf")


# corralation
lncs = deall[which(deall$group=='long_non_coding'),]
autophagygenes = bitr(autophagygenes,fromType = 'SYMBOL',toType = c('ENSEMBL','ENTREZID'),OrgDb = org.Hs.eg.db)

correlation = list()
autogenes = intersect(autophagygenes$ENSEMBL,rownames(deall))
correlation$cor = matrix(,length(autogenes),nrow(lncs))
rownames(correlation$cor) = autogenes
colnames(correlation$cor) = rownames(lncs)
correlation$P = correlation$cor
correlation$merge = matrix(,0,4)
colnames(correlation$merge)<-c("AutophagyGenes","lncRNA","cor","Pvalue")
for (agene in autogenes) {
  
  for(lnc in rownames(lncs)){
    corr = cor.test(rnaExpr[agene,],rnaExpr[lnc,])
    cor = as.numeric(corr$estimate)
    pvalue = as.numeric(corr$p.value)
    correlation$cor[agene,lnc]=cor
    correlation$P[agene,lnc]=pvalue
    correlation$merge=rbind(correlation$merge,c(agene,lnc,cor,pvalue))
  }
}

cor_sig = as.data.frame(correlation$merge,stringsAsFactors = F)
cor_sig$cor=as.numeric(cor_sig$cor)
cor_sig$Pvalue=as.numeric(cor_sig$Pvalue)
cor_sig = cor_sig[which(cor_sig$cor>0.7&cor_sig$Pvalue<0.05),]



# step 5
# (1)
tp1 <- gdcSurvivalAnalysis(gene     = unique(cor_sig$lncRNA), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)
survOutput = coxphTestFun(genes = unique(cor_sig$lncRNA),rna.expr = rnaExpr,metaMatrix = metadata)
# (2)

set.seed(1011)
samples = intersect(colnames(rnaExpr), paste0(metadata$A0_Samples,"-01"))
samp_test = [sample(length(samples),length(samples)/4)
normdata = rnaExpr[unique(cor_sig$lncRNA),samples]
lable = metadata[match(samples,metadata$A0_Samples),c("A1_OS","A2_Event")]


cv.fit.train <- cv.glmnet(norm_trainData1,y,family="binomial",type.measure="deviance",nfolds = 5,keep=T)
# 提取模型
fit.train <- cv.fit.train$glmnet.fit
