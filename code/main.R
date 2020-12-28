rm(list = ls())
library(here)
source(here("code/function.R"))
source(here("code/packages.R"))
source(here("code/method.R"))
dir.create(here("results/"))
path = here("results/")


## get data
dat = readRDS(here("data/TCGA-PRAD_RawData.Rds"))
metaMatrix.RNA = readRDS(here("data/TCGA-PRAD_metaMatrixRNA.Rds"))
metadata = read.table(here("data/Clinical_BCR_XML.merge.txt"), header = T,sep = "\t",stringsAsFactors = F)
exons_gene_lens = readRDS( here("material/exons_gene_lens.Rds"))
rnaExpr.tpm = TPM(countdata = dat,exons_gene_lens = exons_gene_lens)

object = list()
object[["assays"]][["count"]]=dat
object[["assays"]][["TPM"]]=rnaExpr.tpm
meta.merge = data.frame(samples=colnames(dat),
                        samp.name=substr(colnames(dat),1,12),
                        stringsAsFactors = F)
meta.merge = cbind(meta.merge,metadata[match(meta.merge$samp.name,metadata$A0_Samples),])
object[["meta.data"]]=meta.merge
object[["info"]]=meta.merge[,c("samples","samp.name")]


## get autophagygenes
autophagygenes = read.table(here("data/Autophagy genes.txt"),header = F,stringsAsFactors = F)$V1
autophagygenes = unique(autophagygenes)
autophagygenes = bitr(autophagygenes,fromType = 'SYMBOL',toType = c('ENSEMBL','ENTREZID'),OrgDb = org.Hs.eg.db)
object$autophagygenes = autophagygenes

## group for DE analysis 
dir.create(paste0(path,"2.DeAnalysis"))

## split data set
set.seed(10110)
index <- sample(ncol(dat),size = round(ncol(dat)*0.7,0),replace=F)

expr.train = dat[,index]
expr.test = dat[,-index]
object$info$Split = rep("Test",ncol(dat))  
object$info$Split[index] = "Train"
object$info$Group = ifelse(substr(object$info$samples,14,15)=="01","Tumor","Normal")

## DE analysis
index = c(which(object$info$Split=="Train" & object$info$Group=="Tumor"),
          which(object$info$Split=="Train" & object$info$Group=="Normal"))

De.count = object$assays$count[,index]
De.Group = object$info$Group[index]

if(F){
  design <- model.matrix(~ Group+0)
  #design <- model.matrix(~factor(group.1)+0)
  colnames(design) <- c("Normal","Tumor")
  rownames(design)=colnames(logdata)
  
  fit <- lmFit(logdata, design)
  cont.matrix <- makeContrasts(Tumor-Normal, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  DEGAll <- topTable(fit2, adjust="fdr", sort.by="B", number=length(fit2$coefficients))
}

res.de <- DE.limma.CSDN(counts = De.count,compare = "Tumor-Normal",group_list = De.Group,filter = F)
DEGAll = res.de$DE
deall <- DEGAll[which(DEGAll$P.Value<0.05&abs(DEGAll$logFC)>1),]
object[['DeAnalysis']][["count"]] = De.count
object[['DeAnalysis']][["Group"]] = De.Group
object[['DeAnalysis']][["DEGAll"]] = DEGAll
object[['DeAnalysis']][["norm_voom"]] = res.de$norm_voom$E
save(DEGAll,file = paste0(path,"DEGAll.rda"))

# Volcanoplot for De Analysis
volcanoplot(object$DeAnalysis$DEGAll,data.type = here::here("results/2.DeAnalysis/DeAnalysis"))
#BarPlot for De Analysis
DE_BarPlot(deg = deall, angle = 45, data.type = 'RNAseq')

# heatmap for De Analysis
HeatmapPlot_t(exprdata = object$DeAnalysis$norm_voom[rownames(deall),],Group = object$DeAnalysis$Group,
              datatype = here("results/2.DeAnalysis/"),save = "pdf")


lncs = deall[which(deall$group=='long_non_coding'),]
# lncs = biotype[which(biotype$group=="long_non_coding"),]
# lncs = lncs[which(lncs$ensemblID %in% intersect(rownames(rnaExpr),lncs$ensemblID)),]
# rownames(lncs) = lncs$ensemblID

# ARg-lncRNAs corralation 
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
    corr = cor.test(object$DeAnalysis$norm_voom[agene,],object$DeAnalysis$norm_voom[lnc,])
    cor = as.numeric(corr$estimate)
    pvalue = as.numeric(corr$p.value)
    correlation$cor[agene,lnc]=cor
    correlation$P[agene,lnc]=pvalue
    correlation$merge=rbind(correlation$merge,c(agene,lnc,cor,pvalue))
  }
}

cor_merge = as.data.frame(correlation$merge,stringsAsFactors = F)
cor_merge$cor=as.numeric(cor_merge$cor)
cor_merge$Pvalue=as.numeric(cor_merge$Pvalue)
correlation$merge = cor_merge
cor_sig = cor_merge[which(cor_merge$cor>0.5&cor_merge$Pvalue<0.05),]


object[["cor"]][["lncs"]] = lncs
object[["cor"]][["correlation"]] = correlation

## step 5
## (1) unVarcox
if(F){
  rnaExpr = gdcVoomNormalization(counts =object$DeAnalysis$count,filter = F )
  tp1 <- gdcSurvivalAnalysis(gene     = unique(cor_sig$lncRNA), 
                                    method   = 'coxph', 
                                    rna.expr = rnaExpr, 
                                    metadata = metaMatrix.RNA)
  ## KM unVar surv
  survOutput = surv_kmFun(genes = unique(cor_sig$lncRNA),
                          rna.expr = object$DeAnalysis$norm_voom$E,
                          metaMatrix =object$meta.data)
}## 测试GDCRNATools coxph生存分析和单因素KM生存分析

## unVar cox
survOutput=coxphTestFun(genes = unique(cor_sig$lncRNA),
                        rna.expr = object$DeAnalysis$norm_voom,
                        metaMatrix =object$meta.data)

survOutput_sig = survOutput[which(as.numeric(survOutput$pValue)<0.1),]
genes = rownames(survOutput_sig)
surv_KMPlot(gene = genes[1],rna.expr = rnaExpr,metadata = metadata)

## (2) lasso
# multVarCoxph(genes = genes, rna.expr = expr.train, metaMatrix = object$meta.data)

lasso.train = SetLassoData(genes = genes,rna.expr = expr.train,metaMatrix = object$meta.data)
lasso.test = SetLassoData(genes = genes,rna.expr = expr.test,metaMatrix = object$meta.data)

x=as.matrix(lasso.train[,3:ncol(lasso.train)])
y=data.matrix(Surv(lasso.train$Time,lasso.train$Status))
# write.csv(x,file = paste0('Result/5.lasso/','x.csv'))
# write.csv(y,file = paste0('Result/5.lasso/','y.csv'))
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
fit.train <- cvfit$glmnet.fit

pdf("图14. LASSO模型.pdf")
par(mgp = c(2.5,1,0),mar=c(4,4.5,1,1),mai=c(1,1,1,1))
plot(fit.train, xvar="lambda",cex.lab = 2)+
  abline(v = c(log(cvfit$lambda.min), log(cvfit$lambda.1se)),lty=2)+
  text(x = log(cvfit$lambda.min),y = 0.0,
       paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1.2,adj=0.9)+
  text(x = log(cvfit$lambda.1se),y = 0.2,
       paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1.2,adj=0.9)
dev.off()

pdf("图15.LASSO Deviance图.pdf")
par(mgp = c(3,1,0),mar=c(4,4.8,1,1),mai=c(1,1,1,1))
plot(cvfit,cex.lab = 2)+  
  text(x = log(cvfit$lambda.min),y = 0.0,
       paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1.2,adj=0.9)+
  text(x = log(cvfit$lambda.1se),y = 0.2,
       paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1.2,adj=0.9)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(cvfit, s = cvfit$lambda.min)
index <- which(as.numeric(coef) != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("Status","Time",lassoGene)
cox.data=lasso.train[,lassoGene]
lassoSigExp=cbind(id=row.names(cox.data),cox.data)

# fml <- as.formula(paste0('Surv(Time,Status)~',paste0(colnames(cox.data)[-c(1:2)],collapse = '+')))
# f <- coxph(fml, data=cox.data,id = rownames(cox.data))
# riskScore=predict(f,type="risk",newdata=cox.data)


##-------------------------------------------------------------------------------------------------------------
save(object,file = here("cache/object.rda"))
# object[["info"]]=meta.merge[,c("samples","samp.name","A18_Sex","age_at_initial_pathologic_diagnosis","A1_OS",
#                                "A2_Event","A3_T","A4_N","A5_M")]
# colnames(object$info) = c("samples","samp.name","Sex","Age","OS","Event","T","N","M")