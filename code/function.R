#数据框矩阵去科学计数的小数点位数
setdigit<-function(idata,digits = 7){
  
  odata = data.frame(row.names = rownames(idata))
  for(colnm in colnames(idata)){
    col = idata[,colnm]
    message(colnm,": ",class(col),"\n")
    if(class(col)=="numeric" && sum(col%%1)!=0) {
      odata<- cbind(odata,signif(col,digits = digits))
    }else{
      odata<- cbind(odata,col)
    }
  }
  colnames(odata) = colnames(idata)
  return(odata)
}

# 热图
DE.limma.CSDN <- function(counts,group_list,filter = FALSE){
  # 指定TBI-Control，"-" 后面为对照
  dge <- DGEList(counts = counts)
  keep <- rowSums(cpm(dge) > 1) >= 0.5 * length(group_list)
  
  
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(counts)
  contrast.matrix<-makeContrasts("T-N",levels = design)
  
  if (filter == TRUE) {
    dge <- dge[keep, , keep.lib.sizes = TRUE]
  }
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v,design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)
  
  tempOutput = topTable(fit2, coef=1, n=Inf)
  deg = na.omit(tempOutput)
  
  # if (startsWith(rownames(deg)[1], "ENSG")) {
  #   degList <- gtype[match(rownames(deg), gtype$`Gene stable ID`),
  #                      ]
  #   degOutput <- data.frame(symbol = degList$`Gene name`,
  #                           group = degList$`Gene type`, deg)
  # 
  #   return(degOutput)
  # }
  # else {
  #   return(deg)
  # }
  if (startsWith(rownames(deg)[1], "ENSG")) {
    degList <- biotype[match(rownames(deg), biotype$ensemblID),
                       ]
    degOutput <- data.frame(symbol = degList$geneSymbol,
                            group = degList$group, deg)
    return(degOutput)
  }
  else {
    return(deg)
  }
}

#火山图
volcanoplot<-function(deg.plot,fc=2,pval=0.05,data.type=""){
  deg.plot<-cbind(rownames(deg.plot),deg.plot)
  colnames(deg.plot)<-c("ensemblID","symbol","group","logFC","AveExpr","t","PValue","adj.P.Val","B")
  threshold <- c()
  threshold[deg.plot$logFC>log(fc,2) & deg.plot$PValue<pval] <- "up"
  threshold[abs(deg.plot$logFC)<=log(fc,2) | deg.plot$PValue>=pval] <- "nosig"
  threshold[deg.plot$logFC < -log(fc,2) & deg.plot$PValue<pval] <- "down"
  deg.plot$threshold <- as.factor(threshold)
  
  p<-ggplot(deg.plot,aes(x=logFC,y=-log10(PValue)))+
    
    geom_point(aes(color=threshold))+
    scale_colour_manual( values = c("#DC143C","#808080","#00008B"))+
    geom_vline(xintercept = c(-log(fc,2),log(fc,2)),
               color='darkgreen', linetype=3,size=0.6) +
    geom_hline(yintercept = -log(pval,10), color='darkgreen',linetype=3,size=0.6)+
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(),
          panel.background = element_blank())+
    theme(legend.title=element_blank())+
    #theme(legend.key.size=unit(2,'cm'));
    #theme(legend.key.width=unit(2,'cm'));
    theme(legend.text = element_text(size = 30),legend.position = 'right')+ #raw 12
    theme(axis.text=element_text(size=30), # raw 13
          axis.title=element_text(size=30))
  pdf(paste0(data.type,"_Volcano.pdf"),width=15, height=8)
  print(p)
  dev.off()
  
  #print(p)
}

#热图
HeatmapPlot_t<-function(exprdata,Group = group,datatype='',save='pdf'){
  

  annotation_col <- data.frame('Sample' = Group)
  rownames(annotation_col) <- colnames(exprdata)
  #annotation_row<- data.frame(diff=ifelse(deg.sig[order(deg.sig$logFC),]$logFC > 1,"UP","DOWN" ))
  #rownames(annotation_row)<-rownames(degDa)
  pdfname=paste0(datatype,"DegHeatMap.pdf")
  tifname=paste0(datatype,"DegHeatMap.tif")
  #pdf(file=pdfname, width=12, height=8)
  p<-pheatmap(exprdata,
              #main = 'heatmap', # title
              scale = 'row', #column row none 
              annotation_col = annotation_col, 
              #annotation_row = annotation_row, 
              #legend_labels = NA,
              cluster_cols = FALSE, 
              color =colorRampPalette(c("blue", "white","red"))(100),
              #cluster_rows = FALSE,         
              #clustering_method = "complete", # complete average median 
              show_rownames = F, 
              show_colnames = F, 
              #gaps_row = 1169,
              fontsize = 24,
              border_color = NA,
              # cellwidth = 2.5,
              # cellheight = 1,
              angle_col=45)
  if(save=='tif'){
    tiff(tifname,width=200, height=150, units='mm',res=300)
    print(p)
    dev.off()
  }
  if(save=='pdf'){
    pdf(pdfname,width=15, height=15)
    print(p)
    dev.off()
  }
  #print(p)
  
}

DE_BarPlot <- function (deg, angle = 0, data.type) 
{
  deg = deg[which(!is.na(deg$group)),]
  if (data.type == "miRNAs") {
    down <- sum(deg$logFC < 0)
    up <- sum(deg$logFC > 0)
    d <- data.frame(geneClass = c("Up", "Down"), geneNums = c(up, 
                                                              down), Regulation = factor(c("Up-regulated", "Down-regulated"), 
                                                                                         levels = c("Up-regulated", "Down-regulated")))
    if (angle == 0) {
      ggplot(data = d, aes(x = d$geneClass, y = d$geneNums, 
                           fill = Regulation)) + geom_bar(stat = "identity") + 
        scale_x_discrete(limits = d$geneClass) + scale_fill_discrete(name = "") + 
        ylab("No. of Differentially Expressed miRNAs") + 
        xlab("") + theme_bw() + theme(axis.line = element_line(colour = "black"), 
                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                      panel.border = element_rect(colour = "white"), 
                                      panel.background = element_blank(), axis.text.x = element_text(angle = angle, 
                                                                                                     size = 10), axis.text.y = element_text(size = 10))
    }
    else {
      ggplot(data = d, aes(x = d$geneClass, y = d$geneNums, 
                           fill = Regulation)) + geom_bar(stat = "identity") + 
        scale_x_discrete(limits = d$geneClass) + scale_fill_discrete(name = "") + 
        ylab("No. of Differentially Expressed miRNAs") + 
        xlab("") + theme_bw() + theme(axis.line = element_line(colour = "black"), 
                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                      panel.border = element_rect(colour = "white"), 
                                      panel.background = element_blank(), axis.text.x = element_text(angle = angle, 
                                                                                                     vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
    }
  }
  else if (data.type == "RNAseq") {
    gr <- list()
    for (tp in unique(deg$group)) {
      gr[[tp]][["all"]] <- sum(deg$group == tp)
      gr[[tp]]["up"] <- sum(deg$group == tp & deg$logFC > 
                              0)
      gr[[tp]]["down"] <- sum(deg$group == tp & deg$logFC < 
                                0)
    }
    names(gr)[names(gr) == "protein_coding"] <- "Protein coding"
    names(gr)[names(gr) == "long_non_coding"] <- "Long non-coding"
    names(gr)[names(gr) == "pseudogene"] <- "Pseudogene"
    names(gr)[names(gr) == "ncRNA"] <- "Other ncRNA"
    o <- order(unlist(lapply(gr, function(v) v[["all"]])), 
               decreasing = TRUE)
    d <- data.frame(geneClass = rep(names(gr)[o], 2), geneNums = c(unlist(lapply(gr, 
                                                                                 function(v) v[["up"]]))[o], unlist(lapply(gr, function(v) v[["down"]]))[o]), 
                    Regulation = factor(rep(c("Up-regulated", "Down-regulated"), 
                                            each = length(gr)), levels = c("Up-regulated", 
                                                                           "Down-regulated")))
    if (angle == 0) {
      ggplot(data = d, aes(x = d$geneClass, y = d$geneNums, 
                           fill = Regulation)) + geom_bar(stat = "identity") + 
        scale_x_discrete(limits = d$geneClass[seq_len(nrow(d)/2)]) + 
        scale_fill_discrete(name = "") + ylab("No. of Differentially Expressed Genes") + 
        xlab("") + theme_bw() + theme(axis.line = element_line(colour = "black"), 
                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                      panel.border = element_rect(colour = "white"), 
                                      panel.background = element_blank(), axis.text.x = element_text(angle = angle, 
                                                                                                     size = 14), axis.text.y = element_text(size = 14))
    }
    else {
      ggplot(data = d, aes(x = d$geneClass, y = d$geneNums, 
                           fill = Regulation)) + geom_bar(stat = "identity") + 
        scale_x_discrete(limits = d$geneClass[seq_len(nrow(d)/2)]) + 
        scale_fill_discrete(name = "") + ylab("No. of Differentially Expressed Genes") + 
        xlab("") + theme_bw() + theme(axis.line = element_line(colour = "black"), 
                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                      panel.border = element_rect(colour = "white"), 
                                      panel.background = element_blank(), axis.text.x = element_text(angle = angle, 
                                                                                                     vjust = 1, hjust = 1, size = 14), axis.text.y = element_text(size = 14)) + 
        theme(axis.title = element_text(size = 16), legend.text = element_text(size = 14))
    }
  }
}
#cox 生存
coxphTestFun <- function(genes, rna.expr, metaMatrix) {
  
  
  
  samples = intersect(colnames(rna.expr), paste0(metaMatrix$A0_Samples,"-01"))
  exprDa=rna.expr[genes,samples]
  
  clinicalDa=metaMatrix[match(samples,paste0(metaMatrix$A0_Samples,"-01")),]
  daysToDeath <- as.numeric(clinicalDa$A1_OS)
  
  vitalStatus <- clinicalDa$A2_Event=="Dead"
  #daysToDeath = daysToDeath/30

  
  coxphDEGs <- c()
  for (i in seq_len(nrow(exprDa))) {
    DEG <- unlist(exprDa[i,])
    coxtest <- coxph(Surv(daysToDeath, vitalStatus) ~ DEG)
    
    summcph <- summary(coxtest)
    coeffs <- c(summcph$coefficients[,1:2], summcph$conf.int[,3:4], 
                summcph$coefficients[,5])
    coxphDEGs <- rbind(coxphDEGs, coeffs)
    
  }
  rownames(coxphDEGs) <- rownames(exprDa)
  
  colnames(coxphDEGs) <- c('coef','HR','lower95','upper95','pValue')
  coxphDEGs <- data.frame(symbol=biotype$geneSymbol[match(rownames(exprDa), biotype$ensemblID)],
                          coxphDEGs)
  #coxphDEGs$FDR <- p.adjust(coxphDEGs$pValue, method='fdr')
  
  #o <- order(coxphDEGs$pValue)
  #coxphDEGs <- coxphDEGs[o,]
  
  return (coxphDEGs)
}

