TPM <- function(countdata,exons_gene_lens = exons_gene_lens){
  intersect <- intersect(rownames(countdata),names(exons_gene_lens))
  countdata <- countdata[intersect,]
  exons_gene_lens <- unlist(exons_gene_lens[intersect])
  rpk <- countdata*1000/exons_gene_lens
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  return(tpm)
}
SetLassoData<-function(genes, rna.expr, metaMatrix,year=NULL){
  samples = colnames(rna.expr)[grep("-01",colnames(rna.expr))]
  exprDa=rna.expr[genes,samples]
  clinicalDa=metaMatrix[match(samples,metaMatrix$samples),]
  daysToDeath <- as.numeric(clinicalDa$A1_OS)
  
  vitalStatus <- clinicalDa$A2_Event=="Dead"
  #daysToDeath = daysToDeath/30
  if(is.null(year)){
    nonComplt <- is.na(daysToDeath)
  }else{
    nonComplt <- c(is.na(daysToDeath) | (daysToDeath>year*365))
  }
  
  daysToDeath = daysToDeath
  vitalStatus[nonComplt] <- FALSE
  print(length(daysToDeath))
  print(length(vitalStatus))
  print(dim(t(exprDa)))
  res = data.frame(Time = daysToDeath,Status=vitalStatus,t(exprDa),stringsAsFactors = F)
  return(res)
}
