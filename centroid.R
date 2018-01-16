library(TCGAbiolinks)
library(DT)
library(DESeq)
library(pamr)

miRNA_Download = function(proj, dataType){
  query <- GDCquery(project = proj,
                               data.category = "Transcriptome Profiling",
                               data.type = paste(dataType, "Expression Quantification"),
                               experimental.strategy = "miRNA-Seq"
  )
  GDCdownload(query, method = "api", files.per.chunk = 10)
  DF <- GDCprepare(query)
  
  DF_name = DF$miRNA_ID
  DF$miRNA_ID <- NULL
  newDF = DF[c(TRUE, FALSE, FALSE)]
  rownames(newDF) = DF_name
  
  
  return(newDF)
  
}

add_Label = function(DF, n){
  
  DF = as.data.frame(t(DF))
  DF$label = rep(n, nrow(DF))
  return(DF)
}

prepare = function(p, d){
  dd = miRNA_Download(proj = p, dataType = d)
  return(add_Label(dd, n = p))
}

join = function(a,b){
  ab = merge(a,b,all = TRUE)
  return (ab)
}

LUAD = prepare(proj = "TCGA-LUAD", dataType = "miRNA")
LUSC = prepare(proj = "TCGA-LUSC", dataType = "miRNA")
DC = join(LUAD, LUSC)
X = DC[,names(DC) != 'label']
miRNA.data = list(x = t(X), y = DC$label, genenames = colnames(X))
pamr.menu(miRNA.data)
