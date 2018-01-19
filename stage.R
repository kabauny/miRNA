library(TCGAbiolinks)
library(DT)
library(DESeq)
#####################

clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
tumorStage = unique(clinical$tumor_stage)
barcodeList = c()

for(i in 1:length(tumorStage)){
  temp = clinical$tumor_stage == tumorStage[i]
  barcodeList[[i]] = clinical[temp, 'bcr_patient_barcode']
}


miRNA_barcode_Download = function(proj, dataType, bc){
  query <- GDCquery(project = proj,
                    data.category = "Transcriptome Profiling",
                    data.type = paste(dataType, "Expression Quantification"),
                    experimental.strategy = "miRNA-Seq",
                    barcode = bc
  )
  GDCdownload(query, method = "api", files.per.chunk = 10)
  DF <- GDCprepare(query)
  
  DF_name = DF$miRNA_ID
  DF$miRNA_ID <- NULL
  newDF = DF[c(TRUE, FALSE, FALSE)]
  rownames(newDF) = DF_name
  
  
  return(newDF)
  
}

query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = paste("miRNA", "Expression Quantification"),
                  experimental.strategy = "miRNA-Seq",
                  barcode = barcodeList[1]
)
