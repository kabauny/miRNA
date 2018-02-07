library(biomaRt)
tsi = 1
parent = "/Users/zhongningchen/Data/mRNA/shrukenCentroid"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")
mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

ensgtoGeneName = function(tsi){
  ensg_df = read.table(paste(parent, tumorStageList[tsi], sep = "/"), sep = "\t")
  desired_attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "mirbase_id")
  
  ensg = gsub("\\..*","",ensg_df$V1)
  ensgv = ensg_df$V1
  
  geneName = getBM(attributes= desired_attributes, filters = 'ensembl_gene_id', 
                   values = ensg, mart = mart)
  return(geneName)
  
}
tsi = 1
df = ensgtoGeneName(tsi)
ts = tumorStageList[tsi]
write.csv(df, paste("/Users/zhongningchen/Data/mRNA/geneList", ts, sep = "/"), 
          quote = FALSE, row.names = FALSE)


