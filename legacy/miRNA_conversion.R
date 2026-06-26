library(biomaRt)
parent = "/Users/zhongningchen/Data/miRNA/shrukenCentroid"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")
mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

MirtoEnsg = function(tsi){
  mir_df = read.table(paste(parent, tumorStageList[tsi], sep = "/"), sep = "\t")
  mir = mir_df$V1
  
  desired_attributes = c("mirbase_id", "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")

  
  ensg = getBM(attributes= desired_attributes, filters = 'mirbase_id', 
                   values = mir, mart = mart)
  return(ensg)
  
}

tsi = 1
ts = tumorStageList[tsi]
df2 = MirtoEnsg(tsi)

write.csv(df, paste("/Users/zhongningchen/Data/miRNA/geneList", ts, sep = "/"), 
          quote = FALSE, row.names = FALSE)
