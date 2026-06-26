
#core components
ReadmiRNA = function(cancerType, sample, tumorStage){
  #creates dataframe given Cancer type and tumor stage. 
  #p x n
  #p name := miRNA name
  #n name := cancer_stage_index

  f = paste(parent, cancerType, sample, tumorStage, sep = "/")

  manif = paste(f, "MANIFEST.txt", sep = "/")
  manifest = read.table(file = manif, sep = "\t", header = TRUE)
  fileName = paste(f, toString(manifest$filename[1]), sep = "/")
  dft = read.table(file = fileName, sep = "\t", header = TRUE)
  df = data.frame(dft$read_count)
  colnames(df)[1] = paste(cancerType, tumorStage, "1", sep = "_")

  rownames(df) = dft$miRNA_ID
  
  for(i in 2:length(manifest$filename)){
    fileName = paste(f, toString(manifest$filename[i]), sep = "/")
    dft = read.table(file = fileName, sep = "\t", header = TRUE)
    tempSampleName = paste(cancerType, tumorStage, toString(i), sep = "_")
    df[tempSampleName] = dft$read_count
  }
  
  return(df)
}

ReadNormal = function(){

  f = paste(parent, "Normal", sep = "/")
  manif = paste(f, "MANIFEST.txt", sep = "/")
  manifest = read.table(manif, sep = "\t", header = TRUE)
  fileName = paste(f, toString(manifest$filename[1]), sep = "/")
  dft = read.table(fileName, sep = "\t", header = TRUE)
  df = data.frame(dft$read_count)
  
  colnames(df)[1] = paste("Normal", "1", sep = "_")
  
  rownames(df) = dft$miRNA_ID
  
  for(i in 2:length(manifest$filename)){
    fileName = paste(f, toString(manifest$filename[i]), sep = "/")
    dft = read.table(file = fileName, sep = "\t", header = TRUE)
    tempSampleName = paste("Normal", toString(i), sep = "_")
    df[tempSampleName] = dft$read_count
  }
  
  return(df)
}

add_Label = function(DF, n){
  #transpose dataframe to n x p
  #adds label column 
  DF = as.data.frame(t(DF))
  DF$label = rep(n, nrow(DF))
  return(DF)
}

join = function(a,b){
  #merges column together
  ab = merge(a,b,all = TRUE)
  return (ab)
}
##############################################
#cgts = centroid given cancer tyoe and tumor stage

diffgeneList = function(ct, ts){
  q = normal
  q = add_Label(q, "Normal")
  
  df = ReadmiRNA(cancerType = ct, sample = "PrimaryTumor", tumorStage = ts)
  p = add_Label(df, ct)

  pq = rbind(a = p, b = q)
  X = pq[,names(pq) != 'label']
  
  miRNA.data = list(x = t(X), y = pq$label, genenames = colnames(X), geneid = colnames(X))
  
  miRNA.train = pamr.train(miRNA.data)
  miRNA.cv = pamr.cv(miRNA.train, miRNA.data)
  deltaIndex = which.min(miRNA.cv$error)
  delta = miRNA.cv$threshold[deltaIndex]
  miRNA.confusion = pamr.confusion(miRNA.cv, delta)
  
  lg = pamr.listgenes(miRNA.train, miRNA.data, threshold=delta)
  lg = data.frame(lg)

  dl = df[row.names(df) %in% lg$id,]
  
  return(list(dl, df))  
}

acgts = function(ts){
  #adjusted centroid given cancer type and tumor stage
  temp = diffgeneList("LUAD", ts)
  ad = data.frame(temp[1])
  adeno = data.frame(temp[2])
  
  temp = diffgeneList("LUSC", ts)
  sq = data.frame(temp[1])
  squamous = data.frame(temp[2])
  
  ad2 = adeno[row.names(adeno) %in% row.names(sq),]
  sq2 = squamous[row.names(squamous) %in% row.names(ad),]
  
  atemp = ad[!(rownames(ad) %in% rownames(ad2)),]
  stemp = sq[!(rownames(sq) %in% rownames(sq2)),]
  
  ad3 = rbind(ad2, atemp)
  sq3 = rbind(sq2, stemp)
  
  ad4 = add_Label(ad3, "LUAD")
  sq4 = add_Label(sq3, "LUSC")

  adsq = join(a = ad4, b = sq4)

  X = adsq[,names(adsq) != 'label']
  miRNA.data = list(x = t(X), y = adsq$label, genenames = colnames(X), geneid = colnames(X))
  
  miRNA.train = pamr.train(miRNA.data)
  miRNA.cv = pamr.cv(miRNA.train, miRNA.data)
  deltaIndex = which.min(miRNA.cv$error)
  delta = miRNA.cv$threshold[deltaIndex]
  miRNA.confusion = pamr.confusion(miRNA.cv, delta)
  
  lg = pamr.listgenes(miRNA.train, miRNA.data, threshold=delta)
  lg = data.frame(lg)
  
  write.csv(lg$id, paste("/Users/zhongningchen/Data/miRNA/shrukenCentroid", ts, sep = "/"), quote = FALSE, row.names = FALSE)
  return(lg$id)

}

library(pamr)
parent = "/Users/zhongningchen/Data/miRNA"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")
normal = ReadNormal()

ts = tumorStageList[5]
#tempVariable = acgts(ts)

 


