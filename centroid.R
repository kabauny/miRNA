library(pamr)
parent = "/Users/zhongningchen/Data"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iiia", "stage_iiib", "stage_iv")
###############################################
#core components
ReadmiRNA = function(cancerType, sample, tumorStage){
  #creates dataframe given Cancer type and tumor stage. 
  #p x n
  #p name := miRNA name
  #n name := cancer_stage_index
  if (sample != "Normal"){
    f = paste(parent, cancerType, sample, tumorStage, sep = "/")
  } else{
    f = paste(parent, cancerType, sample, sep = "/")
  }
  
  manif = paste(f, "MANIFEST.txt", sep = "/")
  manifest = read.table(file = manif, sep = "\t", header = TRUE)
  fileName = paste(f, toString(manifest$filename[1]), sep = "/")
  dft = read.table(file = fileName, sep = "\t", header = TRUE)
  df = data.frame(dft$read_count)
  
  if (sample != "Normal"){
    colnames(df)[1] = paste(cancerType, tumorStage, "1", sep = "_")
  } else{
    colnames(df)[1] = paste(cancerType, "Normal", "1", sep = "_")
  }
  
  rownames(df) = dft$miRNA_ID
  
  for(i in 2:length(manifest$filename)){
    fileName = paste(f, toString(manifest$filename[i]), sep = "/")
    dft = read.table(file = fileName, sep = "\t", header = TRUE)
    if (sample != "Normal"){
      tempSampleName = paste(cancerType, tumorStage, toString(i), sep = "_")
    } else{
      tempSampleName = paste(cancerType, "Normal", toString(i), sep = "_")
    }
    
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
cgts = function(ts){
  adeno = ReadmiRNA(cancerType = "LUAD", sample = "PrimaryTumor", tumorStage = ts)
  adeno = add_Label(adeno, "LUAD")
  sq = ReadmiRNA(cancerType = "LUSC", sample = "PrimaryTumor", tumorStage = ts)
  sq = add_Label(sq, "LUSC")
  
  adsq = join(a = adeno, b = sq)
  X = adsq[,names(adsq) != 'label']
  miRNA.data = list(x = t(X), y = adsq$label, genenames = colnames(X), geneid = colnames(X))
  return(miRNA.data)
}

normcgts = function(ct, ts){
  p = ReadmiRNA(cancerType = ct, sample = "PrimaryTumor", tumorStage = ts)
  p = add_Label(p, ct)
  #q = ReadmiRNA(cancerType = ct, sample = "Normal")
  q = ReadNormal()
  q = add_Label(q, "Normal")
  
  pq = join(a = p, b = q)
  X = pq[,names(pq) != 'label']
  miRNA.data = list(x = t(X), y = pq$label, genenames = colnames(X), geneid = colnames(X))
  return(miRNA.data)
}

diffgeneList = function(ct, ts){
  miRNA = normcgts(ct, ts)
  miRNA.train = pamr.train(miRNA)
  deltaIndex = which.min(miRNA.train$errors)
  delta = miRNA.train$threshold[deltaIndex]
  lg = pamr.listgenes(miRNA.train, miRNA, threshold=delta)
  lg = data.frame(lg)
  return(lg)
}

adj = function(ct, ts){
  b = diffgeneList(ct, ts)
  df = ReadmiRNA(ct, "PrimaryTumor", ts)
  return(df[row.names(df) %in% b$id,])   
}

acgts = function(ts){
  #adjusted centroid given cancer type and tumor stage
  adeno = ReadmiRNA(cancerType = "LUAD", sample = "PrimaryTumor", tumorStage = ts)
  squamous = ReadmiRNA(cancerType = "LUSC", sample = "PrimaryTumor", tumorStage = ts)
  
  ad = adj("LUAD", ts)
  sq = adj("LUSC", ts)
  
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
  
  print(length(ad3))
  print(length(sq3))
  return(miRNA.data)
}

ac = acgts(tumorStageList[6])
pamr.menu(ac)
