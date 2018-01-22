library(pamr)
parent = "/Users/zhongningchen/Data"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iiia", "stage_iiib", "stage_iv")

ReadmiRNA = function(cancerType, sample, tumorStage){
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

add_Label = function(DF, n){
  
  DF = as.data.frame(t(DF))
  DF$label = rep(n, nrow(DF))
  return(DF)
}

join = function(a,b){
  ab = merge(a,b,all = TRUE)
  return (ab)
}

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
  q = ReadmiRNA(cancerType = ct, sample = "Normal")
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

geneList = diffgeneList("LUAD", tumorStageList[1])
LUAD.stage_ia = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[1])


a = LUAD.stage_ia[row.names(LUAD.stage_ia) == geneList$id]
