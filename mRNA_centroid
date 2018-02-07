
###############################################
#core components
mod_filename = function(dir){
  
  df = read.table(dir, sep = "\t", header = TRUE)
  f = paste(df$id, df$filename, sep = "/")
  return (data.frame(f))

  }

ReadmiRNA = function(cancerType, tumorStage){
  #creates dataframe given Cancer type and tumor stage. 
  #p x n
  #p name := miRNA name
  #n name := cancer_stage_index
  f = paste(parent, cancerType, sep = "/")
  
  manif = paste(f, "Manifest", paste(tumorStage, "txt", sep = "."), sep = "/")
  manifest_subset = mod_filename(manif)
  colnames(manifest_subset) = "filename"
  
  directory = paste(f, "Full", sep = "/")
  fn = paste(directory, toString(manifest_subset$filename[1]), sep = "/")
  dft = read.table(gzfile(fn))
  
  dft = head(dft, -5)
  
  df = data.frame(dft$V2)
  rownames(df) = dft$V1
  colnames(df) = paste(cancerType, tumorStage, toString(1), sep = "_")


  for(i in 2:length(manifest_subset$filename)){
    fn = paste(directory, toString(manifest_subset$filename[i]), sep = "/")
    dft = read.table(gzfile(fn))
    dft = head(dft, -5)
    tempSampleName = paste(cancerType, tumorStage, toString(i), sep = "_")
    df[tempSampleName] = dft$V2
  }
  
  return(df)
}

ReadNormal = function(){
  
  normal1 = ReadmiRNA("LUAD", "normal")
  normal2 = ReadmiRNA("LUSC", "normal")
  n = merge(normal1, normal2, by =  0, all = TRUE)
  
  rownames(n) = n$Row.names
  n$Row.names = NULL
  return(n)
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

diffgeneList = function(ct, ts){

  q = normal
  q = add_Label(q, "Normal")
  
  df = ReadmiRNA(ct, ts)
  p = add_Label(df, ct)

  pq = rbind(p, q)
  X = pq[,names(pq) != 'label']
  mRNA.data = list(x = t(X), y = pq$label, genenames = colnames(X), geneid = colnames(X))
  mRNA.train = pamr.train(mRNA.data)
  mRNA.cv = pamr.cv(mRNA.train, mRNA.data)
  deltaIndex = which.min(mRNA.cv$error)
  delta = mRNA.cv$threshold[deltaIndex]
  mRNA.confusion = pamr.confusion(mRNA.cv, delta)

  lg = pamr.listgenes(mRNA.train, mRNA.data, threshold=delta)
  lg = data.frame(lg)

  dl = df[row.names(df) %in% lg$id,]
  return(list(dl, df))
}
compareStage = function(ct){
  ia = data.frame(diffgeneList(ct, tumorStageList[1])[1])
  ib = data.frame(diffgeneList(ct, tumorStageList[2])[1])
  iia = data.frame( diffgeneList(ct, tumorStageList[3])[1])
  iib = data.frame(diffgeneList(ct, tumorStageList[4])[1])
  iii = data.frame(diffgeneList(ct, tumorStageList[5])[1])
  iv = data.frame(diffgeneList(ct, tumorStageList[6])[1])
  
  ia2 = add_Label(ia, tumorStageList[1])
  ib2 = add_Label(ib, tumorStageList[2])
  iia2 = add_Label(iia, tumorStageList[3])
  iib2 = add_Label(iib, tumorStageList[4])
  iii2 = add_Label(iii, tumorStageList[5])
  iv2 = add_Label(iv, tumorStageList[6])
  
  column = list(ia = colnames(ia2), ib = colnames(ib2), 
                iia = colnames(iia2), iib = colnames(iib2), 
                iii = colnames(iii2), iv = colnames(iv2))

  common.all = Reduce(intersect, column)
  commom.stage_i = intersect(column.ia, column.ib)

}
acgts = function(ts){
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
  mRNA.data = list(x = t(X), y = adsq$label, genenames = colnames(X), geneid = colnames(X))
  
  mRNA.train = pamr.train(mRNA.data)
  mRNA.cv = pamr.cv(mRNA.train, mRNA.data)
  deltaIndex = which.min(mRNA.cv$error)
  delta = mRNA.cv$threshold[deltaIndex]
  mRNA.confusion = pamr.confusion(mRNA.cv, delta)
  
  lg = pamr.listgenes(mRNA.train, mRNA.data, threshold=delta)
  lg = data.frame(lg)
  write.csv(lg$id, paste("/Users/zhongningchen/Data/mRNA/shrukenCentroid", ts, sep = "/"), quote = FALSE, row.names = FALSE)
  return(lg$id)
}

library(pamr)
parent = "/Users/zhongningchen/Data/mRNA"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")
normal = ReadNormal()
ts = tumorStageList[5]
ct = "LUAD"









