library(pamr)
parent = "/Users/zhongningchen/Data"

manifest = function(c, n){

  f = paste(parent, c, paste(n, "txt", sep = "."), sep = "/")
  fi = data.frame(read.table(f, sep = "\t", header = TRUE))
  fi$md5 = NULL
  fi$size = NULL
  fi$state = NULL
  fi$loc = paste(fi$id, fi$filename, sep = "/")
  
  return(data.frame(fi$loc))
}

ReadmiRNA = function(ct, loc){
  
  f = paste(parent, ct, "Data_Full", sep = "/")
  
  fileName = paste(f, toString(loc$id[1]), sep = "/")
  dft = read.table(fileName, sep = "\t", header = TRUE)
  df = data.frame(dft$read_count)
  colnames(df)[1] = paste(ct, mutation, "1", sep = "_")
  
  rownames(df) = dft$miRNA_ID
  
  for(i in 2:length(loc$id)){
    fileName = paste(f, toString(loc$id[i]), sep = "/")
    dft = read.table(fileName, sep = "\t", header = TRUE)
    tempSampleName = paste(ct, mutation, toString(i), sep = "_")
  
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

###################################################

mg = function(ct, mutation)

  mutationLoc = manifest(ct, mutation)
  colnames(mutationLoc) = c("id")
  
  mut = ReadmiRNA(ct, mutationLoc)
  mut = add_Label(mut, mutation)
  
  tumLoc_ = manifest(ct, "Tumor")
  colnames(tumLoc_) = c("id")
  tumLoc = data.frame(tumLoc_[!(tumLoc_$id %in% mutationLoc$id),])
  colnames(tumLoc) = c("id")
  
  tum = ReadmiRNA(ct, tumLoc)
  tum = add_Label(tum, "comp")
  
  adsq = join(a = tum, b = mut)
  X = adsq[,names(adsq) != 'label']
  miRNA.data = list(x = t(X), y = adsq$label, genenames = colnames(X), geneid = colnames(X))

  return (miRNA.data)

miRNA = mg("LUAD", "P53")
pamr.menu(miRNA)
