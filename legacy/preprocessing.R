library(rcompanion)
library(MASS)
library(e1071)
library(pamr)
parent = "/Users/zhongningchen/Data/miRNA"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")

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

add_Label = function(DF, n){
  #transpose dataframe to n x p
  #adds label column 
  DF = as.data.frame(t(DF))
  DF$label = rep(n, nrow(DF))
  return(DF)
}

norm = function(Turbidity){
  Box = boxcox(Turbidity ~ 1, lambda = seq(-6,6,0.1))
  Cox = data.frame(Box$x, Box$y)
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]  
  T_box = (Turbidity ^ lambda - 1)/lambda
  return(T_box)
}

norm2 = function(adsq, miR){
  miR = miR + 1
  Box = boxcox(miR ~ adsq$label,
               data = adsq,
               lambda = seq(-6,6,0.1)
  )
  
  Cox = data.frame(Box$x, Box$y)
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  Cox2[1,]
  lambda = Cox2[1, "Box.x"]
  if (lambda != 0){
    return((miR ^ lambda - 1)/lambda)   
  }else{
    return(log(miR))
  } 
  
}

near_zero_var = function(X, i){
  if(all(X[,i] == 0)){
    return(FALSE)
  }else{
    n <- length(X[,i])
    u = as.double(length(unique(X[,i]))/n)
    s = sort(X[,i],partial=n-1)
    diff = as.double(s[n]/s[n-1])
    
    return(!((u < .1) & (diff>20)))
  }
}

skew = function(X, i){
  return((abs(kurtosis(X[,i])) >= 2))
}
{
df1 = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[1])
df2 = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[2])
df3 = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[3])
df4 = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[4])
df5 = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[5])
df6 = ReadmiRNA("LUAD", "PrimaryTumor", tumorStageList[6])

adf1 = add_Label(df1, tumorStageList[1])
adf2 = add_Label(df2, tumorStageList[2])
adf3 = add_Label(df3, tumorStageList[3])
adf4 = add_Label(df4, tumorStageList[4])
adf5 = add_Label(df5, tumorStageList[5])
adf6 = add_Label(df6, tumorStageList[6])
}
adsq_ = rbind(adf1, adf2, adf3, adf4, adf5, adf6)
gene = colnames(adsq_)

#eliminating al near zero variance attributes
NZV_gene = data.frame(matrix(TRUE, 1882, 1))
for(i in 1:1881){
  NZV_gene[i,] = near_zero_var(X = adsq_[,names(adsq_) != 'label'], i)
}
NZV = colnames(adsq_)[NZV_gene == TRUE]
adsq_nzv = adsq_[,colnames(adsq_) %in% NZV]

{
###Finding all the skewed data and apply boxcox
# SkewList = (data.frame((matrix(TRUE, 1512, 1))))
# for(i in 1:1511){
#   SkewList[i,] = skew(X = adsq_nzv[,names(adsq_nzv) != 'label'], i)
# }
# adsq_skew = adsq_nzv[,colnames(adsq_nzv) %in% colnames(adsq_nzv)[SkewList == TRUE]]
# 
# 
# c1 = norm2(adsq = adsq_skew, miR =  adsq_skew[,1])
# df = data.frame(c1)
# colnames(df) = colnames(adsq_nzv)[SkewList == TRUE][1]
# 
# for (i in 2:1504){
#   df[colnames(adsq_nzv)[SkewList == TRUE][i]] = norm2(adsq = adsq_skew, miR = adsq_skew[,i])
#}
}
c1 = norm2(adsq = adsq_nzv, miR =  adsq_nzv[,1])
df = data.frame(c1)
colnames(df) = colnames(adsq_nzv)[1]

for (i in 2:1511){
  df[colnames(adsq_nzv)[i]] = norm2(adsq = adsq_nzv, miR = adsq_nzv[,i])
}
df$label = adsq_nzv$label

m = mean(df[,1])
s = sd(df[,1])
c1 = (adsq[,1] - 1)/s
adsq = data.frame(c1)
colnames(adsq) = colnames(df)[1]

for(i in 1:1511){
  m = mean(df[,i])
  s = sd(df[,i])
  adsq[colnames(df)[i]] = (df[,i] - m)/s
}
adsq$label = df$label

# stage_ia = adsq_nzv[adsq_nzv$label==tumorStageList[1],]
# stage_ib = adsq_nzv[adsq_nzv$label==tumorStageList[2],]
# stage_iia = adsq_nzv[adsq_nzv$label==tumorStageList[3],]
# 
# adsq = rbind(stage_ia, stage_iia)

adsq_prime = adsq
adsq = df
adsq = adsq_prime

adsq = rbind(adsq_prime[adsq_prime$label == tumorStageList[1],], 
             adsq_prime[adsq_prime$label == tumorStageList[4],])
  


X = adsq[,names(adsq) != 'label']
miRNA.data = list(x = t(X), y = adsq$label, genenames = colnames(X), geneid = colnames(X))

miRNA.train = pamr.train(miRNA.data)
miRNA.cv = pamr.cv(miRNA.train, miRNA.data)
deltaIndex = which.min(miRNA.cv$error)
delta = miRNA.cv$threshold[deltaIndex]
miRNA.confusion = pamr.confusion(miRNA.cv, delta)



bart = data.frame(colnames(X))
for (i in 1:(length(adsq) - 1)){
  bart[i,2] = bartlett.test(adsq[,i] ~ adsq$label, adsq)[3]
}

dev = bart$colnames.X.[bart$p.value < 0.05]
boxplot(adsq[,dev[3]] ~ adsq$label, data = adsq)

for (i in 1:length(adsq)){
  boxplot(adsq_nzv[,i] ~ adsq_nzv$label, data = adsq_nzv)
  x = readline(prompt="Press [enter] to continue")
  print (i)
  boxplot(adsq[,i] ~ adsq$label, data = adsq)
  x = readline(prompt="Press [enter] to continue")
  print (i)
  if (x == 1){
    i = 1999
  }
}


boxplot(c1 ~ adsq$label,
        data = adsq)




# for(i in 1:1881){
#   temp = kruskal.test(adsq[,i] ~ as.factor(adsq$label), data = adsq)[[3]]
#   zebra[[i]] = temp
# }
# dog = zebra < 0.005
# diff = gene[dog &!is.na(dog)]
# 
# m = mean(adsq[])
# i = 12
# boxplot(adsq[diff[i]][,1] ~ as.factor(adsq$label), data = adsq)
