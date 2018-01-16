library(pamr)
predictors = read.csv('/Users/zhongningchen/LungCancer/Data')
label = read.csv('/Users/zhongningchen/LungCancer/Label')
predictors = predictors[,2:ncol(predictors)]
label = label[,2:ncol(label)]
miRNA.data = list(x = t(predictors), y = label, genenames = colnames(predictors))
