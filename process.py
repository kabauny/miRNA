#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 08:59:52 2018

@author: zhongningchen
"""
import pandas as pd


parent = '/Users/zhongningchen/LungCancer/miRNA_Data/'
#file_list = ['miRNA_LUAD.csv', 'miRNA_LUSC.csv', 'miRNA_BLCA.csv', 
#             'miRNA_BRCA.csv', 'miRNA_COAD.csv', 'miRNA_PRAD.csv']
ulabel = ['LUAD', 'LUSC']
file_list = ['miRNA_LUAD.csv', 'miRNA_LUSC.csv']

def load():
    dfList = []
    for file in file_list:
        df_temp = pd.read_csv(parent + file)
        del df_temp['Unnamed: 0']
        
        temp_dict = {}
        for i in range(len(df_temp['name'])):
            temp_dict[i] = df_temp['name'][i]
        df_temp = df_temp.rename(index = temp_dict)
        del df_temp['name']
        dfList += [df_temp.transpose()]
        
    return dfList

def combined():
    #loads all data into dataframe dfList
    #for miRNA Lung cancer, it loads LUAD and LUSC as dataframes each
    #and stores to a list called dfList

    dfList = []
    label = []
    for file in file_list:
        df_temp = pd.read_csv(parent + file)
        del df_temp['Unnamed: 0']
        
        temp_dict = {}
        for i in range(len(df_temp['name'])):
            temp_dict[i] = df_temp['name'][i]
        df_temp = df_temp.rename(index = temp_dict)
        
        del df_temp['name']
        

        dfList += [df_temp.transpose()]
        label += [file]*len(df_temp.transpose())
        
    DF = pd.concat(dfList)
    DF = DF.fillna(0)
    

    return DF, pd.DataFrame(label)

parent = '/Users/zhongningchen/LungCancer/miRNA_Data/'
#file_list = ['miRNA_LUAD.csv', 'miRNA_LUSC.csv', 'miRNA_BLCA.csv', 
#             'miRNA_BRCA.csv', 'miRNA_COAD.csv', 'miRNA_PRAD.csv']
ulabel = ['LUAD', 'LUSC']
file_list = ['miRNA_LUAD.csv', 'miRNA_LUSC.csv']

df, y = combined()
df.to_csv('/Users/zhongningchen/LungCancer/Data')
y.to_csv('/Users/zhongningchen/LungCancer/Label')












































