f_micro = "/Users/zhongningchen/Data/miRNA/geneList"
f_RNA = "/Users/zhongningchen/Data/mRNA/geneList"

tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")

micro.stage_ia = read.csv(paste(f_micro, tumorStageList[1], sep = "/"))
micro.stage_ib = read.csv(paste(f_micro, tumorStageList[2], sep = "/"))
micro.stage_iia = read.csv(paste(f_micro, tumorStageList[3], sep = "/"))
micro.stage_iib = read.csv(paste(f_micro, tumorStageList[4], sep = "/"))
micro.stage_iii = read.csv(paste(f_micro, tumorStageList[5], sep = "/"))
micro.stage_iv = read.csv(paste(f_micro, tumorStageList[6], sep = "/"))

RNA.stage_ia = read.csv(paste(f_RNA, tumorStageList[1], sep = "/"))
RNA.stage_ib = read.csv(paste(f_RNA, tumorStageList[2], sep = "/"))
RNA.stage_iia = read.csv(paste(f_RNA, tumorStageList[3], sep = "/"))
RNA.stage_iib = read.csv(paste(f_RNA, tumorStageList[4], sep = "/"))
RNA.stage_iii = read.csv(paste(f_RNA, tumorStageList[5], sep = "/"))
RNA.stage_iv = read.csv(paste(f_RNA, tumorStageList[6], sep = "/"))

# temp1 = data.frame(micro.stage_ia$ensembl_gene_id[micro.stage_ia$ensembl_gene_id %in% RNA.stage_ia$ensembl_gene_id])
# temp2 = data.frame(micro.stage_ib$ensembl_gene_id[micro.stage_ib$ensembl_gene_id %in% RNA.stage_ib$ensembl_gene_id])
# temp3 = data.frame(micro.stage_iia$ensembl_gene_id[micro.stage_iia$ensembl_gene_id %in% RNA.stage_iia$ensembl_gene_id])
# temp4 = data.frame(micro.stage_iib$ensembl_gene_id[micro.stage_iib$ensembl_gene_id %in% RNA.stage_iib$ensembl_gene_id])
# temp5 = data.frame(micro.stage_iii$ensembl_gene_id[micro.stage_iii$ensembl_gene_id %in% RNA.stage_iii$ensembl_gene_id])
# temp6 = data.frame(micro.stage_iv$ensembl_gene_id[micro.stage_iv$ensembl_gene_id %in% RNA.stage_iv$ensembl_gene_id])

temp1 = micro.stage_ia$mirbase_id[micro.stage_ia$mirbase_id %in% micro.stage_ib$mirbase_id]
temp2 = micro.stage_ia$mirbase_id[micro.stage_ia$mirbase_id %in% micro.stage_iia$mirbase_id]
temp3 = micro.stage_ia$mirbase_id[micro.stage_ia$mirbase_id %in% micro.stage_iib$mirbase_id]
temp4 = micro.stage_ia$mirbase_id[micro.stage_ia$mirbase_id %in% micro.stage_iii$mirbase_id]
temp5 = micro.stage_ia$mirbase_id[micro.stage_ia$mirbase_id %in% micro.stage_iv$mirbase_id]


# temp2 = data.frame(micro.stage_ib$ensembl_gene_id[micro.stage_ib$ensembl_gene_id %in% RNA.stage_ib$ensembl_gene_id])
# temp3 = data.frame(micro.stage_iia$ensembl_gene_id[micro.stage_iia$ensembl_gene_id %in% RNA.stage_iia$ensembl_gene_id])
# temp4 = data.frame(micro.stage_iib$ensembl_gene_id[micro.stage_iib$ensembl_gene_id %in% RNA.stage_iib$ensembl_gene_id])
# temp5 = data.frame(micro.stage_iii$ensembl_gene_id[micro.stage_iii$ensembl_gene_id %in% RNA.stage_iii$ensembl_gene_id])
# temp6 = data.frame(micro.stage_iv$ensembl_gene_id[micro.stage_iv$ensembl_gene_id %in% RNA.stage_iv$ensembl_gene_id])
