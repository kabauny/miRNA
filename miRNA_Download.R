library(TCGAbiolinks)
library(DT)
library(DESeq)
########################
query_isoform_LUSC <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Isoform Expression Quantification",
                  experimental.strategy = "miRNA-Seq"
                  )
GDCdownload(query_isoform_LUSC, method = "api", files.per.chunk = 10)
isoform_data_LUSC <- GDCprepare(query_isoform_LUSC)
write.csv(isoform_data_LUSC, 'isoform_data_LUSC.csv')
########################

query_isoform_LUAD <- GDCquery(project = "TCGA-LUAD",
                               data.category = "Transcriptome Profiling",
                               data.type = "Isoform Expression Quantification",
                               experimental.strategy = "miRNA-Seq"
)
GDCdownload(query_isoform_LUAD, method = "api", files.per.chunk = 10)
isoform_data_LUAD <- GDCprepare(query_isoform_LUAD)
write.csv(isoform_data_LUAD, 'isoform_data_LUAD.csv')
########################

query_miRNA_LUSC <- GDCquery(project = "TCGA-LUSC",
                          data.category = "Transcriptome Profiling",
                          data.type = "miRNA Expression Quantification",
                          experimental.strategy = "miRNA-Seq"
)
GDCdownload(query_miRNA_LUSC, method = "api", files.per.chunk = 10)
miRNA_data_LUSC <- GDCprepare(query_miRNA_LUSC)

LUSC_name = miRNA_data_LUSC$miRNA_ID
miRNA_data_LUSC$miRNA_ID <- NULL
LUSC = miRNA_data_LUSC[c(TRUE, FALSE, FALSE)]
LUSC$name = LUSC_name
LUSC$label = rep('LUSC', nrow(LUSC))

write.csv(LUSC, 'miRNA_LUSC.csv')

########################
query_miRNA_LUAD <- GDCquery(project = "TCGA-LUAD",
                             data.category = "Transcriptome Profiling",
                             data.type = "miRNA Expression Quantification",
                             experimental.strategy = "miRNA-Seq"
)
GDCdownload(query_miRNA_LUAD, method = "api", files.per.chunk = 10)
miRNA_data_LUAD <- GDCprepare(query_miRNA_LUAD)

LUAD_name = miRNA_data_LUAD$miRNA_ID
miRNA_data_LUAD$miRNA_ID <- NULL
LUAD = miRNA_data_LUAD[c(TRUE, FALSE, FALSE)]
LUAD$name = LUAD_name
LUAD$label = rep('LUAD', nrow(LUAD))

write.csv(LUAD, 'miRNA_LUAD.csv')
