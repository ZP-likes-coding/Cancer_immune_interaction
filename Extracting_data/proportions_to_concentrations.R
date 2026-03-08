library(SeuratObject)


tismo_data <-readRDS('../../Data/TISMO_vivo_compiled.rds')
study_ids <- tismo_data$Study_ID

filtered_tismo_data <-  subset(tismo_data,subset = Study_ID =='GSE124821')
filtered_tismo_data[['']]

metadata <-  names(filtered_tismo_data@meta.data)
quanTIseq_labels <- metadata[grepl('quanTIseq',metadata)]
Idents(filtered_tismo_data) <- quanTIseq_labels
quanTIseq_data <- subset(filtered_tismo_data,idents= 'B_quanTIseq' )
layer <- 
