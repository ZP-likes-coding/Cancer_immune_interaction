library(SeuratObject)


tismo_data <-readRDS('/Users/manoharan/Desktop/MIB MP/TISMO_vivo_compiled.rds')
study_ids <- tismo_data$Study_ID

filtered_tismo_data <-  subset(tismo_data,subset = Study_ID =='GSE124821')

metadata <-  names(filtered_tismo_data@meta.data)
quanTIseq_labels <- metadata[grepl('quanTIseq',metadata)]

quanTIseq_labels <- c("Dendritic_quanTIseq","NK_quanTIseq","T.CD8_quanTIseq","Tregs_quanTIseq")
quanTIseq_labels <- append(quanTIseq_labels,c('Mouse_treatment','Condition','Timepoint','GSM_ID', 'Cell_Line'))
#Idents(filtered_tismo_data) <- quanTIseq_labels
#quanTIseq_data <- subset(filtered_tismo_data,idents= 'B_quanTIseq' )

quanTIseq_data <-  filtered_tismo_data[[quanTIseq_labels]]


masses = data.frame(lymphocyte = 2e-10,dendritic = 6e-10,NK=6e-10)
show(masses)

cell_density <- 1.0e8
show(masses)
show(quanTIseq_data)
# interfering 
rowSums(quanTIseq_data[1:10, c("Dendritic_quanTIseq", "NK_quanTIseq", "T.CD8_quanTIseq", "Tregs_quanTIseq")])





masses[['dendritic']]

quanTIseq_data$Dendritic_quanTIseq <- quanTIseq_data$Dendritic_quanTIseq*(masses[['dendritic']][1])
quanTIseq_data$Dendritic_quanTIseq <- quanTIseq_data$Dendritic_quanTIseq*cell_density

quanTIseq_data$NK_quanTIseq <- quanTIseq_data$NK_quanTIseq*(masses[['NK']][1])
quanTIseq_data$NK_quanTIseq <- quanTIseq_data$NK_quanTIseq*cell_density

quanTIseq_data$T.CD8_quanTIseq <- quanTIseq_data$T.CD8_quanTIseq*(masses[['lymphocyte']][1])
quanTIseq_data$T.CD8_quanTIseq <- quanTIseq_data$T.CD8_quanTIseq*cell_density

quanTIseq_data$Tregs_quanTIseq <- quanTIseq_data$Tregs_quanTIseq*(masses[['lymphocyte']][1])
quanTIseq_data$Tregs_quanTIseq <- quanTIseq_data$Tregs_quanTIseq*cell_density

show(quanTIseq_data)

# saving the quanTIseq_data 
write.csv(quanTIseq_data, '~/Desktop/GSE124821_QSdata.csv')

        
