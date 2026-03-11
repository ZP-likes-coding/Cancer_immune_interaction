install.packages("SeuratObject")

obj <- readRDS("/Users/manoharan/Downloads/TISMO_vivo_compiled.rds")
str(obj)
class(obj)
slotNames(obj)
obj@assays
obj[["RNA"]]
rownames(obj)[1:10]
colnames(obj)[1:10]

head(obj@meta.data)
colnames(obj@meta.data)

