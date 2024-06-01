setwd("C://Users//penguin//Desktop//学习资料//高通量测序02//大作业//scRNA")
mt <- Matrix::readMM("GSM5975238_scR104_matrix.mtx")
dim(mt)
#[1] 31053  2846
#行为基因，列为细胞

format(object.size(mt), units = "Mb")
#[1] "67.9 Mb"

bc <- read.table("GSM5975238_scR104_barcodes.tsv")
genes <- read.table("GSM5975238_scR104_features.tsv", sep = "\t")
dim(bc)
#[1] 2846      1
dim(genes)
#[1] 31053     3

row.names(mt) <- genes$V1
colnames(mt) <- bc$V1

saveRDS(mt,"scRNA_matrix.Rds")