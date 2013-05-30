rm(list = ls())

#data <- read.table("pfam_stat.csv",header=TRUE)
data <- read.table("pdb_stat_2013-02.csv",header=TRUE)
#data <- read.table("tmp.csv",header=TRUE, sep="\t")
data2 <- read.table("../hmmscan_Pfam-SwissProt/pfam_stat_2013-02_select_max3.csv", header=TRUE)
#selection <- data[which(data[,4] >= 3 & data[,7] > 500),]
#selection <- data[which(data[,7] >= 7 & data[,8] > 500 & data[,2] > 25 & data[,9] == 'None'),]
selection <- data[which(data[,1] %in% data2[,1]),]
selection <- selection[which(selection[,7] > 1),]
selection.sorted <- selection[order(selection[,14],selection[,15]),]
#write.table(selection.sorted, file = "pdb_stat_2013-02_select.csv", quote=FALSE, sep="\t",row.names=FALSE)
write.table(selection.sorted, file = "pdb_stat_2013-02_select2.csv", quote=FALSE, sep="\t",row.names=FALSE)

selection.unique <- selection.sorted[!duplicated(selection.sorted[,11]),]
write.table(selection.unique[,c(11,13,15)], file="query.txt", quote=FALSE, sep="\t",row.names=FALSE)
