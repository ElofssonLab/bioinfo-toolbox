library(dplyr)

data <-read.table("foo.txt",fill=TRUE)
colnames(data) <- c("target","ali","num","mindist","maxdist","length","model","TM","Pcons","cns","noe","ProQ2D","ProQ3D")

data.pcons <- data %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
data.noe <- data %>% group_by(target) %>% slice(which(noe == max(noe)))
data.proq3d <- data %>% group_by(target) %>% slice(which(ProQ3D == max(ProQ3D)))

data.rank <- data %>% group_by(target) %>% slice(which(model == "fa_1"))


mean(data.pcons$TM)
mean(data.noe$TM)
mean(data.proq3d$TM)
mean(data.rank$TM)
