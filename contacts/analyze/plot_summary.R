library(dplyr)
library(ggplot2)


                                        #data<-read.table("summary-all.csv",header=TRUE,sep=",",fill=TRUE,stringsAsFactors = FALSE)

data<-read.csv("summary-all.csv",header=TRUE,stringsAsFactors=F,strip.white=T)

                                        #colnames(data) <- c("target","ali","num","mindist","maxdist","length","model","TM","Pcons","cns","noe","ProQ2D","ProQ3D")


                                        # Best selection methods..

temp <- data %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons <- temp %>%  filter (! duplicated(target))

temp <- data %>% group_by(target) %>% slice(which(noe == max(noe)))
noe <- temp %>%  filter (! duplicated(target))

temp <- data %>% group_by(target) %>% slice(which(ProQ3D == max(ProQ3D)))
proq3d <- temp  %>% filter (! duplicated(target))

temp <- data %>% group_by(target) %>% slice(which(TM == max(TM)))
tm <- temp %>%  filter (! duplicated(target))
                                        #data.rank <- data %>% group_by(target) %>% slice(which(model == "fa_1"))

summary(as.numeric(pcons$TM))
summary(as.numeric(noe$TM))
summary(as.numeric(proq3d$TM))
summary(as.numeric(tm$TM))
                                        #summary(as.numeric(data.rank$TM))


                                        # Check max och mindist
mindist20<-data[which(data$mindist == 20) ,]
mindist30<-data[which(data$mindist == 30) ,]
mindist40<-data[which(data$mindist == 40) ,]
mindist50<-data[which(data$mindist == 50) ,]
maxdist20<-data[which(data$maxdist == 20) ,]
maxdist30<-data[which(data$maxdist == 30) ,]
maxdist40<-data[which(data$maxdist == 40) ,]
maxdist50<-data[which(data$maxdist == 50) ,]

temp <- mindist20 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.min20 <- temp %>%  filter (! duplicated(target))
temp <- mindist30 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.min30 <- temp %>%  filter (! duplicated(target))
temp <- mindist40 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.min40 <- temp %>%  filter (! duplicated(target))
temp <- mindist50 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.min50 <- temp %>%  filter (! duplicated(target))
temp <- maxdist20 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.max20 <- temp %>%  filter (! duplicated(target))
temp <- maxdist30 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.max30 <- temp %>%  filter (! duplicated(target))
temp <- maxdist40 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.max40 <- temp %>%  filter (! duplicated(target))
temp <- maxdist50 %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.max50 <- temp %>%  filter (! duplicated(target))

temp <- mindist20 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.min20 <- temp %>%  filter (! duplicated(target))
temp <- mindist30 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.min30 <- temp %>%  filter (! duplicated(target))
temp <- mindist40 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.min40 <- temp %>%  filter (! duplicated(target))
temp <- mindist50 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.min50 <- temp %>%  filter (! duplicated(target))
temp <- maxdist20 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.max20 <- temp %>%  filter (! duplicated(target))
temp <- maxdist30 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.max30 <- temp %>%  filter (! duplicated(target))
temp <- maxdist40 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.max40 <- temp %>%  filter (! duplicated(target))
temp <- maxdist50 %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.max50 <- temp %>%  filter (! duplicated(target))

summary(as.numeric(pcons.min20$TM))
summary(as.numeric(pcons.min30$TM))
summary(as.numeric(pcons.min40$TM))
summary(as.numeric(pcons.min50$TM))
summary(as.numeric(pcons.max20$TM))
summary(as.numeric(pcons.max30$TM))
summary(as.numeric(pcons.max40$TM))
summary(as.numeric(pcons.max50$TM))

summary(as.numeric(tm.min20$TM))
summary(as.numeric(tm.min30$TM))
summary(as.numeric(tm.min40$TM))
summary(as.numeric(tm.min50$TM))
summary(as.numeric(tm.max20$TM))
summary(as.numeric(tm.max30$TM))
summary(as.numeric(tm.max40$TM))
summary(as.numeric(tm.max50$TM))







                                        # We need to find common subset
noconstr.tmp<-data[which(data$mindist == 0) ,]
constr.tmp<-data[which(data$mindist > 0) ,]
constr.list <-constr.tmp %>%  filter (! duplicated(target)) 
noconstr.list <-noconstr.tmp %>%  filter (! duplicated(target)) 
constr<-constr.tmp[constr.tmp$target %in% noconstr.list$target,]
noconstr<-noconstr.tmp[noconstr.tmp$target %in% constr.list$target,]



                                        # commons subset


temp <- noconstr %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.noconstr <- temp %>%  filter (! duplicated(target))
temp <- constr %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.constr <- temp %>%  filter (! duplicated(target))
temp <- noconstr %>% group_by(target) %>% slice(which(ProQ3D == max(ProQ3D)))
proq3d.noconstr <- temp %>%  filter (! duplicated(target))
temp <- constr %>% group_by(target) %>% slice(which(ProQ3D == max(ProQ3D)))
proq3d.constr <- temp %>%  filter (! duplicated(target))
temp <- noconstr %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.noconstr <- temp %>%  filter (! duplicated(target))
temp <- constr %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.constr <- temp %>%  filter (! duplicated(target))

summary(pcons.noconstr$TM)
summary(pcons.constr$TM)
summary(proq3d.noconstr$TM)
summary(proq3d.constr$TM)
summary(tm.noconstr$TM)
summary(tm.constr$TM)







                                        # Not really neededto merge all to one dataset.
methods <- list(pcons,noe,proq3d,tm)
names(methods) <- c("pcons.", "noe.", "proq3d.", "tm.") #just manually naming, presumably your list has names
mycolnames <- colnames(pcons)
mycolnames1 <- lapply(names(methods), function(x) paste0(x, mycolnames))
for(i in 1:length(methods)){
    colnames(methods[[i]]) <- mycolnames1[[i]]
    colnames(methods[[i]])[1] <- "target" #put Var back in so you can merge
}
## Change column names by pasting name of dataframe in list with standard column names. - using ugly mix of `lapply` and a `for` loop:
merge.all <- function(x, y)
    merge(x, y, all=TRUE, by="target")

m <- Reduce(merge.all, methods)
                                        # Now we have the data in a single table.
plot(m$pcons.TM,m$proq3d.TM)

