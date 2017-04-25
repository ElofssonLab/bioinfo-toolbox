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


                                        #TEst pcons
noconstr.tmp<-data[which(data$mindist == 0) ,]
constr.tmp<-data[which(data$mindist > 0) ,]

constr.list <-c(proq3d.constr$target)
noconstr.list <-c(proq3d.noconstr$target)

constr<-constr.tmp[constr.tmp$target %in% noconstr.list,]
noconstr<-noconstr.tmp[noconstr.tmp$target %in% constr.list,]




temp <- noconstr %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.noconstr <- temp %>%  filter (! duplicated(target))
                                        # commons subset


temp <- constr %>% group_by(target) %>% slice(which(Pcons == max(Pcons)))
pcons.constr <- temp %>%  filter (! duplicated(target))
summary(pcons.noconstr$TM)
summary(pcons.constr$TM)

temp <- noconstr %>% group_by(target) %>% slice(which(ProQ3D == max(ProQ3D)))
proq3d.noconstr <- temp %>%  filter (! duplicated(target))
temp <- constr %>% group_by(target) %>% slice(which(ProQ3D == max(ProQ3D)))
proq3d.constr <- temp %>%  filter (! duplicated(target))
summary(proq3d.noconstr$TM)
summary(proq3d.constr$TM)


temp <- noconstr %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.noconstr <- temp %>%  filter (! duplicated(target))
temp <- constr %>% group_by(target) %>% slice(which(TM == max(TM)))
tm.constr <- temp %>%  filter (! duplicated(target))
summary(tm.noconstr$TM)
summary(tm.constr$TM)


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

