# Description:
# Do chisq test given GO frequency
# 

# Created 2013-09-12, updated 2013-09-12, Nanjiang Shu 

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

args <- commandArgs(trailingOnly = TRUE)

numFile <- length(args)
cat ("numFile =", numFile, "\n")
filelist <- c()
for (i in 1:numFile){
    filelist <- c(filelist, args[i])
}
taglist <- gsub(".gofreq.txt", "", filelist)

cat ("file list=", filelist, "\n")
cat ("tag list=", taglist, "\n")


tablelist <- list()

for (i in 1:numFile){
    tablelist[[i]] <- read.table(filelist[i], header=F, sep="\t")
}

sumvector <- c()
for (i in 1:numFile){
    sumvector <- c(sumvector, sum(tablelist[[i]]$V3))
}


# 1. Analyzing the distribution of "molecular function" for proteins at various
# topology variation categories

pvaluelist <- c()
for (i in 1:numFile){
    for (j in 1:numFile) {
        x <- tablelist[[i]]$V3
        y <- tablelist[[j]]$V3
        m <- rbind(v1=x, v2=y)
        m <- t(m)
        mfiltered <- m[ rowSums(m)>=30, ]
        rst <- chisq.test(mfiltered)
        pvaluelist <- c(pvaluelist, rst$p.value)
        cat ("#====================================\n")
        cat (taglist[i] , " -- " , taglist[j], "\n")
        print (rst$observed)
        print (rst)
        warnings()
    }
}

# show p-value of comparisons between different distributions
print("P-value matrix between different topology variation categories.")
print("Alternative hypothesis, the distribution of \"molecular function\" is different.")
mtx <- matrix(pvaluelist, nrow=numFile, ncol=numFile, dimnames=list(taglist, taglist))
print(mtx)
print(sprintf("%.3g",mtx))

# 2. Analyzing whether a certain topology variation category is overrepresented
# in a GeneOntology "molecular function", every category is compared to the
# background
# in this case, the last input file should be the background

hypolist <- c("two.sided", "greater", "less")
SEL_GOID_LIST <- c("GO:0060089", "GO:0009055", "GO:0005488", "GO:0005215",
"GO:0005198", "GO:0003824")

sel_tablelist <- list()
for (i in 1:numFile){
    tb <- tablelist[[i]]
    sel_tablelist[[i]] <- tb[tb$V1 %in% SEL_GOID_LIST, ]
}

bg_table <- sel_tablelist[[numFile]] # background table
go_termvector <- as.vector(bg_table$V2)
go_termvector <- trim(go_termvector)
num_row <- nrow(bg_table)
num_hypo <- length(hypolist)
pvaluelist <- list()
for (igo in 1:num_row){
    cat ("GOTERM:",go_termvector[igo],"\n")
    pvaluevector <- c()
    for (i in 1:(numFile-1)){
        cat ("TAG:",taglist[i],"\n")
        tb <- sel_tablelist[[i]]
        v1 <- c(tb$V3[igo], bg_table$V3[igo])
        v2 <- c(sumvector[i], sumvector[numFile])
        for (ihypo in 1:num_hypo){
            rst <- prop.test(v1, v2, alternative = hypolist[ihypo])
            pvaluevector <- c(pvaluevector, rst$p.value)
            cat ("v1=",v1, "v2=",v2, "hypo=", hypolist[ihypo],"\n")
            print (rst$observed)
            print (rst)
            warnings()
        }
    }
    #mtx <- matrix(pvaluevector, nrow=numFile-1, ncol=num_hypo, dimnames=list(taglist[1:numFile-1], hypolist))
    mtx <- matrix(pvaluevector, nrow=num_hypo, ncol=numFile-1, dimnames=list(hypolist, taglist[1:numFile-1]))
    pvaluelist[[igo]] <- mtx

}
for (igo in 1:num_row){
    mtx <- pvaluelist[[igo]]
    cat("\n")
    cat (go_termvector[igo], "\n")
    print ("%e", mtx)

}

