# Filename: mergeGOfreqdata.R
# Description:
# merge gofreq.txt data
# calculate the frequency of GOTERMs for a certain topology category,
# normalized by the total counts of that GOTERM
# Calculate also the confidence interval

# Created 2013-09-16, updated 2013-10-08, Nanjiang Shu 

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


# parsing arguments
args <- commandArgs(trailingOnly = TRUE)
numArg <- length(args)
cat ("numArg =", numArg, "\n")
filevector <- c()
basenamevector <- c()
tagvector <- c()
stemnamevector <- c()
allFile <- ""
for (i in 1:numArg){
    bsname <- basename(args[i])
    stemname <-  gsub(".gofreq.txt", "", bsname)
    parts <- strsplit(stemname, "\\.")[[1]]
    tag <- parts[length(parts)]
    if (length(grep("all", tolower(tag))>0)){
        allFile <- args[i]
    }else{
        filevector <- c(filevector, args[i])
        basenamevector <- c(basenamevector, bsname)
        stemnamevector <- c(stemnamevector, stemname)
        tagvector <- c(tagvector, tag)
    }
}

numFile <- length(filevector)
idxINV <- -1
idxMixed <- -1
idxSP2TM <- -1
idxTM2SP <- -1
idxTM2GAP <- -1
idxTM2SEQ <- -1

for (i in 1:numFile){
    if (tagvector[i] == "INV"){
        idxINV <- i
    }else if (tagvector[i] == "Mixed"){
        idxMixed <- i
    }else if (tagvector[i] == "SP2TM"){
        idxSP2TM <- i
    }else if (tagvector[i] == "TM2SP"){
        idxTM2SP <- i
    }else if (tagvector[i] == "TM2GAP"){
        idxTM2GAP <- i
    }else if (tagvector[i] == "TM2SEQ"){
        idxTM2SEQ <- i
    }
}

tagvector <- c(tagvector, c("DiffAll", "DiffWithoutTM2SPandTM2SEQ"))
general_stemname <- gsub(tagvector[idxINV], "", stemnamevector[idxINV])
for (i in (numFile+1):(numFile+2)){
    stemnamevector <- c(stemnamevector, sprintf("%s%s", general_stemname, tagvector[i]))
}

cat ("File list: ", filevector, "\n")
cat ("Tag list: ", tagvector, "\n")
cat ("stemname list: ", stemnamevector, "\n")

# check whether table for the background is input
if ( allFile == ""){
    stop("gofreq for the background is not input! Exit.")
}
if ( numFile < 1 ){
    stop("No input file. Exit.")
}

outpath <- dirname(filevector[1]) #set output the same as the data file

tablelist <- list()
sumvector <- c()
for (i in 1:numFile){
    tablelist[[i]] <- read.table(filevector[i], header=F, sep="\t")
}
alltable <- read.table(allFile, header=F, sep="\t")  
sumAll <- sum(alltable$V3)


# create two more tables for DiffAll and DiffWithoutTM2SPandTM2SEQ,
# respectively
# Get count for DiffAll
v_sum1 <- c(0)
for (j in c(idxINV, idxMixed, idxTM2GAP, idxTM2SEQ, idxSP2TM, idxTM2SP)){
    if (j != -1){
        v_sum1 <- v_sum1 +  tablelist[[j]]$V3
    }
}
cat ("v_sum1=", v_sum1, "\n")
v_sum2 <- 0
for (j in c(idxINV, idxMixed, idxTM2GAP)){
    if (j != -1){
        v_sum2 <- v_sum2 +  tablelist[[j]]$V3
    }
}
cat ("v_sum2=", v_sum2, "\n")
tablelist[[numFile+1]] <- as.data.frame(t(rbind(as.vector(alltable$V1), as.vector(alltable$V2), as.integer(v_sum1))))
tablelist[[numFile+2]] <- as.data.frame(t(rbind(as.vector(alltable$V1), as.vector(alltable$V2), as.integer(v_sum2))))


# calculate total counts for every topology variation group
for (i in 1:(numFile+2)){
#    cat ("#### i = ", i, "vector=", as.integer(as.vector(tablelist[[i]]$V3)), "\n")
    sum <- sum(as.integer(as.vector(tablelist[[i]]$V3)))
#    cat ("sum=", sum, "\n")
    sumvector <- c(sumvector, sum)
}

#print ("HERE")


totalcountvector <- alltable$V3
names(totalcountvector) <- alltable$V1 # now, totalcountvector works like a dictionary in python
# select GOs with at least 100 counts
TH_count <- sum(alltable$V3)*0.0035
SEL_GOID_VECTOR <- alltable[alltable$V3>=TH_count,]$V1 
#print (SEL_GOID_VECTOR)
gotermvector <- alltable$V2
gotermvector <- trim(gotermvector)
names(gotermvector) <- alltable$V1


#hypolist <- c("two.sided", "greater", "less")
#SEL_GOID_LIST <- c("GO:0060089", "GO:0009055", "GO:0005488", "GO:0005215", "GO:0005198", "GO:0003824")

sel_tablelist <- list()
for (i in 1:(numFile+2)){
    tb <- tablelist[[i]]
    #print (tb)
    sel_tablelist[[i]] <- tb[tb$V1 %in% SEL_GOID_VECTOR, ]
    #print (sel_tablelist[[i]])
}

# =====================
# 1. output table with confidence interval
# 2. output table with p-value compared to the background, using prop.test
conflevel <- 0.95  #level of confidence interval
matrixlist_confint <- list()
matrixlist_pvalue <- list()
for (i in 1:(numFile+2)){ #{{{
    cat ("Topology:", tagvector[i],"\n")
    v_goterm <- c()
    v_goid <- c()
    v_freq <- c()
    v_count <- c()
    v_total <- c()
    v_min <- c()
    v_max <- c()

    v_pvalue <- c()

    avg <- sumvector[i]/sumAll
    tb <- sel_tablelist[[i]]
#    print (tb)
    num_row <- nrow(tb)
    for (j in 1:num_row){
        goid <- as.vector(tb[j, ]$V1)
        count <- as.integer(as.vector(tb[j, ]$V3))
        totalcount <- totalcountvector[goid]
        cat ("goid =", goid, "count =", count, "totalcount =", totalcount, "\n")
        rst_confint <- prop.test(count, totalcount, conf.level = conflevel)

        freq <- rst_confint$estimate
#        print (rst_confint)
        v_goid <- c(v_goid, goid)
        v_goterm <- c(v_goterm, gotermvector[goid])
        v_freq <- c(v_freq, rst_confint$estimate)
        v_count <- c(v_count, count)
        v_min <- c(v_min, rst_confint$conf.int[1])
        v_max <- c(v_max, rst_confint$conf.int[2])
        v_total <- c(v_total, totalcount)

        hypo <- "two.sided"
        if (freq < avg){
            hypo <- "less"
        } else if (freq > avg){
            hypo <- "greater"
        }else{
            hypo <- "two.sided"
        }
        v1 <- c(count, sumvector[i])
        v2 <- c(totalcount, sumAll)
        rst_prop <- prop.test(v1,v2, alternative = hypo)
        v_pvalue <- c(v_pvalue, rst_prop$p.value)
    }
    mtx_confint <-  rbind(GOID=v_goid, GOTERM=v_goterm, Count=v_count, Freq=v_freq, ci_min=v_min, ci_max=v_max)
    mtx_pvalue <- rbind(GOID=v_goid, GOTERM=v_goterm, Num=v_count, NumAll=v_total, Freq=v_freq, NumGroup=rep(sumvector[i], num_row), NumTotal=rep(sumAll, num_row), avgFreq=rep(sumvector[i]/sumAll, num_row), Pvalue=v_pvalue)
#    print(mtx)
    matrixlist_confint[[i]] <- t(mtx_confint)
    matrixlist_pvalue[[i]] <- t(mtx_pvalue)
}
#}}}

# =====================

# ======================
# output the result,
# 1. frequency with confidence interval
#outpath <- "goana"
for (i in 1:(numFile+2)){
    cat ("Topology:", tagvector[i],"\n")
    cat("   Matrix with confidence interval\n")
    mtx <- matrixlist_confint[[i]]
    print (mtx)
    outfile <- sprintf("%s/%s.merged.withconfint.txt", outpath, stemnamevector[i])
    write.table(mtx, outfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

    cat("\n   Matrix with P-value\n")
    mtx <- matrixlist_pvalue[[i]]
    print (mtx)
    outfile <- sprintf("%s/%s.merged.withpvalue.txt", outpath, stemnamevector[i])
    write.table(mtx, outfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

