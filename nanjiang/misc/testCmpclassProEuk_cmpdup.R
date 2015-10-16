# Filename: testCmpclassProEuk_cmpdup.R
# Description:
# calculate the p value of cmpclass fraction between prokaryota and eukaryota

# Created 2013-09-16, updated 2013-09-16, Nanjiang Shu 
# first input is prokar
# second input is eukar

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

args <- commandArgs(trailingOnly = TRUE)

numFile <- length(args)
cat ("numFile =", numFile, "\n")
filevector <- c()
for (i in 1:numFile){
    filevector <- c(filevector, args[i])
}
tagvector <- c("prokar", "eukar")

cat ("File list: ", filevector, "\n")
cat ("Tag list: ", tagvector, "\n")

# check whether the last file is all.gofreq.txt
if (numFile != 2){
    stop("numFile != 2")
}

tablelist <- list()
for (i in 1:numFile){
    tablelist[[i]] <- read.table(filevector[i], header=F)
}

prokar <- tablelist[[1]]
eukar <- tablelist[[2]]

prokar_count <- prokar
eukar_count <- eukar

hypolist <- c("two.sided", "greater", "less")
num_hypo <- length(hypolist)
num_row <- nrow(prokar)

pvaluelist <- list()
taglist <- c("IDT",  "INV", "DUP", "TM2GAP", "TM2SEQ",     "TM2SP",
"Mixed", "DiffAll", "Diff without MT2SEQ and TM2SP")

for (icol in 3:11){
    pvaluevector <- c()

    if (icol <= 9){
        cnt_prokar <- as.integer(prokar[[icol]]*prokar$V11/100)
        cnt_eukar <-  as.integer(eukar[[icol]]*eukar$V11/100)   
    }else if (icol == 10){ #DiffAll
        cnt_prokar <- rep(0,num_row)
        cnt_eukar <- rep(0,num_row)
        for (j in 4:9){
            cnt_prokar <- cnt_prokar + prokar[[j]]
            cnt_eukar <-  cnt_eukar + eukar[[j]]
        }
        cnt_prokar <- as.integer(cnt_prokar*prokar$V11/100)
        cnt_eukar <- as.integer(cnt_eukar*eukar$V11/100)
    }else if (icol == 11){ #DiffAllWithout
        cnt_prokar <- rep(0,num_row)
        cnt_eukar <- rep(0,num_row)
        for (j in c(4,5,6,9)){
            cnt_prokar <- cnt_prokar + prokar[[j]]
            cnt_eukar <-  cnt_eukar + eukar[[j]]
        }
        cnt_prokar <- as.integer(cnt_prokar*prokar$V11/100)
        cnt_eukar <- as.integer(cnt_eukar*eukar$V11/100)
    }

    total_prokar <- prokar$V11
    total_eukar <- eukar$V11

    for (j in 1:num_row){
        v1 <- c(cnt_prokar[j], cnt_eukar[j])
        v2 <- c(total_prokar[j], total_eukar[j])
        for (ihypo in 1:num_hypo){
            rst <- prop.test(v1, v2, alternative = hypolist[ihypo])
            pvaluevector <- c(pvaluevector, rst$p.value)
            cat ("v1=",v1, "v2=",v2, "hypo=", hypolist[ihypo],"\n")
            print (rst$observed)
            print (rst)
            warnings()
        }
    }
    #mtx <- matrix(pvaluevector, nrow=num_row, ncol=num_hypo, dimnames=list(prokar$V2, hypolist))
    mtx <- matrix(pvaluevector, nrow=num_hypo, ncol=num_row, dimnames=list(hypolist, prokar$V2))
    mtx <- t(mtx)
    pvaluelist[[icol]] <- mtx
}


for (icol in 3:11){
    mtx <- pvaluelist[[icol]]
    cat("\n")
    cat (taglist[icol-2], "\n")
    print (mtx)
}

