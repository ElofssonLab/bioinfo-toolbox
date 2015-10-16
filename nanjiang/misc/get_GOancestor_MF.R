# Description:
#   Get ancenstors of GO term for ontology "Molecular function", given a list
#   of GO IDs. This code is modified from get_GOancestor.R



# bugs still exist in GO tree. e.g. for GOID = "GO:0004129", even the
# geneontology web server give out the wrong tree.

# For ontology "Molecular function", the first level keeps only those point to
# GO:0003674

# ChangeLog 2013-09-19 /*{{{*/
#   output to stderr should be done with the following syntax
#   sink(stderr()); cat("message\n"); sink()
#   From 2013-09-19, this code is maintained under /data3/wk/MPTopo/src
# /*}}}*/
usage <- "
Usage: Rscript  get_GOancestor_MF.R GOIDFile > OUTFILE

Created 2013-09-11, updated 2013-09-19, Nanjiang Shu
"

library(GO.db)

GO_MF <- "GO:0003674"

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

GOIDListfile <- args[1]

goidlist <- read.table(GOIDListfile, header=F)
goidlist = t(goidlist)
#cat ("Number of GO idlist is:", length(goidlist), "\n\n")
# GOMFANCESTOR[[term]]
# goidlist

numID = length(goidlist)
for (i in 1:numID){
#    cat (i, "\n")
     goid <- goidlist[i]
     lst <- list()

     cntlevel <- 1
     lst[[cntlevel]] <- c(goid)
     #cat("level", cntlevel, lst[[cntlevel]], "\n") #debug

     addedGOIDSet = c(goid)
     prtlist <- GOMFPARENTS[[goid]]
     prtlist <- setdiff(prtlist, addedGOIDSet)

     if (length(prtlist) < 1){
        sink(stderr());
        cat(goid, ";", "NO PARENTS\n")
        sink()
        next
     }

     cntlevel <- cntlevel + 1
     lst[[cntlevel]] <- prtlist

     addedGOIDSet <- union(addedGOIDSet, prtlist)
     nprtlist <- length(prtlist)

# Search iteratively the parent GOTERM until "all" is reached
     while((length(prtlist) >= 1) && (is.element("all", prtlist) == FALSE)){
# Note: this solved the bad ancestor tree problem, stop if all is reached
        unionList <- c()
        for (tp in prtlist){
            newprtlist <- GOMFPARENTS[[tp]]
            newprtlist <- setdiff(newprtlist, addedGOIDSet)
            unionList <- union(newprtlist, unionList)
            addedGOIDSet <- union(addedGOIDSet, newprtlist)
        }
        prtlist <- unionList
        cntlevel <- cntlevel + 1
        lst[[cntlevel]] <- prtlist
     }

# check the list
     li1 <- lst[[cntlevel]]
     li2 <- lst[[cntlevel-1]]
     li3 <- lst[[cntlevel-2]]
     if (is.element("all", li1)  && is.element(GO_MF, li2)){
        newli <- c()
        for (goid in li3){
            newprtlist <- GOMFPARENTS[[goid]]
            if (is.element(GO_MF, newprtlist)){
                newli <- c(newli, goid)
            }
        }
        lst[[cntlevel]] <- c("all")
        lst[[cntlevel-1]] <- c(GO_MF)
        lst[[cntlevel-2]] <- newli
        cat(lst[[1]])
        for (j in 2:cntlevel){
            cat(";", lst[[j]])
        }
        cat ("\n")
     }else{
        sink(stderr());
        cat(goid, ";", "BAD ANCESTOR TREE\n")
        print(lst)
        sink()
     }
}

