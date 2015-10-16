# Get all GO terms
# usage: get_all_GOTERM.R

library(GO.db)

li <- as.list(GOTERM)
numID <- length(li)

for (i in 1:numID){
#    cat (i, "\n")
     goid <- GOID(li[[i]])
     term <- GOTERM[[goid]]
     termname <- attr(term, "Term")
     cat (goid , "\t", termname, "\n")
}

