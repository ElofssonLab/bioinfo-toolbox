library(bio3d)
library(Biostrings)

args <- (commandArgs(TRUE))
args

if(length(args) != 3){
    #stop("Usage: Rscript plot_contact_map.R <PDB_ID> <DOMAIN_ID> <DOMAIN_LEN> <DI_FILE> <PFAM_OFFSET> <PDB_OFFSET>")
    stop("Usage: Rscript plot_contact_map.R <DOMAIN_ID> <DOMAIN_LEN> <DI_FILE>")
}

DOMAIN_ID <- args[1]
DOMAIN_LEN <- as.numeric(args[2])
DI_FILE <- args[3]


# read predicted contacts
di_scores <- read.table(DI_FILE,sep=",")
di_scores_sorted <- di_scores[order(-di_scores[,3]),]
di_scores_sorted <- di_scores_sorted[which(abs(di_scores_sorted[,1] - di_scores_sorted[,2]) > 1),]
contacts <- di_scores_sorted[1:(DOMAIN_LEN * 3),]

# plot contact map with CA distance cutoff = 8 Å 
pdf(paste(DOMAIN_ID, "_ContactMap_nopdb.pdf", sep=""))
#plot.dmat(ref.cont, col=c("white","grey"), nlevels=2, key=FALSE, flip=FALSE, pch=15, xlab="residue index", ylab="residue index", main=paste(DOMAIN_ID, " - ", PDB_ID, sep=""))
plot(contacts[,2], contacts[,1], pch=3, col="red", cex=0.8, xlab="residue index", ylab="residue index", main=DOMAIN_ID)
dev.off()

# compute TP and FP
#contact_map <- matrix(0,nrow=143,ncol=143)
#for(i in 1:nrow(contacts)){contact_map[contacts[i,1],contacts[i,2]] <- 1}
#
#mask <- which(t(contact_map) == 1) %in% which(ref.cont==1)
#
#TP <- length(which(mask))/nrow(contacts)
#FP <- length(which(!mask))/nrow(contacts)
#
#TP
#FP
