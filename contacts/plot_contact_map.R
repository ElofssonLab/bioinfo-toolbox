library(bio3d)
library(Biostrings)

args <- (commandArgs(TRUE))
args

if(length(args) != 5){
    #stop("Usage: Rscript plot_contact_map.R <PDB_ID> <DOMAIN_ID> <DOMAIN_LEN> <DI_FILE> <PFAM_OFFSET> <PDB_OFFSET>")
    stop("Usage: Rscript plot_contact_map.R <PDB_ID> <DOMAIN_ID> <DOMAIN_LEN> <DI_FILE> <PFAM_OFFSET>")
}

PDB_ID <- args[1]
DOMAIN_ID <- args[2]
DOMAIN_LEN <- as.numeric(args[3])
DI_FILE <- args[4]
PFAM_OFFSET <- as.numeric(args[5])
#PDB_OFFSET <- as.numeric(args[6])

# read observed PDB contacts
ref.prot <- read.pdb(PDB_ID)
inds <- atom.select(ref.prot,"calpha")

# calculate PDB_OFFSET
# TODO: align domain profile to PDB sequence to get exact offset
ref.s1 <- ref.prot$seqres # get SEQRES sequence
ref.s2 <- seq.pdb(ref.prot) # get actual sequence from coordinates
#ref.s1 <- paste(aa321(ref.prot$seqres),collapse="") # get SEQRES sequence
#ref.s2 <- paste(seq.pdb(ref.prot),collapse="") # get actual sequence from coordinates
#ref.aln <- seqaln(seqbind(ref.s1,ref.s2))
#PDB_OFFSET <- min(which(ref.aln$ali[1,] == ref.aln$ali[2,])) - 1
PDB_OFFSET <- length(ref.s1) - length(ref.s2)
PDB_OFFSET <- 0


# read predicted contacts
di_scores <- read.table(DI_FILE,sep=",")
di_scores_sorted <- di_scores[order(-di_scores[,3]),]
di_scores_sorted <- di_scores_sorted[which(abs(di_scores_sorted[,1] - di_scores_sorted[,2]) > 1),]
contacts <- di_scores_sorted[1:(DOMAIN_LEN * 1.5),]

# plot contact map with CA distance cutoff = 8 Å 
pdf(paste(DOMAIN_ID, "_ContactMap_8_1.5.pdf", sep=""))
ref.cont <- cmap(ref.prot$xyz[inds$xyz], dcut=8, scut=0)
plot.dmat(ref.cont, col=c("white","grey"), nlevels=2, key=FALSE, flip=FALSE, pch=15, xlab="residue index", ylab="residue index", main=paste(DOMAIN_ID, " - ", PDB_ID, sep=""))
points(contacts[,2] + PFAM_OFFSET - PDB_OFFSET, contacts[,1] + PFAM_OFFSET - PDB_OFFSET, pch=3, col="red", cex=0.8)
rect(PFAM_OFFSET - PDB_OFFSET, PFAM_OFFSET - PDB_OFFSET, PFAM_OFFSET - PDB_OFFSET + DOMAIN_LEN, PFAM_OFFSET - PDB_OFFSET + DOMAIN_LEN)
dev.off()

# plot contact map with CA distance cutoff = 12 Å 
pdf(paste(DOMAIN_ID, "_ContactMap_12_1.5.pdf", sep=""))
ref.cont.12 <- cmap(ref.prot$xyz[inds$xyz], dcut=12, scut=0)
plot.dmat(ref.cont.12, col=c("white","grey"), nlevels=2, key=FALSE, flip=FALSE, pch=15, xlab="residue index", ylab="residue index", main=paste(DOMAIN_ID, " - ", PDB_ID, sep=""))
points(contacts[,2] + PFAM_OFFSET - PDB_OFFSET, contacts[,1]+ PFAM_OFFSET - PDB_OFFSET, pch=3, col="red", cex=0.8)
rect(PFAM_OFFSET - PDB_OFFSET, PFAM_OFFSET - PDB_OFFSET, PFAM_OFFSET - PDB_OFFSET + DOMAIN_LEN, PFAM_OFFSET - PDB_OFFSET + DOMAIN_LEN)
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
