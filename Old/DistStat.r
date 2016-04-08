"
Calculate distribution of distances between 
nucleosome and methylation peaks.

"
# Package for bigWig filetype
library(rtracklayer)
# Load functions to convert GRanges to vectors and find peaks
source("/Users/maratmufteev/Desktop/CSC2431/Computations/MyFunc.r")

# Choose chromosome
chr <- "chr1"

# Pointer to the bigWig file with nucleosome data
bwfNuc <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Computations/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bigWig")

# Pointer to the bigWig file with methylation data
bwfMet <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Computations/GSM1368907_TW246_Gm12878_MeDIP.bigWig")

# Choose region of chromosome
max <- seqlengths(seqinfo(bwfNuc)[chr])
selection <- GRanges(seqnames=chr, ranges=IRanges(start=1, end=max))

# Load data from the selection region in bigWig
trackNuc <- import(bwfNuc, selection=selection)
trackMet <- import(bwfMet, selection=selection)

# Convert GRanges to vectors
vecNuc <- toVector(trackNuc)
vecMet <- toVector(trackMet)

# Save data as vectors
write.table(vecNuc, "/Users/maratmufteev/Desktop/CSC2431/Computations/dataNuc.txt", sep="\t")
write.table(vecMet, "/Users/maratmufteev/Desktop/CSC2431/Computations/dataMet.txt", sep="\t")

# Load data as vectors
vecNuc <- read.table("/Users/maratmufteev/Desktop/CSC2431/Computations/dataNuc.txt", sep="\t")
vecMet <- read.table("/Users/maratmufteev/Desktop/CSC2431/Computations/dataMet.txt", sep="\t")

# Smooth data in vectors
smNuc <- smoothVector(vecNuc, bandwidth=50, kernel="normal")
smMet <- smoothVector(vecMet, bandwidth=200, kernel="box")

# Find peaks in vectors
peakNuc <- findPeaks(smNuc, peakHeight=0.8)
peakMet <- findPeaks(smMet, peakHeight=100)








