source("/Users/maratmufteev/Desktop/CSC2431/Code/RawSig/MyFunc.r")
library(rtracklayer)

# Choose chromosome
chr <- "chr1"

# Pointer to the bigWig file with nucleosome data
bwfNuc <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Data/Signals/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bigWig")

# Choose region of chromosome
max <- seqlengths(seqinfo(bwfNuc)[chr])
max <- 10000000
selection <- GRanges(seqnames=chr, ranges=IRanges(start=1, end=max))

# Load data from the selection region in bigWig
trackNuc <- import(bwfNuc, selection=selection)

st <- proc.time()
# Convert GRanges to vectors
vecNuc <- toVector(trackNuc)
en <- proc.time()

