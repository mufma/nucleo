# Package for bigWig filetype
library(rtracklayer)

# Portion of data
selection <- GRanges(seqnames = "chr1", ranges = IRanges(start =  1, end = 14000))

# Pointer to the bigWig file
bwf <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Computations/GSM1368907_TW245_Gm12878_MeDIP.bigWig")

#bwf <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Computations/GSM920557_hg19_wgEncodeSydhNsomeK562Sig.bigWig")

# Load data from the selection region
out <- import(bwf, selection=selection)

# Get scores and plot them
data = score(out)
dev.new(width=10, height=4)
plot(data, main="Methylation sites", xlab ="Nucleotide, (bp)", ylab="Score", col="blue", type='l', cex=0.5)
