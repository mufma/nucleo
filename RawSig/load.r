# Package for bigWig filetype
library(rtracklayer)
# Load functions to convert GRanges to vectors and find peaks
source("/Users/maratmufteev/Desktop/CSC2431/Code/RawSig/func.r")

# Choose chromosome
chr <- "chr1"

# Pointer to the bigWig file with nucleosome data
bwfNuc <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Data/Signals/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bigWig")

# Pointer to the bigWig file with methylation data
bwfMet <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Data/Signals/GSM1368907_TW246_Gm12878_MeDIP.bigWig")

# Convert GRanges to vectors and save to txt files
getVectors(chr, bwfMet, bwfNuc)

for (part in 1:10) {
	# Find peaks and save to txt files
	getPeaks(chr, sprintf("dataNuc%d.txt",part), sprintf("smNuc%d.txt",part), sprintf("peakNuc%d.txt",part), 50, "normal", 0.5)
	getPeaks(chr, sprintf("dataMet%d.txt",part), sprintf("smMet%d.txt",part), sprintf("peakMet%d.txt",part), 200, "box", 10)
}
