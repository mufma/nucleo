"Functions smoothVector and findPeaks are from func.r"

# Load functions to convert GRanges to vectors and find peaks
source("/Users/maratmufteev/Desktop/CSC2431/Code/RawSig/func.r")

# Choose chromosome
chr <- "chr1"
part <- 1

getPeaks(chr, "dataNuc1.txt", "smNuc1.txt", "peakNuc1.txt", 50, "normal", 0.5)
getPeaks(chr, "dataMet1.txt", "smMet1.txt", "peakMet1.txt", 200, "box", 10)

