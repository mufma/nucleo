library(rtracklayer)

# Load table from txt file
tabData <- read.table("/Users/maratmufteev/Desktop/CSC2431/Genome/chr1genes.txt")

# Pointer to the bigWig file
bwf <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Computations/GSM1368906_TW245_K562_MeDIP.bigWig")

# Choose number of genes
#max <- length(tabData[4])
max <- 5
# Get starts and ends of genes in chr1
geneStart <- as.numeric(levels(tabData[c(2:max),4]))[tabData[c(2:max),4]]
geneEnd <- as.numeric(levels(tabData[c(2:max),5]))[tabData[c(2:max),5]]


for (num in 1:length(geneStart)) {
	# Specify gene region
	st <- geneStart[num]
	en <- geneEnd[num]
	# Choose region of chromosome
	selection <- GRanges(seqnames = "chr1", 		ranges = IRanges(start =  st, end = 		en))
	# Load data from the selection region in 	bigWig
	out <- import(bwf, selection=selection)
	# Open device to plot data into pdf
	pdf(file=sprintf("/Users/maratmufteev/Desktop/CSC2431/Genome/Pictures/picture%g.pdf", num))
	plot(score(out), main="Methylation sites", xlab ="Nucleotide, (bp)", ylab="Score", col="blue", type='l', cex=0.5)
	# Close pdf device
	graphics.off()
}




