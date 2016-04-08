# Package for bigWig filetype
library(rtracklayer)

# Pointer to the bigWig file
bwf <- BigWigFile("/Users/maratmufteev/Desktop/CSC2431/Computations/GSM920558_hg19_wgEncodeSydhNsomeGm12878Sig.bigWig")

# Choose chromosome
chr <- "chr1"

# Choose region of chromosome
max <- seqlengths(seqinfo(bwf)[chr])
selection <- GRanges(seqnames=chr, ranges=IRanges(start=1, end=max))

# Load data from the selection region in bigWig
track <- import(bwf, selection=selection)

# Get lengths and values of platos
x <- rle(score(track))

# Get peaks of the data
peaks <- ranges(track[which(rep(diff(sign(diff(c(-Inf, x$values, -Inf)))) == -2, times=x$lengths))])

# Free memory
rm(track, x, selection, max, bwf)

# Get midpoints
midpoints <- (start(peaks) + end(peaks))/2

# Get distancies between nucleosomes
dist <- diff(midpoints)

# Get histogram data of distancies in logscale
mydata_hist <- hist(log(diff(midpoints)))

# Create dataframe with distance data
df <- data.frame(chromosome = factor( rep(c(chr), each=length(dist)) ), distance = dist)

# Show and save density plot
pdf(file=sprintf("/Users/maratmufteev/Desktop/CSC2431/Computations/Pictures/picture%g.pdf", num))
ggplot(df, aes(x=distance, fill=chromosome)) + geom_density(alpha=.8) + scale_x_log10()

# Close figure
graphics.off()

