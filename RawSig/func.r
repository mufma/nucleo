# Load library for for distance functions
library(binhf)

"Convert GRanges data to the vector format"
toVector <- function(track){ 
	dataVector = rep(0, max(end(track)))
	# Get start,end,width and score of IRanges in data track
	st = start(track)
	en = end(track)
	wd = width(track)
	sc = score(track)
	for (i in 1:length(track)) {
		# Get indexes of ith range
		vec = st[i]:en[i]
		# Get score of the plateu, repeat according to width
		data = rep(sc[i], wd[i])
		# Save to vector
		dataVector[vec] = data
	}
	return(dataVector)
}

"Pipeline to vectorize data"
vectorizeData <- function(bwf, selection, chr, savefile){
	
	# Load nucleosome data from the selection region in bigWig
	cat("Loading data from bigWig files ...", "\n")
	track <- import(bwf, selection=selection)
	cat("Converting GRanges to vectors ...", "\n")
	vec <- toVector(track)
	
	# Remove nucleosome tracks
	cat("Removing tracks ...", "\n")
	rm(track)

	# Save data
	cat("Saving data to the txt files ...", "\n")
	write.table(vec, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/%s", chr, savefile), sep="\t")

	# Removing vectors
	cat("Removing vectors", "\n")
	rm(vec)
	
	return(0)
}

"Smooth data and save to vector"
smoothVector <- function(dataVector, bandwidth, kernel){
	# Get range of locations on the genome
	genLoc <- 1:length(dataVector)
	# Smooth vector with signals
	dataList <- ksmooth(genLoc,dataVector, kernel=kernel, bandwidth=bandwidth)
	# Convert list of signals to vector of signals
	dataSmooth <- as.vector(dataList[[2]])
	return(dataSmooth)
}

"Find peaks in the vector"
findPeaks <- function(dataSmooth, peakHeight){
	# Use smoothed vector with signals
	x <- dataSmooth
	# Find levels of constant value, calculate lengths
	r <- rle(x)
	# Find peaks
	peaks <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
	# Select large enough peaks
	peaks <- peaks[x[peaks]>peakHeight]
	return(peaks)
}

"Primary algorithm to calculate distancies"
calcDist <- function(peakMet, peakNuc, numMet) {
	st = proc.time()
	# Number of methylation peaks to consider
	num = numMet
	dist = c()
	for (i in 1:num) {
		sample = peakNuc[peakNuc<(peakMet[i]+10000) & peakNuc>(peakMet[i]-10000)]
		temp = rep(0,length(sample))
		for (j in 1:length(sample)){
			# Calculate distance between meth. peak i and nucleosome j
			temp[j] = peakMet[i] - sample[j]
			}
		# Choose minimal absolute distance nucleosome
		dist = c(dist,temp[which.min(abs(temp))])
	}
	en = proc.time()
	cat("Time elapsed:", en-st, "\n")
	return(dist)
}

"First algorithm to find distancies between peaks"
calcDistOne <- function(peakMet, peakNuc, numMet) {
	st = proc.time()
	# Number of methylation peaks to consider
	num = numMet
	dist = c()
	for (i in 1:num) {
		temp = rep(0,length(peakNuc))
		for (j in 1:length(peakNuc)){
			# Calculate distance between meth. peak i and nucleosome j
			temp[j] = peakMet[i] - peakNuc[j]
			}
		# Choose minimal absolute distance nucleosome
		dist = c(dist,temp[which.min(abs(temp))])
	}
	en = proc.time()
	cat("Time elapsed:", en-st, "\n")
	return(dist)
}

"Second algorithm to find distancies between peaks"
calcDistTwo <- function(peakMet, peakNuc, elements) {
	st = proc.time()
	num = min(length(peakMet), length(peakNuc))
	num = elements
	v1 = peakMet[1:num]
	v2 = peakNuc[1:num]
	dM = c()
	dv = v1 - v2
	dM = cbind(dM, dv)
	for (i in 1:num) {
		v2 <- shift(v2, 1)
		dv = v1 - v2
		dM = cbind(dM, dv)
	}
	dist = c()
	for (i in 1:num) {
		val = dM[,i]
		dist = c(dist, val[which.min(abs(val))])
	}
	en = proc.time()
	cat("Time elapsed:", en-st, "\n")
	return(dist)
}

"Get distancies to files"
getDist <- function(chr, part, numMetPeaks){
	
	# Load data as dataframes
	cat("Loading data for nucleosomes and methylation ...", "\n")
	peakN <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/peakNuc%d.txt",chr,part), sep="\t")
	peakM <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/peakMet%d.txt",chr,part), sep="\t")

	# Make vectors out of dataframes
	cat("Creating vectors out of dataframes ...", "\n")
	peakNuc <- as.vector(peakN[,])
	peakMet <- as.vector(peakM[,])

	# Remove dataframes
	rm(peakN, peakM)

	cat("Calculating distancies ...", "\n")
	dist <- calcDist(peakMet, peakNuc, numMetPeaks)
	
	# Removing peak vectors
	rm(peakMet, peakNuc)
return(dist)
}

"Get peaks data to files"
getPeaks <- function(chr, filename, savesmVec, savePeak, bandwidth, kernel, peakThreshold){
	# Load data
	cat("Loading data for nucleosomes and methylation ...", "\n")
	vecdf <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/%s",chr,filename), sep="\t")
	
	# Make vectors out of dataframes
	cat("Creating vectors out of dataframes ...", "\n")
	vec <- as.vector(vecdf[,])
	
	# Remove dataframes
	rm(vecdf)
	
	# Smooth data in vectors
	cat("Smoothing data in vectors ...", "\n")
	smVec <- smoothVector(vec, bandwidth=bandwidth, kernel=kernel)
	
	# Remove data vectors
	rm(vec)
	
	# Save smoothed vectors in txt files
	cat("Saving smoothed vectors to the txt files ...", "\n")
	write.table(smVec, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/%s",chr,savesmVec), sep="\t")
	
	# Find peaks in vectors
	cat("Searching for peaks ...", "\n")
	peak <- findPeaks(smVec, peakHeight=peakThreshold)
	
	# Remove smoothed vectors
	rm(smVec)

	# Save peaks in txt files
	cat("Saving peaks to the txt files ...", "\n")
	write.table(peak, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/%s",chr,savePeak), sep="\t")
	
	# Remove peak vector
	rm(peak)
	
	return(0)
}

"Get vectors data to files"
getVectors <- function(chr, bwfMet, bwfNuc) {
	# Full length of the chromosome
	total <- seqlengths(seqinfo(bwfNuc)[chr])
	border <- seq(from=1,to=total,length.out=11)
	for (part in 1:10) {
		cat("Processing part:", part, "\n")
	
		# Select region of the chromosome
		min <- border[part]
		max <- border[part+1]
		
		# Specify region of the chromosome
		selection <- GRanges(seqnames=chr, ranges=IRanges(start=min, end=max))
	
		# Vectorize data
		vectorizeData(bwfNuc, selection, chr, sprintf("dataNuc%d.txt",part))
		vectorizeData(bwfMet, selection, chr, sprintf("dataMet%d.txt",part))
	}
	return(0)
}


