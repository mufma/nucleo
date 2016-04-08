




"Usage of ggplot2"
# Choose distancies close to origin
distance <- dist[abs(dist)<200]

# Convert to dataframe for ggplot
distdf <- data.frame(distance)

# Plot dataframe with ggplot
ggplot(distdf, aes(x=distance)) + geom_density(fill='red', alpha=0.5, adjust=0.5)

"Calculate distancies"
size=50000
dist <- calcDist(peakMet, peakNuc, size)
subpeakmet <- peakMet[1:size]

"Get statistics of nucleotide words"
st = peakMet[1]-width
en = peakMet[1]+width
piece = seqFull[st:en]
p1 = count(piece,1)
p2 = count(piece,2)
p3 = count(piece,3)
tot = c(p1,p2,p3,smMet[peakMet[1]])
x = tot
width=100
for (peak in 2:size){
	st = peakMet[peak]-width
	en = peakMet[peak]+width
	piece = seqFull[st:en]
	p1 = count(piece,1)
	p2 = count(piece,2)
	p3 = count(piece,3)
	tot = c(p1,p2,p3,smMet[peakMet[peak]])
	x = rbind(x,tot)
}

"Plot densities on the same figures"
par(mfrow=c(2,2))
piece = seq(from=1,to=50000,length.out=5)
for (ind in 1:4) {
	st = piece[ind]
	en = piece[ind+1]
	D = dist[st:en]
	d<-density(D[abs(D)<100],bw=4)
	if (ind==3) name = plot(d, main=sprintf("Piece %d",ind), xlab="Distance(bp)")
	if (ind==4) name = plot(d, main=sprintf("Piece %d",ind), xlab="Distance(bp)")
	if (ind==1) name = plot(d, main=sprintf("Piece %d",ind))
	if (ind==2) name = plot(d, main=sprintf("Piece %d",ind))	
}


"Multiple plots at single figure. Requires multiplot function."
p1 <- qplot((peakMet[1]-200):(peakMet[1]+200),smMet[(peakMet[1]-200):(peakMet[1]+200)],geom='smooth',xlab='Genomic Location (bp)',ylab='Signal',main='Methylation')
p2 <- qplot((peakMet[452]-200):(peakMet[452]+200),smMet[(peakMet[452]-200):(peakMet[452]+200)],geom='smooth',xlab='Genomic Location (bp)',ylab='Signal',main='Methylation')
p3 <- qplot((peakMet[5636]-200):(peakMet[5636]+200),smMet[(peakMet[5636]-200):(peakMet[5636]+200)],geom='smooth',xlab='Genomic Location (bp)',ylab='Signal',main='Methylation')
multiplot(p1,p2,p3)


"Load data for methylation or nucleosome"
cat("Loading data for nucleosomes and methylation ...", "\n")
vecM <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/dataMet%d.txt",chr,part), sep="\t")
cat("Creating vectors out of dataframes ...", "\n")
vecMet <- as.vector(vecM[,])



"Old workflow FindPeaks"
# Load data as dataframes
cat("Loading data for nucleosomes and methylation ...", "\n")
vecN <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/dataNuc%d.txt",chr,part), sep="\t")
vecM <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/dataMet%d.txt",chr,part), sep="\t")

# Make vectors out of dataframes
cat("Creating vectors out of dataframes ...", "\n")
vecNuc <- as.vector(vecN[,])
vecMet <- as.vector(vecM[,])

# Remove dataframes
rm(vecN, vecM)

# Smooth data in vectors
cat("Smoothing data in vectors ...", "\n")
smNuc <- smoothVector(vecNuc, bandwidth=50, kernel="normal")
smMet <- smoothVector(vecMet, bandwidth=200, kernel="box")

# Remove data vectors
rm(vecNuc, vecMet)

# Save smoothed vectors in txt files
cat("Saving smoothed vectors to the txt files ...", "\n")
write.table(smNuc, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/smNuc%d.txt",chr,part), sep="\t")
write.table(smMet, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/smMet%d.txt",chr,part), sep="\t")

# Find peaks in vectors
cat("Searching for peaks ...", "\n")
peakNuc <- findPeaks(smNuc, peakHeight=0.5)
peakMet <- findPeaks(smMet, peakHeight=10)

# Remove smoothed vectors
rm(smNuc, smMet)

# Save peaks in txt files
cat("Saving peaks to the txt files ...", "\n")
write.table(peakNuc, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/peakNuc%d.txt",chr,part), sep="\t")
write.table(peakMet, sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/peakMet%d.txt",chr,part), sep="\t")

# Remove peaks data
rm(peakNuc, peakMet)
