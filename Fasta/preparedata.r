# Load library to work with sequence data
library(seqinr)

# Load data about sequences
faData <- read.fasta("/Users/maratmufteev/Desktop/CSC2431/Data/Sequence/chr1.fa")

# Create vector of sequence data
seqFull <- faData$chr1

# Load peaks as dataframes
cat("Loading peaks data for nucleosomes and methylation ...", "\n")
peakN <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/peakNuc%d.txt","chr1",1), sep="\t")
peakM <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/Peaks/%s/peakMet%d.txt","chr1",1), sep="\t")

# Make peaks vectors out of dataframes
cat("Creating peaks vectors out of dataframes ...", "\n")
peakNuc <- as.vector(peakN[,])
peakMet <- as.vector(peakM[,])

cat("Removing peaks dataframes","\n")
rm(peakN,peakM)
gc()
gc()

# Load smoothed data as dataframes
cat("Loading peaks data for nucleosomes and methylation ...", "\n")
smN <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/smNuc%d.txt","chr1",1), sep="\t")
smM <- read.table(sprintf("/Users/maratmufteev/Desktop/CSC2431/Data/VectorSignals/%s/smMet%d.txt","chr1",1), sep="\t")

# Make smooth vectors out of dataframes
cat("Creating peaks vectors out of dataframes ...", "\n")
smNuc <- as.vector(smN[,])
smMet <- as.vector(smM[,])

cat("Removing smoothed dataframes","\n")
rm(smN,smM)
gc()
gc()
