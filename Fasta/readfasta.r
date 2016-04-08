library(seqinr)

# Load sequence from fasta file
faData <- read.fasta("/Users/maratmufteev/Desktop/CSC2431/Data/Sequence/chr1.fa")

# Choose sequence from chr1
seqFull <- faData$chr1

# Slice sequence
seqData <- seqFull[1000000:1000300]

# Count number of nucleotides, dinucleotides, tinucleotides
oneNucl <- count(seqData,1)
twoNucl <- count(seqData,2)
threeNucl <- count(seqData,3)