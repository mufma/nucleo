# Load table from txt file
tabData <- read.table("/Users/maratmufteev/Desktop/CSC2431/Genome/chr1genes.txt")

# Create numeric type out of gene start data
max <- length(tabData[[4]])
levData <- tabData[c(2:max),4]
stGene <- as.numeric(levels(levData))[levData]

# Create numeric type out of gene end data
levData <- tabData[c(2:max),5]
enGene <- as.numeric(levels(levData))[levData]

# Don't forget the ranges of loaded meth and nucl data
num = 4
en <- enGene[num]
st <- stGene[num]
smM <- smMet[st:en]
smN <- smNuc[st:en]
genLoc <- st:en
dfG <- data.frame(genLoc, smM, smN)

# Show methylation and nucleosomes on the gene
ggplot(dfG, aes(genLoc)) + geom_line(aes(y=smM), colour="red") + geom_line(aes(y=smN), colour="green")
