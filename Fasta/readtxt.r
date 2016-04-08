# Load table from txt file
tabData <- read.table("/Users/maratmufteev/Desktop/CSC2431/Genome/chr1genes.txt")

# Slice rows
rowsData <- tabData[c(2:10),]

# Slice columns
colData <- tabData[1:3]

# Create numeric type out of column data
max <- length(tabData[[4]])
levData <- tabData[c(2:max),4]
vecData <- as.numeric(levels(levData))[levData]

# Show histogram
plot(density(vecData))
#hist(vecData)