#subsample size - how many discrete measurements to convert the frequency table back into
subsample<-10000
percentiles<-c(0,0.01,0.02,0.05,0.10,0.25,0.50,0.75,0.80,0.85,0.90,0.95,0.98,0.99,0.999,1.00)

#Set working directory to file directory with data
setwd("D:/ImageJ Macros/Fleiszig Lab/Bacterial velocity measurement/input2/Output")

#Load the frequency table
freqMatrix<-read.table("Hist Matrix.csv", header = TRUE, sep = ",")

#Create a vector of the bin values
binVector<-freqMatrix[,"Speed"]

#Create a value and label vector from the freq table to perform a Kruskal-Wallis test
if(0){
	freqVector <- vector(mode="numeric", length=0) #Initialize empty vector
	labelVector <- vector(mode="numeric", length=0) #Initialize empty vector
	for(a in 2:ncol(freqMatrix)){
		freqTemp <- rep(binVector, freqMatrix[,a]) #Repeat each value for the number of times measured in the freq
		labelTemp <- rep(colnames(freqMatrix)[a], length(freqTemp)); #Generate a label vector oequal length
		
		#Concenate the temporary vectors onto the final vectors
		freqVector <- c(freqVector, freqTemp)
		labelVector <- c(labelVector, labelTemp)
	}
	
	labelVector <- as.factor(labelVector) #Convert character vector to label (group) type vector
	kruskal.test(freqVector, labelVector)
}


#Create a distance matrix off of the KS D-statistic
distMatrix <- matrix(, nrow = ncol(freqMatrix)-2, ncol = ncol(freqMatrix)-2) #Initialize a matrix that is nSamples x nSamples
colnames(distMatrix) <- colnames(freqMatrix)[3:ncol(freqMatrix)]
rownames(distMatrix) <- colnames(freqMatrix)[3:ncol(freqMatrix)]

#Also create matrix of percentiles
quantMatrix <- matrix(, nrow = length(percentiles)+1, ncol = ncol(freqMatrix)-2) #Initialize a matrix that is nSamples x nSamples
colnames(quantMatrix) <- colnames(freqMatrix)[3:ncol(freqMatrix)]
rownames(quantMatrix) <- c(percentiles,"mean")

for(a in 3:ncol(freqMatrix)){
	for(b in 3:ncol(freqMatrix)){
		freqTemp1 <- rep(binVector, freqMatrix[,a]) #Repeat each value for the number of times measured in the freq
		freqTemp2 <- rep(binVector, freqMatrix[,b]) #Repeat each value for the number of times measured in the freq
		distMatrix[a-2,b-2]<-ks.test(freqTemp1, freqTemp2)$statistic
	}
	quantMatrix[,a-2]<-c(quantile(freqTemp1, probs = percentiles),mean(freqTemp1))
}

write.csv(distMatrix, file = "Distance Matrix.csv")
write.csv(quantMatrix, file = "Quantile Matrix.csv")

plot(hclust(as.dist(distMatrix)), method = "average", hang = -1) #Create UPGMA ("average") denfrogram with level branches (hang = -1)

#Make notched box plots of the 95 percentile
boxplotMatrix <- matrix(NA, nrow = 10, ncol = 3) #Initialize a matrix that is nConditions x max nMeasurements
colnames(boxplotMatrix)<-c("PAO1", "LatA", "Noco")
rowSel = 15
r = 1;
for(a in 1:ncol(quantMatrix)){
	if(a < 7){
		boxplotMatrix[r, 1]<-quantMatrix[rowSel,a]
		r<-r+1;
	}
	else if(a<15){
		if(a == 7) r = 1
		boxplotMatrix[r, 2]<-quantMatrix[rowSel,a]
		r<-r+1;
	}
	else{
		if(a == 15) r = 1
		boxplotMatrix[r, 3]<-quantMatrix[rowSel,a]
		r<-r+1;
	}
}
pvalueMatrix <- matrix(NA, nrow = 3, ncol = 3) #Initialize a matrix that is nConditions x max nMeasurements
for(a in 1:3){
	for(b in 1:3){
		pvalueMatrix[a,b]<-wilcox.test(boxplotMatrix[,a], boxplotMatrix[,b], alternative = "greater")$p.value
	}
}		
boxplot(boxplotMatrix, las = 2, names = colnames(boxplotMatrix), notch = FALSE, outline = FALSE, ps = 1, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1, bty="n")
mtext("90th percentile speed (μm/s)", side=2, line = 3, cex = 1)
pvalueMatrix





