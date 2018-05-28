 #!/usr/bin/Rscript

#####################################################################
####     Biolog Data analysis for the experimental evolution       ##
#####################################################################

#Libraries

library("readxl")
library("tools")
library("stringr")
library("ggplot2")
library("ggfortify")
library("FactoMineR")
library("stats")
library("dendextend")
library("colorspace")

#Input

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Load the biolog log file which contains all the strains, names,
# nicknames, date plate was read, length of incubation time before
# plate read, temperature treatment of strain, and if run was successful or not
logFile=read_xlsx("biologLogFile.xlsx", col_names = TRUE)

# Create list of all strains that biolog was done for,
# thus has an .xls file containing the biolog plate reader data
biologFiles=list.files(pattern = "*.xls")
biologFiles=biologFiles[biologFiles != "biologLogFile.xlsx"]

# Read all biolog files for all the strains
biologFileNames=lapply(biologFiles,file_path_sans_ext)

# Create a data frame containing all the biolog data for all the strains
biologVectorsDataFrame2=as.data.frame(matrix(ncol=length(biologFiles),nrow=99))
# Set columns of data frame as corresponding filenames
names(biologVectorsDataFrame2)=biologFileNames
row.names(biologVectorsDataFrame2)=c("Strain","Treatment","Condition",sprintf("Well%d",1:96))

# Set strain name in the data frame row "Strain"
strainNames=str_split_fixed(biologFileNames,"-",3)[,1]
strainTypes=logFile$`Strain Type`[match(strainNames,logFile$Nickname)]
biologVectorsDataFrame2[1,]=strainTypes

# set strain treatments in the data frame row "Treatment"
strainTreatments=logFile$Treatment[match(strainNames,logFile$Nickname)]
biologVectorsDataFrame2[2,]=strainTreatments

# set conditions (temperatures of strains) in the data frame row "Condition"
biologVectorsDataFrame2[3,]=as.integer(str_split_fixed(biologFileNames,"-",3)[,2])
# set control strains condition to 37 instead of 606 (ex: control1-606-37-1.xls filename and we want 37 and not 606)
biologVectorsDataFrame2[3,biologVectorsDataFrame2[3,]==606] <- as.integer(37)

for(i in 1:length(biologFiles)) {
  # Import biologFile as xls file, and convert to matrix form then convert 
  # matrix to vector and add this vector to appropriate column in data frame
  biologVectorsDataFrame2[4:99,i]=as.vector(data.matrix(read_excel(biologFiles[i]))[1:8,2:13])
}

normalizedBiologDataFrame=biologVectorsDataFrame2[c(1:3,5:75,77:99),] 
# for loop to normalize all strains and their replicates (number of columns
# in biologVectorsDataFrame2 is equal to number of strains)
for(i in 1:length(biologVectorsDataFrame2)) {
  # normalize carbon/nitrogen assays' entries 5:75 in biologVectorsDataFrame2 (entries 4:74 in BiologDF2)
  # to negative control (entry 4 in biologVectorsDataFrame2)
  normalizedBiologDataFrame[4:74,i] <- as.numeric(normalizedBiologDataFrame[4:74,i]) - as.numeric(biologVectorsDataFrame2[4,i])
  
  # normalize chemical assays' wells 74:96 on biolog plates (entries 77:99 in biologVectorsDataFrame2) 
  # to positive control (well 73; entry 75 in biologVectorsDataFrame2)
  normalizedBiologDataFrame[75:97,i] <- as.numeric(biologVectorsDataFrame2[76,i]) - as.numeric(normalizedBiologDataFrame[75:97,i])
 
}


##### PCA of all biolog data ######
BiologDF2mat <- t(data.matrix(normalizedBiologDataFrame[4:97,])) #if transpose a data.frame containing any characters
    #then all numeric values will also be converted to strings; so must convert to data.matrix first
pcaBiologDF2 <- princomp(BiologDF2mat)
# Plot PCA results with strains colored by Treatment
dev.new()
autoplot(pcaBiologDF2, data = t(normalizedBiologDataFrame), colour = "Treatment", shape = FALSE, label = TRUE, label.size = 3, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, main = "PCA plot of all Biolog Data colored according to Treatment")
# Plot PCA results with strains colored by Type
dev.new()
autoplot(pcaBiologDF2, data = t(normalizedBiologDataFrame), colour = "Strain", shape = FALSE, label = TRUE, label.size = 3, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, main = "PCA plot of all Biolog Data colored according to Type")


##### Hierarchical clustering of all biolog data #######
# Create dendrogram of biolog data for all strains using euclidean distance and ward method for clustering
dend <- as.dendrogram(hclust(dist(BiologDF2mat, method = "euclidean"), method = "ward.D"))
# Set labels of dendrogram to strain names
labels(dend) <- names(normalizedBiologDataFrame)[order.dendrogram(dend)]
# Color dendrogram labels and branches according to treatment of strains
treatments <- factor(strainTreatments)
treatmentsColorsVec <- colorspace::rainbow_hcl(length(unique(treatments)), c = 70, l = 50)
colorTreatments <- treatmentsColorsVec[treatments]
labels_colors(dend) <- colorTreatments[order.dendrogram(dend)]
dend <- color_branches(dend, col = colorTreatments[order.dendrogram(dend)])
# Plot dendrogram with colored labels and branches matching Treatment
dev.new()
plot(dend)
title(main = "Dendrogram of all Biolog Data colored according to Treatment")

# Create dendrogram of biolog data for all strains using euclidean distance and ward method for clustering
dendTypes <- as.dendrogram(hclust(dist(BiologDF2mat, method = "euclidean"), method = "ward.D"))
# Set labels of dendrogram to strain names
labels(dendTypes) <- names(normalizedBiologDataFrame)[order.dendrogram(dendTypes)]
# Color dendrogram labels and branches according to Strain Type
types <- factor(strainTypes)
typesColorsVec <- colorspace::rainbow_hcl(length(unique(types)), c = 70, l = 50)
colorTypes <- typesColorsVec[types]
labels_colors(dendTypes) <- colorTypes[order.dendrogram(dendTypes)]
dendTypes <- color_branches(dendTypes, col = colorTypes[order.dendrogram(dendTypes)])
# Plot dendrogram with colored labels and branches matching Strain Type
dev.new()
plot(dendTypes)
title(main = "Dendrogram of all Biolog Data colored according to Strain Type")
