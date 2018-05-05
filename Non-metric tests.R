# NON-METRIC DATA ANALYSIS

#Import dataset X - keep quotations

NMData <- data.frame(X)

# Relevant Packages
library("dplyr", lib.loc="C:/Program Files/R/R-3.3.0/library")
library("psych", lib.loc="C:/Program Files/R/R-3.3.0/library")

# To view multiple variables in order:
# NMSex <- arrange(NMData, Genotype_Treatment, Sex) # Data arranges by strain then sex
# NMDam <- arrange(NMData, Genotype_Treatment, Mother) # Data arranged by strain then mother

# For all levels
sapply(NMData[,11:ncol(NMData)], unique, na.rm=TRUE)

# Sex and Dam in table
Grouped_Sex <- group_by(NMData, Genotype_Treatment, Sex)
Sex_by_Group <- summarise(Grouped_Sex, count = n())

Grouped_Dam <- group_by(NMData, Genotype_Treatment, Mother)
Dam_by_Group <- summarise(Grouped_Dam, count = n())

# Strain
Grouped_Strain <- group_by(NMData, Genotype_Treatment)
Strain <- summarise(Grouped_Strain, count = n())

# 2-way frequency Tables

Nas_fus <- table(NMData$Genotype_Treatment, NMData$Nasal.fusion..not.B.)
summary(Nas_fus) # For chi-squared


Nas <- table(NMData$Genotype_Treatment, NMData$Nasal.fusion..not.B.)
SF <- table(NMData$Genotype_Treatment, NMData$Squamosal.frontal.fusion)
SP <- table(NMData$Genotype_Treatment, NMData$Squamosal.parietal.fusion)
TH <- table(NMData$Genotype_Treatment, NMData$Post.tympanic.hook.of.squamosal.parietal.fusion)
BB <- table(NMData$Genotype_Treatment, NMData$Basisphenoid.basioccipital.fusion..not.B.)
BP <- table(NMData$Genotype_Treatment, NMData$Basisphenoid.presphenoid.fusion..not.B.)
DF <- table(NMData$Genotype_Treatment, NMData$Dorsal.frontal.fusion..not.bilateral.)
OP <- table(NMData$Genotype_Treatment, NMData$Occipital.pteriotic.fusion)




