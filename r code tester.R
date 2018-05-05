install.packages("installr")

library("geomorph", lib.loc="C:/Program Files/R/R-3.1.0/library") #load geomorph
#Import dataset: change as dataset changes
Mandibles <- read.delim("C:/Users/Morphometrics lab/Desktop/Mandibles.txt", quote="")
View(Mandibles)
names <- Mandibles[,1]
coords <- arrayspecs(Mandibles[,12:ncol(Mandibles)], 34, 3)#creates 3D array, p=no. landmarks, k=no. dimensions
dimnames(coords)[[3]] <- names
classifiers <- Mandibles[,1:11] #change latter column
dim(coords) # check to see that there are 3 numbers
land.pairs <- read.delim("C:/Users/Morphometrics lab/Desktop/land.pairs.txt", header=FALSE, quote="") #loads landmark pair matrix
View(land.pairs) #Table of landmark pair matrix
#General Procrustes Analysis
GPA <- gpagen(coords)
View(as.table(GPA$Csize))
write.table(as.table(GPA$Csize), "C:/CACZcent.txt", sep="\t") # export cent
#Bilateral symmetry
gdf <- geomorph.data.frame(coords = GPA$coords, ind = names, Cent = GPA$Csize)
Blat <- bilat.symmetry(coords, ind = ind, object.sym = TRUE, land.pairs = land.pairs, data = gdf)
Blat$shape.anova
plot(Blat, warpgrids = TRUE, mesh = NULL)
dim(Blat$asymm.shape)

#Splitting groups (initial)
two.coords <- two.d.array(GPA$coords) # converts symm, Procrustes coordinates into 2D
split <- split(as.data.frame(two.coords), classifiers$Genotype_Treatment, drop=FALSE) # Splits two coordinate dataset up into the groups
splitNames <- split(names, classifiers$Genotype_Treatment, drop = FALSE)
splitCent <- split(as.data.frame(GPA$Csize), classifiers$Genotype_Treatment, drop = FALSE)
#T.test on Csize: repeat per group
t.test(splitCent$`CAST-EiJ`, splitCent$CastxCzech) # tests centroid size- repeat per groups
#Size boxplots: http://www.cookbook-r.com/Graphs/
library("ggplot2", lib.loc="C:/Program Files/R/R-3.1.0/library")
GDF <- data.frame(ind=names, Centroid_Size = GPA$Csize, Strain = classifiers$Genotype_Treatment)
BP <- ggplot(GDF, aes(x= Strain, y = Centroid_Size, fill = Strain)) + geom_boxplot() + theme_bw() + guides(fill=FALSE)
BP  <- BP + scale_x_discrete(limits= c("CAST-EiJ", "WSB", "CZECHI-EiJ", "CastxCzech", "CastxWSB", "WSBxCzech", "WSBxCzech*WSB", "CastXWSB*Cast", "CastxCzech*Czech", "CastxCzech*F2", "CastxWSB*F2"))
ggsave(filename = "BP.tiff", plot = BP, dpi = 1500) #saves in documents
#can set colours with scale_fill_manual. 
#Bilat symmetry: repeat section for CACZ, CZECHI
gpaCast <- arrayspecs(split$`CAST-EiJ`, 34, 3) #into 3D
dimnames(gpaCast)[[3]] <- splitNames$`CAST-EiJ`
blatCast <- bilat.symmetry(gpaCast, ind = splitNames$`CAST-EiJ`, object.sym = TRUE, land.pairs = land.pairs)
summary(blatCast)

#PCA
library("ggplot2", lib.loc="C:/Program Files/R/R-3.1.0/library")
PCA <-plotTangentSpace(Blat$symm.shape, verbose = TRUE)
scores <- as.data.frame(PCA$pc.scores)
scores$Strain <- classifiers$Genotype_Treatment
PCplot <- ggplot(data = scores, aes(x=PC1, y = PC2, colour = Strain)) + geom_point()
PCP <- PCplot + stat_ellipse(geom="polygon", alpha=0.3, aes(fill=Strain))+ theme_bw()
PCP
ggsave(filename = "PCP.tiff", plot = PCP, dpi = 1500)
#For filled polygons, a little complex: need to split and recombine dataset with chull function and ddply
no.missing <- na.omit(scores) #chulling won't work with missing data
find_hull <- function(no.missing) no.missing[chull(no.missing$X1, no.missing$X2),] # adding hull
library("plyr", lib.loc="C:/Program Files/R/R-3.1.0/library") #need plyr
hulls <- ddply(no.missing, "Strain", find_hull)
PCplot
PCplot + geom_polygon(data=hulls, alpha=0.3, aes(fill=Strain))
#Thin Plate Splines for PCA
PC1mesh <- plotRefToTarget(PCA$pc.shapes$PC1min, PCA$pc.shapes$PC1max, method = "TPS")
PC2mesh <-  plotRefToTarget(PCA$pc.shapes$PC2min, PCA$pc.shapes$PC2max, method = "TPS")

#allomtery
All <- procD.allometry(Blat$symm.shape ~ gdf$Cent, f2 = NULL, f3 = NULL, logsz = TRUE, iter = 99)
#Use All$RSC for residuals
#Group allometries
gdf <- geomorph.data.frame(coords = Blat$symm.shape, ind = names, Cent = GPA$Csize, species = classifiers$Genotype_Treatment)
All <- procD.allometry(coords ~ Cent, f2 = ~species, f3 = NULL, logsz = TRUE, data = gdf, iter = 99)
#Dont use grouped allometries for regression
gdf <- data.frame(CAC = All$CAC, Strain = classifiers$Genotype_Treatment, Cent = GPA$Csize)
AllP <- ggplot(gdf, aes(x=Cent, y= CAC, color=Strain)) +geom_point() + geom_smooth(method=lm, se =FALSE, size=1) + theme_bw() + theme(legend.position = "bottom")
ggsave(filename = "AllP.tiff", plot = AllP, dpi = 1500)
#Using Regression Projection and log Cent
gdf <- data.frame(Regression = All$Reg.proj, Strain = classifiers$Genotype_Treatment, LogCentroidSize = log(GPA$Csize))
AllP <- ggplot(gdf, aes(x=LogCentroidSize, y= Regression, color=Strain)) +geom_point() + geom_smooth(method=lm, se =FALSE, size=1) + theme_bw() + theme(legend.position = "bottom")

#PCA on residuals
PCAresids <-plotTangentSpace(arrayspecs(All$RSC, 34, 3), verbose = TRUE)
View(PCAresids$pc.summary$importance)
scores <- as.data.frame(PCAresids$pc.scores)
scores$Strain <- classifiers$Genotype_Treatment
PCplot <- ggplot(data = scores, aes(x=PC1, y = PC2, colour = Strain)) + geom_point()
PCP <- PCplot + stat_ellipse(geom="polygon", alpha=0.3, aes(fill=Strain))+ theme_bw()
PCP
ggsave(filename = "PCP.tiff", plot = PCP, dpi = 1500)

#Morphological disparity without size
gdf <- geomorph.data.frame(ind = names, coords = Blat$symm.shape, Strain = classifiers$Genotype_Treatment)
MD <- morphol.disparity(coords ~ Strain, groups= ~Strain, iter = 99, data = gdf)

#Procrustes Distance from mean shape
p.dist <- function (a) {
  #### a is an array of ALIGNED Procrusted coordinates (e.g. after GPA). An array is of the form (p x k x n) where p is number of landmarks, k is dimensionality of the data (2D or 3D), n is number of individuals.
  #### Function p.dist calculates and returns procrustes distance(s) between each individual and mean shape of the data.
  p <- dim (a)[3]
  ms <- apply (a, c(1, 2), mean)
  dists <- vector ("numeric", p)
  for(i in 1:p) {
    dists[i] <- sqrt (sum ((a[,,i] - ms)^2)) # this is the formula for PD
  }
  dists
}
## FYI, if necessery you can transform input data to an array, using "arrayspecs" function from the "geomorph" package.For more information look at the package's manual: http://cran.r-project.org/web/packages/geomorph/geomorph.pdf
### Below is a check of the p.dist function
shape <- GPA$coords  # Perform GPA and extract shape variables
dim (shape) # 46 is p, 3 is k, 5 is number of individuals
p.dist (shape) # vector of procrustes distances between each individual and mean shape of the data

#USING MORPHOJ
# 1) Import Centroid Size, symmetric shape, regression score, PCA scores, PCA on regression scores
# 2) Import classifiers
# 3) For PC scores: load into excel and change to "number", otherwise cannot handle shit
# 4) ggplot as per above for boxplots and PCA graphs (change to data.frame first!!!)
# 5) For PC changes, upload PC min and max, and change using as.matrix(), then can plotRefToTarget
# 6) For proc dist: add symm data (through excel), split up coords, create arrays, use above formula to work out distances to mshape

