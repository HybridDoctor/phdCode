#Save untransformed univariate data as UT text file:
UT <- read.delim("C:/Users/Morphometrics lab/Desktop/January 2015/Article 1 folder/Results/UT.txt", quote="") #import into r

View(UT)


d <- data.frame(UT$Strain, UT[8:56]) #create data frame

#Shapiro test:
library("car", lib.loc="C:/Program Files/R/R-3.1.0/library") #open car package

ST <- with(d, aggregate(d[,-1], list(d[,1]), FUN = function(x) 
  shapiro.test(x)$p.value))  #shapiro test with table output
View(ST)

write.table(ST, file="C:/Users/Morphometrics lab/Desktop/January 2015/Article 1 folder/ST.txt", sep="\t")

#Levene test:
LT=rep(NA,49) #creating space to save LT

for (i in 2:50) {
  LT[i] = leveneTest(d[[i]], d$UT.Strain)$"Pr(>F)"[1]
}  #loop Levene test from median for each measurement

View(LT)

#ANOVA:
AT=rep(NA,49) #creating space to save AT
for (i in 2:50) {
  
  AT[i] = summary(aov(d[[i]] ~ d$UT.Strain))[[1]][["Pr(>F)"]][[1]]
  
}  #loop anova for each measurement
View(AT)

#MANOVA
MT <- manova(as.matrix(d[,2:50])~d$UT.Strain, data=d)
summary.manova(MT) #Pillai-Bartlett
summary.manova(MT, test="Wilks")
summary(MT)$stats[1, "Pr(>F)"] #p-value only


#Working on each separately:
  dsplit <- split(d[2:50], d[1]) #split up data-frame

source("http://bit.ly/dasi_inference") # load inference package

#Pairwise Wilcox test:
  WTC.CM=rep(NA,49) #creating space to save WT
for (i in 2:50) {
  WTC.CM[i] = pairwise.wilcox.test(d[[i]], d$UT.Strain, p.adj="none")$"p.value"[1,1]}  #loop Wilcox test for each column getting first p-value
View(WTC.CM)
write.table(WTC.CM, "c:/WTC.CM.txt", sep="\t")