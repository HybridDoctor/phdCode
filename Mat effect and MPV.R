# Import Datasets

rm(list=ls())
setwd("C:/Users/user/Google Drive/Past work/JHE1 data and results/September 2017 hand in/Outputs")

mand <- read.table("C:/Users/user/Google Drive/Past work/JHE1 data and results/September 2017 hand in/mand distances with sex and dam.txt", header = TRUE)
cran <- read.table("C:/Users/user/Google Drive/Past work/JHE1 data and results/September 2017 hand in/cran distances with sex and dam.txt", header = TRUE)

library("plyr", lib.loc="C:/Program Files/R/R-3.4.1/library") #load plyr
library("tidyr", lib.loc="C:/Program Files/R/R-3.4.1/library") #load 
library("boot", lib.loc="C:/Program Files/R/R-3.4.1/library")

### MATERNAL EFFECT

# 1. Find averages of hybrids with different parents (AxB AND BxA)

  Mand <- unite(mand, "Strain_Dam", c("Strain", "Dam"), sep="_", remove = FALSE) # new column
  Cran <- unite(cran, "Strain_Dam", c("Strain", "Dam"), sep="_", remove = FALSE) # new column
  
  Mand <- unite(Mand, "Strain_Sex", c("Strain", "Sex"), sep="_", remove = FALSE) # new column
  Cran <- unite(Cran, "Strain_Sex", c("Strain", "Sex"), sep="_", remove = FALSE)

  Mand.means <- ddply(Mand, c("Strain"), numcolwise(mean))
  Cran.means <- ddply(Cran, c("Strain"), numcolwise(mean))
  
  Mand.mat.means <- ddply(Mand, c("Strain_Dam"), numcolwise(mean)) # data frame of means
  Cran.mat.means <- ddply(Cran, c("Strain_Dam"), numcolwise(mean))
  
  Mand.sex.means <- ddply(Mand, c("Strain_Sex"), numcolwise(mean)) # data frame of means
  Cran.sex.means <- ddply(Cran, c("Strain_Sex"), numcolwise(mean))

# 2. Compare the two hybrids of each set to see if maternal effect

  splitM <- split(Mand, Mand$Strain_Dam) # split up group
  splitC <- split(Cran, Cran$Strain_Dam) # split up group

  # CASTxCZE and CZExCAST
  CACZM <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CACZM[i] <- t.test(splitM$`CASxCZE_CAST-EiJ`[,i], splitM$`CASxCZE_CZECHI-EiJ`[,i])[["p.value"]]
    }
  CACZM
  
  CACZC <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CACZC[i] <- t.test(splitC$`CACZF1_CAST-EiJ`[,i], splitC$`CACZF1_CZECHI-EiJ`[,i])[["p.value"]]
  }
  CACZC

  # CASTxWSB and WSBxCAST
  CAWSM <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CAWSM[i] <- t.test(splitM$`CASxWSB_CAST-EiJ`[,i], splitM$`CASxWSB_WSB`[,i])[["p.value"]]
  }
  CAWSM
  
  CAWSC <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CAWSC[i] <- t.test(splitC$`CAWSF1_CAST-EiJ`[,i], splitC$`CAWSF1_WSB`[,i])[["p.value"]]
  }
  CAWSC

  # CZExWSB and WSBxCZE
  CZWSM <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZWSM[i] <- t.test(splitM$`WSBxCZE_CZECHI-EiJ`[,i], splitM$`WSBxCZE_WSB`[,i])[["p.value"]]
  }
  CZWSM
  
  CZWSC <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZWSC[i] <- t.test(splitC$`WSCZF1_CZECHI-EiJ`[,i], splitC$`WSCZF1_WSB`[,i])[["p.value"]]
  }
  CZWSC

 ## TABLE of hybrid comparisons for maternal effect!!!
  
  Mand_Comparisons <- data.frame(CACZM, CAWSM, CZWSM)
  Cran_Comparisons <- data.frame(CACZC, CAWSC, CZWSC)
  View(Mand_Comparisons)
  View(Cran_Comparisons)

### MPV ANALYSIS

SplitMand <- split(Mand, Mand$Strain) # split up group
SplitCran <- split(Cran, Cran$Strain) # split up group

samplemean <- function(x, d) {
  return(mean(x[d]))
}

MPVl <- function(a, b) {
  x <- boot(data=a, statistic = samplemean, R = 9999)
  m <- boot.ci(x, type="basic")$basic[1,4]
  y <- boot(data=b, statistic = samplemean, R = 9999)
  n <- boot.ci(y, type="basic")$basic[1,4]
  z <- 0.5 * (m + n)
  return(z)
}

MPVu <- function(a, b) {
  x <- boot(data=a, statistic = samplemean, R = 9999)
  o <- boot.ci(x, type="basic")$basic[1,5]
  y <- boot(data=b, statistic = samplemean, R = 9999)
  p <- boot.ci(y, type="basic")$basic[1,5]
  w <- 0.5 * (o + p)
  return(w)
}


#Mand

  #CACZ MPVs: lower and upper
  CASCZEM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CASCZEM_MPVl[i] <- MPVl(SplitMand$CAST[,i], SplitMand$CZECHI[,i])
  }
  CASCZEM_MPVl
  
  CASCZEM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CASCZEM_MPVu[i] <- MPVu(SplitMand$CAST[,i], SplitMand$CZECHI[,i])
  }
  CASCZEM_MPVu

  #CAWS MPVs: lower and upper
  CASWSBM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CASWSBM_MPVl[i] <- MPVl(SplitMand$CAST[,i], SplitMand$WSB[,i])
  }
  CASWSBM_MPVl
  
  CASWSBM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CASWSBM_MPVu[i] <- MPVu(SplitMand$CAST[,i], SplitMand$WSB[,i])
  }
  CASWSBM_MPVu

  #CZWS MPVs: lower and upper
  CZEWSBM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZEWSBM_MPVl[i] <- MPVl(SplitMand$CZECHI[,i], SplitMand$WSB[,i])
  }
  CZEWSBM_MPVl
  
  CZEWSBM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZEWSBM_MPVu[i] <- MPVu(SplitMand$CZECHI[,i], SplitMand$WSB[,i])
  }
  CZEWSBM_MPVu 

# Crans
  #CACZ MPVs: lower and upper  
  CASCZEC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CASCZEC_MPVl[i] <- MPVl(SplitCran$CAST[,i], SplitCran$CAST[,i])
  }
  CASCZEC_MPVl
  
  CASCZEC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CASCZEC_MPVu[i] <- MPVu(SplitCran$CAST[,i], SplitCran$CAST[,i])
  }
  CASCZEC_MPVu
  
  #CAWS MPVs: lower and upper
  CASWSBC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CASWSBC_MPVl[i] <- MPVl(SplitCran$CZECHI[,i], SplitCran$CZECHI[,i])
  }
  CASWSBC_MPVl
  
  CASWSBC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CASWSBC_MPVu[i] <- MPVu(SplitCran$CZECHI[,i], SplitCran$CZECHI[,i])
  }
  CASWSBC_MPVu
  
  #CZWS MPVs: lower and upper
  CZEWSBC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZEWSBC_MPVl[i] <- MPVl(SplitCran$CZECHI[,i], SplitCran$WSB[,i])
  }
  CZEWSBC_MPVl
  
  CZEWSBC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZEWSBC_MPVu[i] <- MPVu(SplitCran$CZECHI[,i], SplitCran$WSB[,i])
  }
  CZEWSBC_MPVu 

  ## Tables of the MPV Cis
  
  MandMPVs <- data.frame(CASCZEM_MPVl, CASCZEM_MPVu, CASWSBM_MPVl, CASWSBM_MPVu, CZEWSBM_MPVl, CZEWSBM_MPVu)
  View(MandMPVs)
  CranMPVs <- data.frame(CASCZEC_MPVl, CASCZEC_MPVu, CASWSBC_MPVl, CASWSBC_MPVu, CZEWSBC_MPVl, CZEWSBC_MPVu)
  View(CranMPVs) 
  
### More maternal effect
# 3. Find MPVs of split and recombined groups

  SplitSexMand <- split(Mand, Mand$Strain_Sex) # split up group
  SplitSexCran <- split(Cran, Cran$Strain_Sex) # split up group  
  
# Mand  
    
  # MPV of CAST females and CZE males
  CACZM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CACZM_MPVl[i] <- MPVl(SplitSexMand$CAST_F[,i], SplitSexMand$CZECHI_M[,i])
  }
  CACZM_MPVl
  
  CACZM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CACZM_MPVu[i] <- MPVu(SplitSexMand$CAST_F[,i], SplitSexMand$CZECHI_M[,i])
  }
  CACZM_MPVu
  
  # MPV of CZE females and CAST males
  CZCAM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZCAM_MPVl[i] <- MPVl(SplitSexMand$CZECHI_F[,i], SplitSexMand$CAST_M[,i])
  }
  CZCAM_MPVl
  
  CZCAM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZCAM_MPVu[i] <- MPVu(SplitSexMand$CZECHI_F[,i], SplitSexMand$CAST_M[,i])
  }
  CZCAM_MPVu  
  
  # MPV of CAST females and WSB males
  CAWSM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CAWSM_MPVl[i] <- MPVl(SplitSexMand$CAST_F[,i], SplitSexMand$WSB_M[,i])
  }
  CAWSM_MPVl
  
  CAWSM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CAWSM_MPVu[i] <- MPVu(SplitSexMand$CAST_F[,i], SplitSexMand$WSB_M[,i])
  }
  CAWSM_MPVu
  
  # MPV of WSB females and CAST males  
  WSCAM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    WSCAM_MPVl[i] <- MPVl(SplitSexMand$WSB_F[,i], SplitSexMand$CAST_M[,i])
  }
  WSCAM_MPVl
  
  WSCAM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    WSCAM_MPVu[i] <- MPVu(SplitSexMand$WSB_F[,i], SplitSexMand$CAST_M[,i])
  }
  WSCAM_MPVu  
  
  # MPV of CZE females and WSB males
  CZWSM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZWSM_MPVl[i] <- MPVl(SplitSexMand$CZECHI_F[,i], SplitSexMand$WSB_M[,i])
  }
  CZWSM_MPVl
  
  CZWSM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    CZWSM_MPVu[i] <- MPVu(SplitSexMand$CZECHI_F[,i], SplitSexMand$WSB_M[,i])
  }
  CZWSM_MPVu  
  
  # MPV of WSB females and CZE males
  WSCZM_MPVl <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    WSCZM_MPVl[i] <- MPVl(SplitSexMand$WSB_F[,i], SplitSexMand$CZECHI_M[,i])
  }
  WSCZM_MPVl
  
  WSCZM_MPVu <- ncol(Mand)
  for (i in 7:ncol(Mand)){
    WSCZM_MPVu[i] <- MPVu(SplitSexMand$WSB_F[,i], SplitSexMand$CZECHI_M[,i])
  }
  WSCZM_MPVu
  
# Cran  
  # MPV of CAST females and CZE males
  CACZC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CACZC_MPVl[i] <- MPVl(SplitSexCran$CAST_F[,i], SplitSexCran$CZECHI_M[,i])
  }
  CACZC_MPVl
  
  CACZC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CACZC_MPVu[i] <- MPVu(SplitSexCran$CAST_F[,i], SplitSexCran$CZECHI_M[,i])
  }
  CACZC_MPVu
  
  # MPV of CZE females and CAST males
  CZCAC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZCAC_MPVl[i] <- MPVl(SplitSexCran$CZECHI_F[,i], SplitSexCran$CAST_M[,i])
  }
  CZCAC_MPVl
  
  CZCAC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZCAC_MPVu[i] <- MPVu(SplitSexCran$CZECHI_F[,i], SplitSexCran$CAST_M[,i])
  }
  CZCAC_MPVu  
  
  # MPV of CAST females and WSB males
  CAWSC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CAWSC_MPVl[i] <- MPVl(SplitSexCran$CAST_F[,i], SplitSexCran$WSB_M[,i])
  }
  CAWSC_MPVl
  
  CAWSC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CAWSC_MPVu[i] <- MPVu(SplitSexCran$CAST_F[,i], SplitSexCran$WSB_M[,i])
  }
  CAWSC_MPVu
  
  # MPV of WSB females and CAST males  
  WSCAC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    WSCAC_MPVl[i] <- MPVl(SplitSexCran$WSB_F[,i], SplitSexCran$CAST_M[,i])
  }
  WSCAC_MPVl
  
  WSCAC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    WSCAC_MPVu[i] <- MPVu(SplitSexCran$WSB_F[,i], SplitSexCran$CAST_M[,i])
  }
  WSCAC_MPVu  
  
  # MPV of CZE females and WSB males
  CZWSC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZWSC_MPVl[i] <- MPVl(SplitSexCran$CZECHI_F[,i], SplitSexCran$WSB_M[,i])
  }
  CZWSC_MPVl
  
  CZWSC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    CZWSC_MPVu[i] <- MPVu(SplitSexCran$CZECHI_F[,i], SplitSexCran$WSB_M[,i])
  }
  CZWSC_MPVu  
  
  # MPV of WSB females and CZE males
  WSCZC_MPVl <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    WSCZC_MPVl[i] <- MPVl(SplitSexCran$WSB_F[,i], SplitSexCran$CZECHI_M[,i])
  }
  WSCZC_MPVl
  
  WSCZC_MPVu <- ncol(Cran)
  for (i in 7:ncol(Cran)){
    WSCZC_MPVu[i] <- MPVu(SplitSexCran$WSB_F[,i], SplitSexCran$CZECHI_M[,i])
  }
  WSCZC_MPVu    

 ## Tables of recombined MPVs!!!
  MandMPVs_recom <- data.frame(CACZM_MPVl, CACZM_MPVu, CZCAM_MPVl,CZCAM_MPVu, WSCZM_MPVl, WSCZM_MPVu, CZWSM_MPVl,CZWSM_MPVu,CAWSM_MPVl, CAWSM_MPVu, WSCAM_MPVl,WSCAM_MPVu)
  View(MandMPVs_recom)
  CranMPVs_recom <- data.frame(CACZC_MPVl, CACZC_MPVu, CZCAC_MPVl,CZCAC_MPVu, WSCZC_MPVl, WSCZC_MPVu, CZWSC_MPVl,CZWSC_MPVu,CAWSC_MPVl, CAWSC_MPVu, WSCAC_MPVl,WSCAC_MPVu)
  View(CranMPVs_recom)
    
# 4. Compare hybrids to split/recombined MPV
#Mand  
  # CASTxCZE
    #Is hybrid smaller than MPVl
    CACZsmallM <- Mand.mat.means[2,2:ncol(Mand.mat.means)] < MandMPVs_recom[7:nrow(MandMPVs_recom),1]
    #Is hybrid larger than MPVu
    CACZlargeM <- Mand.mat.means[2,2:ncol(Mand.mat.means)] > MandMPVs_recom[7:nrow(MandMPVs_recom),2]
  
  # CZExCAST
    #Is hybrid smaller than MPVl
    CZCAsmallM <- Mand.mat.means[3,2:ncol(Mand.mat.means)] < MandMPVs_recom[7:nrow(MandMPVs_recom),3]
    #Is hybrid larger than MPVu
    CZCAlargeM <- Mand.mat.means[3,2:ncol(Mand.mat.means)] > MandMPVs_recom[7:nrow(MandMPVs_recom),4]
    
  # CASTxWSB
    #Is hybrid smaller than MPVl
    CAWSsmallM <- Mand.mat.means[4,2:ncol(Mand.mat.means)] < MandMPVs_recom[7:nrow(MandMPVs_recom),9]
    #Is hybrid larger than MPVu
    CAWSlargeM <- Mand.mat.means[4,2:ncol(Mand.mat.means)] > MandMPVs_recom[7:nrow(MandMPVs_recom),10]
  
  # WSBxCAST
    #Is hybrid smaller than MPVl
    WSCAsmallM <- Mand.mat.means[5,2:ncol(Mand.mat.means)] < MandMPVs_recom[7:nrow(MandMPVs_recom),11]
    #Is hybrid larger than MPVu
    WSCAlargeM <- Mand.mat.means[5,2:ncol(Mand.mat.means)] > MandMPVs_recom[7:nrow(MandMPVs_recom),12]
    
  # CZExWSB
    #Is hybrid smaller than MPVl
    CZWSsmallM <- Mand.mat.means[8,2:ncol(Mand.mat.means)] < MandMPVs_recom[7:nrow(MandMPVs_recom),7]
    #Is hybrid larger than MPVu
    CZWSlargeM <- Mand.mat.means[8,2:ncol(Mand.mat.means)] > MandMPVs_recom[7:nrow(MandMPVs_recom),8]
  
  # WSBxCZE  
    #Is hybrid smaller than MPVl
    WSCZsmallM <- Mand.mat.means[9,2:ncol(Mand.mat.means)] < MandMPVs_recom[7:nrow(MandMPVs_recom),5]
    #Is hybrid larger than MPVu
    WSCZlargeM <- Mand.mat.means[9,2:ncol(Mand.mat.means)] > MandMPVs_recom[7:nrow(MandMPVs_recom),6]
    
#Cran  
    # CASTxCZE
    #Is hybrid smaller than MPVl
    CACZsmallC <- Cran.mat.means[1,2:ncol(Cran.mat.means)] < CranMPVs_recom[7:nrow(CranMPVs_recom),1]
    #Is hybrid larger than MPVu
    CACZlargeC <- Cran.mat.means[1,2:ncol(Cran.mat.means)] > CranMPVs_recom[7:nrow(CranMPVs_recom),2]
    
    # CZExCAST
    #Is hybrid smaller than MPVl
    CZCAsmallC <- Cran.mat.means[2,2:ncol(Cran.mat.means)] < CranMPVs_recom[7:nrow(CranMPVs_recom),3]
    #Is hybrid larger than MPVu
    CZCAlargeC <- Cran.mat.means[2,2:ncol(Cran.mat.means)] > CranMPVs_recom[7:nrow(CranMPVs_recom),4]
    
    # CASTxWSB
    #Is hybrid smaller than MPVl
    CAWSsmallC <- Cran.mat.means[4,2:ncol(Cran.mat.means)] < CranMPVs_recom[7:nrow(CranMPVs_recom),9]
    #Is hybrid larger than MPVu
    CAWSlargeC <- Cran.mat.means[4,2:ncol(Cran.mat.means)] > CranMPVs_recom[7:nrow(CranMPVs_recom),10]
    
    # WSBxCAST
    #Is hybrid smaller than MPVl
    WSCAsmallC <- Cran.mat.means[5,2:ncol(Cran.mat.means)] < CranMPVs_recom[7:nrow(CranMPVs_recom),11]
    #Is hybrid larger than MPVu
    WSCAlargeC <- Cran.mat.means[5,2:ncol(Cran.mat.means)] > CranMPVs_recom[7:nrow(CranMPVs_recom),12]
    
    # CZExWSB
    #Is hybrid smaller than MPVl
    CZWSsmallC <- Cran.mat.means[8,2:ncol(Cran.mat.means)] < CranMPVs_recom[7:nrow(CranMPVs_recom),7]
    #Is hybrid larger than MPVu
    CZWSlargeC <- Cran.mat.means[8,2:ncol(Cran.mat.means)] > CranMPVs_recom[7:nrow(CranMPVs_recom),8]
    
    # WSBxCZE  
    #Is hybrid smaller than MPVl
    WSCZsmallC <- Cran.mat.means[9,2:ncol(Cran.mat.means)] < CranMPVs_recom[7:nrow(CranMPVs_recom),5]
    #Is hybrid larger than MPVu
    WSCZlargeC <- Cran.mat.means[9,2:ncol(Cran.mat.means)] > CranMPVs_recom[7:nrow(CranMPVs_recom),6]
    
 ##Tables of where hybrids < MPVl and > MPVu   
    Hybrid_fitM <- data.frame(CAWSsmallM[1,], CAWSlargeM[1,], WSCAsmallM[1,], WSCAlargeM[1,], CACZsmallM[1,], CACZlargeM[1,], CZCAsmallM[1,], CZCAlargeM[1,], CZWSsmallM[1,], CZWSlargeM[1,], WSCZsmallM[1,],WSCZlargeM[1,])
    Hybrid_fitC <- data.frame(CAWSsmallC[1,], CAWSlargeC[1,], WSCAsmallC[1,], WSCAlargeC[1,], CACZsmallC[1,], CACZlargeC[1,], CZCAsmallC[1,], CZCAlargeC[1,], CZWSsmallC[1,], CZWSlargeC[1,], WSCZsmallC[1,],WSCZlargeC[1,])
  
### Table of hybrid results/proportions

    
    
    
    
### Comparisons
  # Count if Mand_Comparisons and Cran_comparisons are <0.05
  sum(Mand_Comparisons$CACZM < 0.05, na.rm = TRUE)
  sum(Mand_Comparisons$CAWSM < 0.05, na.rm = TRUE)
  sum(Mand_Comparisons$CZWSM < 0.05, na.rm = TRUE)
  sum(Cran_Comparisons$CACZC < 0.05, na.rm = TRUE)
  sum(Cran_Comparisons$CAWSC < 0.05, na.rm = TRUE)
  sum(Cran_Comparisons$CZWSC < 0.05, na.rm = TRUE)
 
  # Hybrids intermediate: logic tests
  CASXCZEM_inter <- ((Mand.means[2,]>Mand.means[1,])&(Mand.means[2,]<Mand.means[4,])) | ((Mand.means[2,]<Mand.means[1,])&(Mand.means[2,]>Mand.means[4,]))
  CASxWSBM_inter <- ((Mand.means[3,]>Mand.means[1,])&(Mand.means[3,]<Mand.means[5,])) | ((Mand.means[3,]<Mand.means[1,])&(Mand.means[3,]>Mand.means[5,]))
  CZE_WSBM_inter <- ((Mand.means[6,]>Mand.means[4,])&(Mand.means[6,]<Mand.means[5,])) | ((Mand.means[6,]<Mand.means[4,])&(Mand.means[6,]>Mand.means[5,]))
  CASXCZEC_inter <- ((Cran.means[1,]>Cran.means[2,])&(Cran.means[1,]<Cran.means[4,])) | ((Cran.means[1,]<Cran.means[2,])&(Cran.means[1,]>Cran.means[4,]))
  CASxWSBC_inter <- ((Cran.means[3,]>Cran.means[2,])&(Cran.means[3,]<Cran.means[5,])) | ((Cran.means[3,]<Cran.means[2,])&(Cran.means[3,]>Cran.means[5,]))
  CZE_WSBC_inter <- ((Cran.means[6,]>Cran.means[4,])&(Cran.means[6,]<Cran.means[5,])) | ((Cran.means[6,]<Cran.means[4,])&(Cran.means[6,]>Cran.means[5,]))
  
  Hybrid_inter_Mand <- data.frame(CASXCZEM_inter[1,], CASxWSBM_inter[1,], CZE_WSBM_inter[1,])
  Hybrid_inter_Cran <- data.frame(CASXCZEC_inter[1,], CASxWSBC_inter[1,], CZE_WSBC_inter[1,])
  
  # Larger than both - non significant
  CASXCZEM_larger <- (Mand.means[2,]>Mand.means[1,])&(Mand.means[2,]>Mand.means[4,]) 
  CASxWSBM_larger <- (Mand.means[3,]>Mand.means[1,])&(Mand.means[3,]>Mand.means[5,])
  CZE_WSBM_larger <- (Mand.means[6,]>Mand.means[4,])&(Mand.means[6,]>Mand.means[5,])
  CASXCZEC_larger <- (Cran.means[1,]>Cran.means[2,])&(Cran.means[1,]>Cran.means[4,])
  CASxWSBC_larger <- (Cran.means[3,]>Cran.means[2,])&(Cran.means[3,]>Cran.means[5,])
  CZE_WSBC_larger <- (Cran.means[6,]>Cran.means[4,])&(Cran.means[6,]>Cran.means[5,])
  
  Hybrid_larger_Mand <- data.frame(CASXCZEM_larger[1,], CASxWSBM_larger[1,], CZE_WSBM_larger[1,])
  Hybrid_larger_Cran <- data.frame(CASXCZEC_larger[1,], CASxWSBC_larger[1,], CZE_WSBC_larger[1,])
  
  #Smaller than both - non-significant
  CASXCZEM_smaller <- (Mand.means[2,]<Mand.means[1,])&(Mand.means[2,]<Mand.means[4,]) 
  CASxWSBM_smaller <- (Mand.means[3,]<Mand.means[1,])&(Mand.means[3,]<Mand.means[5,])
  CZE_WSBM_smaller <- (Mand.means[6,]<Mand.means[4,])&(Mand.means[6,]<Mand.means[5,])
  CASXCZEC_smaller <- (Cran.means[1,]<Cran.means[2,])&(Cran.means[1,]<Cran.means[4,])
  CASxWSBC_smaller <- (Cran.means[3,]<Cran.means[2,])&(Cran.means[3,]<Cran.means[5,])
  CZE_WSBC_smaller <- (Cran.means[6,]<Cran.means[4,])&(Cran.means[6,]<Cran.means[5,])
  
  Hybrid_smaller_Mand <- data.frame(CASXCZEM_smaller[1,], CASxWSBM_smaller[1,], CZE_WSBM_smaller[1,])
  Hybrid_smaller_Cran <- data.frame(CASXCZEC_smaller[1,], CASxWSBC_smaller[1,], CZE_WSBC_smaller[1,])
  
  # Heterosis 
  CACZ_het_Mand <- Mand.means[2,2:22]>MandMPVs[7:27,2]
  CAWS_het_Mand <- Mand.means[3,2:22]>MandMPVs[7:27,4]
  CZWS_het_Mand <- Mand.means[6,2:22]>MandMPVs[7:27,6]
  CACZ_het_Cran <- Cran.means[1,2:39]>CranMPVs[7:44,2]
  CAWS_het_Cran <- Cran.means[3,2:39]>CranMPVs[7:44,4]
  CZWS_het_Cran <- Cran.means[6,2:39]>CranMPVs[7:44,6]
  
  Hybrid_het_Mand <- data.frame(CACZ_het_Mand[1,], CAWS_het_Mand[1,], CZWS_het_Mand[1,])
  Hybrid_het_Cran <- data.frame(CACZ_het_Cran[1,], CAWS_het_Cran[1,], CZWS_het_Cran[1,])
  
  # Dysgenesis 
  CACZ_dys_Mand <- Mand.means[2,2:22]<MandMPVs[7:27,1]
  CAWS_dys_Mand <- Mand.means[3,2:22]<MandMPVs[7:27,3]
  CZWS_dys_Mand <- Mand.means[6,2:22]<MandMPVs[7:27,5]
  CACZ_dys_Cran <- Cran.means[1,2:39]<CranMPVs[7:44,1]
  CAWS_dys_Cran <- Cran.means[3,2:39]<CranMPVs[7:44,3]
  CZWS_dys_Cran <- Cran.means[6,2:39]<CranMPVs[7:44,5]
  
  Hybrid_dys_Mand <- data.frame(CACZ_dys_Mand[1,], CAWS_dys_Mand[1,], CZWS_dys_Mand[1,])
  Hybrid_dys_Cran <- data.frame(CACZ_dys_Cran[1,], CAWS_dys_Cran[1,], CZWS_dys_Cran[1,])
  

####OUTPUTS
  
  write.csv(Mand.means, "MandMeansPerStrain.csv") 
  write.csv(Cran.means, "CranMeansPerStrain.csv")
  
  write.csv(Mand.mat.means,"MandMeansPerDam.csv") 
  write.csv(Cran.mat.means,"CranMeansPerDam.csv")
  
  write.csv(Mand.sex.means,"MandMeansPerSex.csv")
  write.csv(Cran.sex.means,"CranMeansPerSex.csv")
  
  write.csv(Mand_Comparisons,"MandDamTests.csv") 
  write.csv(Cran_Comparisons,"CranDamTests.csv")
  
  write.csv( MandMPVs,"MandMPVs.csv")
  write.csv(CranMPVs,"CranMPVs.csv")
  
  write.csv(MandMPVs_recom,"MandMPVsrecom.csv") 
  write.csv(CranMPVs_recom,"CranMPVsrecom.csv")

  write.csv(Hybrid_fitM,"MandMeansPerStrain.csv")
  write.csv(Hybrid_fitC,"MandMeansPerStrain.csv")
 
  write.csv(Hybrid_het_Cran, "hybridHetCran.csv")
  write.csv(Hybrid_het_Mand, "hybridHetMand.csv")
  
  
# Quick stats
  # Intermediate
  sum(Hybrid_inter_Mand[2:22,1]=="TRUE", na.rm = TRUE) #exlude GEOMEAN
  sum(Hybrid_inter_Mand[2:22,2]=="TRUE", na.rm = TRUE)
  sum(Hybrid_inter_Mand[2:22,3]=="TRUE", na.rm = TRUE)
  sum(Hybrid_inter_Cran[2:39,1]=="TRUE", na.rm = TRUE)
  sum(Hybrid_inter_Cran[2:39,2]=="TRUE", na.rm = TRUE)
  sum(Hybrid_inter_Cran[2:39,3]=="TRUE", na.rm = TRUE)
  
  # Larger, non-significant
  sum(Hybrid_larger_Mand[2:22,1]=="TRUE", na.rm = TRUE)
  sum(Hybrid_larger_Mand[2:22,2]=="TRUE", na.rm = TRUE)
  sum(Hybrid_larger_Mand[2:22,3]=="TRUE", na.rm = TRUE)
  sum(Hybrid_larger_Cran[2:39,1]=="TRUE", na.rm = TRUE)
  sum(Hybrid_larger_Cran[2:39,2]=="TRUE", na.rm = TRUE)
  sum(Hybrid_larger_Cran[2:39,3]=="TRUE", na.rm = TRUE)
  
  # Smaller, non-significant
  sum(Hybrid_smaller_Mand[2:22,1]=="TRUE", na.rm = TRUE)
  sum(Hybrid_smaller_Mand[2:22,2]=="TRUE", na.rm = TRUE)
  sum(Hybrid_smaller_Mand[2:22,3]=="TRUE", na.rm = TRUE)
  sum(Hybrid_smaller_Cran[2:39,1]=="TRUE", na.rm = TRUE)
  sum(Hybrid_smaller_Cran[2:39,2]=="TRUE", na.rm = TRUE)
  sum(Hybrid_smaller_Cran[2:39,3]=="TRUE", na.rm = TRUE)  
  
  #Heterosis
  sum(Hybrid_het_Mand$CACZ_het_Mand.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_het_Mand$CAWS_het_Mand.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_het_Mand$CZWS_het_Mand.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_het_Cran$CACZ_het_Cran.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_het_Cran$CAWS_het_Cran.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_het_Cran$CZWS_het_Cran.1...== "TRUE", na.rm = TRUE)
  
  #Dysgenesis
  sum(Hybrid_dys_Mand$CACZ_dys_Mand.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_dys_Mand$CAWS_dys_Mand.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_dys_Mand$CZWS_dys_Mand.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_dys_Cran$CACZ_dys_Cran.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_dys_Cran$CAWS_dys_Cran.1...== "TRUE", na.rm = TRUE)
  sum(Hybrid_dys_Cran$CZWS_dys_Cran.1...== "TRUE", na.rm = TRUE)
  
  # Mat-effect: hybrid-hybrid comparisons
  sum(Mand_Comparisons$CACZM[7:27]<0.05, na.rm = TRUE)
  sum(Mand_Comparisons$CAWSM[7:27]<0.05, na.rm = TRUE)
  sum(Mand_Comparisons$CZWSM[7:27]<0.05, na.rm = TRUE)
  sum(Cran_Comparisons$CACZC[7:44]<0.05, na.rm = TRUE)
  sum(Cran_Comparisons$CAWSC[7:44]<0.05, na.rm = TRUE)
  sum(Cran_Comparisons$CZWSC[7:44]<0.05, na.rm = TRUE)
  
  # Mat effect: hybrid to MPV comparison
    # Heterosis
  sum(Hybrid_fitM$CACZlargeM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$CZCAlargeM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$CAWSlargeM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$WSCAlargeM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$CZWSlargeM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$WSCZlargeM.1...[1:21]=="TRUE", na.rm=TRUE)
  
  sum(Hybrid_fitC$CACZlargeC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$CZCAlargeC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$CAWSlargeC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$WSCAlargeC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$CZWSlargeC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$WSCZlargeC.1...[1:38]=="TRUE", na.rm=TRUE)
  
    # Dysgenesis
  sum(Hybrid_fitM$CACZsmallM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$CZCAsmallM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$CAWSsmallM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$WSCAsmallM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$CZWSsmallM.1...[1:21]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitM$WSCZsmallM.1...[1:21]=="TRUE", na.rm=TRUE)
  
  sum(Hybrid_fitC$CACZsmallC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$CZCAsmallC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$CAWSsmallC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$WSCAsmallC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$CZWSsmallC.1...[1:38]=="TRUE", na.rm=TRUE)
  sum(Hybrid_fitC$WSCZsmallC.1...[1:38]=="TRUE", na.rm=TRUE)
  
  
  
  
  