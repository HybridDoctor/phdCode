setwd("C:/Users/Morphometrics lab/Google Drive/Methodology and Results/Section 2/WD") #Working directory for paper 2
# Can remove items from WD using rm()

#Visual:
Colour <- c("WSB" = "royalblue", 
            "CastxWSB" = "firebrick3", 
            "CastXWSB*Cast" = "plum", 
            "CastxWSB*F2" = "salmon1", 
            "CAST-EiJ" = "deepskyblue2",
            "CZECHI-EiJ" = "seagreen3",
            "CastxCzech" = "firebrick3",
            "CastxCzech*Czech" = "plum",
            "CastxCzech*F2" = "salmon1",
            "WSBxCzech" = "firebrick3",
            "WSBxCzech*WSB" = "plum")

#cranial colour
c("WSB" = "royalblue", 
  "CASxWSB" = "firebrick3", 
  "CASxWSB_CAS" = "plum", 
  "CASxWSB_F2" = "salmon1", 
  "CAST" = "deepskyblue2",
  "CZECHI" = "seagreen3",
  "CASxCZE" = "firebrick3",
  "CASxCZE_CAS" = "plum",
  "CASxCZE_CZE" = "plum",
  "CASxCZE_F2" = "salmon1",
  "WSBxCZE" = "firebrick3",
  "PANCEVO" = "violetred",
  "SPRET" = "deeppink",
  "WSBxSPR" = "red")

# to change names use ctl + f

Mandibular.Data <- read.delim("C:/Users/Morphometrics lab/Google Drive/Methodology and Results/Section 2/WD/Mandibular Data.txt", quote="")
 View(Mandibular.Data)
 
 levels(Mandibular.Data$Genotype_Treatment) #make sure levels are correct
 
 # Load Packages necessary
 library(dplyr)
 library(geomorph)
 library(tidyr)
 library(ggplot2)
 
Mand <- tbl_df(Mandibular.Data) # To work easily with dplyr

# Individual strains using dplyr
CastM <- filter(Mand, Genotype_Treatment == "Cast-EiJ")

# Combined strains into hybrid groups using dplyr

#GPA
names <- Mandibular.Data[,1]
coords <- arrayspecs(Mandibular.Data[,14:ncol(Mandibular.Data)], 17, 3) #creates 3D array, p=no. landmarks, k=no. dimensions
dimnames(coords)[[3]] <- names
classifiers <- Mandibular.Data[,1:13] #change latter column- check column numbers!!!
dim(coords) # check to see that there are 3 numbers
GPA <- gpagen(coords)

gdf <- geomorph.data.frame(coords = coords, 
                           ind = names, 
                           side = classifiers$Side, 
                           replicate = classifiers$Replication)
Mand.sym <- bilat.symmetry(A = coords, 
                           ind = ind, 
                           side = side, 
                           replicate = replicate, 
                           object.sym = FALSE, 
                           data = gdf)

Mand.sym$shape.anova



    
#PCA mesh
    library("geomorph", lib.loc="C:/Program Files/R/R-3.1.0/library") #load geomorph
    # Import procrustes aligned averaged data from MorphoJ
    CAWS <- as.data.frame(`CZExWSB,.averaged`)
    coords <- arrayspecs(CAWS[,17:ncol(CAWS)], 17, 3)
    names <- CAWS[,1]
    classifiers <- CAWS[,1:16] #change latter column
    dim(coords)
    PCA <-plotTangentSpace(coords, verbose = TRUE) #check if comparable
    
    # Import links for mandible outline
    
    PC1mesh <- plotRefToTarget(PCA$pc.shapes$PC1min, #reference
                               PCA$pc.shapes$PC1max, #target 
                               method = "points", 
                               links = links, 
                               gridPars=gridPar(link.col = "colourMin", #reference colour
                                                tar.link.col="colourMax", #target colour
                                                tar.link.lwd=2))
    
    PC1mesh <- plotRefToTarget(PCA$pc.shapes$PC1min, #reference
                               PCA$pc.shapes$PC1max, #target 
                               method = "TPS", 
                               links = links)
    
    #4) Colours for all plots: CAST(deepskyblue2), 
    #WSB (royalblue), 
    #CZECHI(seagreen3), 
    #F1 (firebrick3), 
    #B1 (plum), 
    #F2 (salmon1)   
    
#Allometric regression: groups and together
    library("geomorph", lib.loc="C:/Program Files/R/R-3.1.0/library") #load geomorph
    
    # Import procrustes aligned averaged data from MorphoJ
    CAWS <- as.data.frame(`CZExWSB,.averaged`)
    coords <- arrayspecs(CAWS[,17:ncol(CAWS)], 17, 3)
    names <- CAWS[,1]
    classifiers <- CAWS[,1:16] #change latter column
    dim(coords)
    
    gdfCAWS <- geomorph.data.frame(coords = coords, 
                               ind = names, 
                               Cent = classifiers$Centroid.Size, 
                               LogCS = classifiers$Log.Centroid.Size,
                               Strain = classifiers$Genotype_Treatment)
    
    AllCAWS <- procD.allometry(coords ~ Cent, # Effect of size, species and size within species
                           f2 = NULL, # or ~Strain,
                           f3 = NULL, 
                           logsz = TRUE, 
                           data = gdfCAWS, 
                           iter = 9999,
                           RRPP = TRUE)
    summary(AllCAWS)
    
    # Significance of each group
    
    library("dplyr", lib.loc="C:/Program Files/R/R-3.3.0/library")
    
    CastCACZ <- filter(CACZ, Genotype_Treatment == "CAST-EiJ")
    CzeCACZ <- filter(CACZ, Genotype_Treatment == "CZECHI-EiJ")
    
    ### ???
    
#Morphological Disparity without size
    
    library("geomorph", lib.loc="C:/Program Files/R/R-3.1.0/library") #load geomorph
    
    gdf <- geomorph.data.frame(ind = names, 
                               coords = coords, 
                               Cent = classifiers$Centroid.Size,
                               Strain = classifiers$Genotype_Treatment)
    
    MD <- morphol.disparity(coords ~ 1, # by Strain instead?
                            groups= ~Strain, 
                            iter = 999, 
                            data = gdf)
    summary(MD)
    
#Morphological disparity accounting for allometry
    #As above
    
    MD <- morphol.disparity(coords ~ Cent, #use Cent + Strain instead?
                            groups= ~Strain, 
                            iter = 999, 
                            data = gdf)
    summary(MD)
    
    
# TPS    
  # Import 3D data for PC1 and PC2
    
    library("geomorph", lib.loc="C:/Program Files/R/R-3.1.0/library") #load geomorph
    
    #IMPORTED DATA FROM MORPHO J
    
    # Import centroid size
    
    
    # Import PCA data- name CZWSPCA, CZWSPCA and CZWSPCA, then use as.data.frame
    
    library("ggplot2", lib.loc="C:/Program Files/R/R-3.3.0/library")
    
    CZWSPCA <- as.data.frame(CZWSPCA)
    PCplot <- ggplot(data = CZWSPCA, aes(x= PC1, y = PC2, colour = Genotype_Treatment)) + 
      geom_point(shape=20, size=3)
    
    #Make hulls
    find_hull <- function(CZWSPCA) CZWSPCA[chull(CZWSPCA$PC1, CZWSPCA$PC2),]
    library("plyr", lib.loc="C:/Program Files/R/R-3.1.0/library") #need plyr
    hulls <- ddply(CZWSPCA, "Genotype_Treatment", find_hull)
    PCplot + geom_polygon(data=hulls, alpha=0.3, aes(fill=Genotype_Treatment))
    #Stylizing PCA
    PCP <- PCplot + 
      geom_polygon(data=hulls, alpha=0.3, size = 0.5, aes(fill=Genotype_Treatment)) + 
      #or use for CVA: PCP <- PCplot + stat_ellipse(geom = "polygon", alpha=0.3, size = 0.5, aes(fill=Genotype_Treatment)) +
      theme_minimal() +
      theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = 'none',
            panel.background = element_rect(color = "grey40", size = 0.5),
            axis.text = element_text(face = "italic", size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 14, face = "bold", vjust = 1)
      )
    PCname <- PCP +
      labs(x = "PC1 (perc)", 
           y = "PC2 (perc)", 
           title = "PCA plot of CAST (light blue) vs WSB (dark blue)") +
      scale_color_manual(values = c("WSB" = "royalblue", 
                                    "CastxWSB" = "firebrick3", 
                                    "CastXWSB*Cast" = "plum", 
                                    "CastxWSB*F2" = "salmon1", 
                                    "CAST-EiJ" = "deepskyblue2",
                                    "CZECHI-EiJ" = "seagreen3",
                                    "CastxCzech" = "firebrick3",
                                    "CastxCzech*Czech" = "plum",
                                    "CastxCzech*F2" = "salmon1",
                                    "WSBxCzech" = "firebrick3",
                                    "WSBxCzech*Czech" = "plum")) +
      scale_fill_manual(values = c("WSB" = "royalblue", 
                                   "CastxWSB" = "firebrick3", 
                                   "CastXWSB*Cast" = "plum", 
                                   "CastxWSB*F2" = "salmon1", 
                                   "CAST-EiJ" = "deepskyblue2",
                                   "CZECHI-EiJ" = "seagreen3",
                                   "CastxCzech" = "firebrick3",
                                   "CastxCzech*Czech" = "plum",
                                   "CastxCzech*F2" = "salmon1",
                                   "WSBxCzech" = "firebrick3",
                                   "WSBxCzech*Czech" = "plum"))
    PCname
    ggsave(filename = "PCname.tiff", plot = PCname, dpi = 700)   
    
    #Things to do:
    #1) x/y labels: PC and variance (eg. "PC1(31%)")
    #2) Take out legend in all but one (or put legend below or above???)
    #3) Points should be larger, and geom lines slightly thicker
    #4) Colours for all plots: CAST(deepskyblue2), 
    #WSB (royalblue), 
    #CZECHI(seagreen3), 
    #F1 (firebrick3), 
    #B1 (plum), 
    #F2 (salmon1) 
    #5) thicker label lines? lines thicker at 0.0    

# Regression analysis
    
    #Split up groups
    
    
    # work out lm for each group
    
    
    # use ANOVAs to compare the 4 or 5 lms
    
    
  # Plots
    # Combine reg score with centroid size datasets
    
    PCplot <- ggplot(data = CACZ, aes(x= Centroid.Size, y = regScore, colour = Genotype_Treatment)) + 
      geom_point(shape=20, size=3) 
    
    PCP <- PCplot + 
      theme_minimal() +
      theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = 'none',
            panel.background = element_rect(color = "grey40", size = 0.5),
            axis.text = element_text(face = "italic"),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 14, face = "bold", vjust = 1)
      )
    
    PCname <- PCP +
      labs(x = "Centroid Size", 
           y = "Regression Score", 
           title = "Regression plot of CAST (blue) and CZECHI (green) and hybrids") +
      scale_color_manual(values = c("WSB" = "royalblue", 
                                    "CastxWSB" = "firebrick3", 
                                    "CastXWSB*Cast" = "plum", 
                                    "CastxWSB*F2" = "salmon1", 
                                    "CAST-EiJ" = "deepskyblue2",
                                    "CZECHI-EiJ" = "seagreen3",
                                    "CastxCzech" = "firebrick3",
                                    "CastxCzech*Czech" = "plum",
                                    "CastxCzech*F2" = "salmon1",
                                    "WSBxCzech" = "firebrick3",
                                    "WSBxCzech*WSB" = "plum")) +
      scale_fill_manual(values = c("WSB" = "royalblue", 
                                   "CastxWSB" = "firebrick3", 
                                   "CastXWSB*Cast" = "plum", 
                                   "CastxWSB*F2" = "salmon1", 
                                   "CAST-EiJ" = "deepskyblue2",
                                   "CZECHI-EiJ" = "seagreen3",
                                   "CastxCzech" = "firebrick3",
                                   "CastxCzech*Czech" = "plum",
                                   "CastxCzech*F2" = "salmon1",
                                   "WSBxCzech" = "firebrick3",
                                   "WSBxCzech*WSB" = "plum")) +
      geom_smooth(method = "lm", 
                  se = FALSE, 
                  fullrange = TRUE, 
                  size = 1)
    PCname
    
    ggsave(filename = "PCname.tiff", plot = PCname, dpi = 700) 
    
    # Mesh: import start and end matrices
    
    CAWSstart <- as.matrix(CAWSstart) # change to matrix
    CAWSend <- as.matrix(CAWSend)
    plotRefToTarget(CAWSstart, CAWSend, links = links, mag = 0.4)
    
    
# CENTROID SIZE
    
    # Anova
    #Import averaged Cent size of individuals for each group
    CACZAn <- aov(Centroid.Size ~ Genotype_Treatment, data = CACZCent)
    summary(CACZAn)
    
    # Tukey tests
    TukeyCACZ <- TukeyHSD(CACZAn)
    TukeyCACZ
    
    # Levene's test
    LevCACZ <- leveneTest(CACZCent$Centroid.Size, 
                          group = CACZCent$Genotype_Treatment)
    LevCACZ
    
    # Descriptive stats
    library(psych)
    describeBy(CACZCent$Centroid.Size, CACZCent$Genotype_Treatment)
    
    
 # Box plots
    BP <- ggplot(CACZCent, 
                 aes(x= Genotype_Treatment, 
                     y = Centroid.Size, 
                     fill = Genotype_Treatment)) + 
          geom_boxplot() + 
          theme_minimal() 
    BPname <- BP +
      labs(y = "Centroid Size", 
           x = "") +
      scale_fill_manual(values = c("CAST-EiJ" = "deepskyblue2",
                                   "CZECHI-EiJ" = "seagreen3",
                                   "CastxCzech" = "firebrick3",
                                   "CastxCzech*Czech" = "plum",
                                   "CastxCzech*F2" = "salmon1")) + 
      scale_x_discrete(limits= c("CAST-EiJ", 
                                 "WSB", 
                                 "CZECHI-EiJ", 
                                 "CastxCzech", 
                                 "CastxWSB", 
                                 "WSBxCzech", 
                                 "WSBxCzech*WSB", 
                                 "CastXWSB*Cast", 
                                 "CastxCzech*Czech", 
                                 "CastxCzech*F2", 
                                 "CastxWSB*F2"),
                    labels = c("CAST", 
                               "CZECHI", 
                               "F1", 
                               "B1", 
                               "F2"))  +
      theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = 'none',
            panel.background = element_rect(color = "grey40", size = 0.5),
            axis.text = element_text(face = "italic", size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 14, face = "bold", vjust = 2)
      ) + 
      stat_summary(fun.y = mean, 
                   geom = "point", 
                   shape = 18, 
                   size = 6)
   
     BPname
    

    
# Procrustes variances
   
    # import dataset, 
    # split up into groups: 
     split <- split(as.data.frame(Whole_group_Crania), 
                    Whole_group_Crania$Genotype_Treatment, 
                    drop=FALSE)
     
    # convert to coords and array
     coordsCAST <- arrayspecs(split$CAST[,5:ncol(split$CAST)], 57, 3)
     coordsCZECHI
    
    # calculate procrustes distances to mean of each group and F1:
    # work out procrustes variances for each group (use geomorph: morphol.disparity):
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

     dim (coordsCAST) # 46 is p, 3 is k, 5 is number of individuals
     CASTdist <- p.dist (coordsCAST) # vector of procrustes distances between each individual and mean shape of the data
     
    # repeat on all groups
     coordsX <- arrayspecs(split$X[,5:ncol(split$X)], 57, 3)
     dim(coordsX)
     Xdist <- p.dist (coordsX)
     Xmean <- mean(Xdist)
     
     
     
    # bargraphs (x3- 1 for each group: with error bars?)
    
      #make two vectors and coerce into data.frame
    CAWSP <- data.frame(name = c("CAST", "WSB", "F1", "B1", "F2"), 
                        PV = c(0.00082, 0.00096, 0.00092, 0.00139, 0.00168)) # mean in each group
    
    CAWSbar <- ggplot(CAWSP, 
                      aes(x=name, 
                          y=PV, 
                          fill = name)) +
      geom_bar(colour = "grey40", 
               stat = "identity") +
      labs(y = "", 
           x = "") + 
      theme_minimal() +
      scale_fill_manual(values = c("CAST" = "deepskyblue2",
                                   "WSB" = "royalblue",
                                   "F1" = "firebrick3",
                                   "B1" = "plum",
                                   "F2" = "salmon1")) + 
      scale_x_discrete(limits= c("CAST",
                                 "WSB", 
                                 "F1", 
                                 "B1", 
                                 "F2")) +
      theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = 'none',
            panel.background = element_rect(color = "grey40", size = 0.5),
            axis.text = element_text(face = "italic", size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 14, face = "bold", vjust = 2)
      )

    
# Classification
    
    
    # convert to coords and array
    coordsCAST <- arrayspecs(split$CAST[,5:ncol(split$CAST)], 57, 3) #array deats, continue for all groups
    
    msX <- apply (coordsX, c(1, 2), mean) #work out mean shape of P or F1
    
    p <- dim (coordsY)[3] # save dimensionality of p if unknown
    dists <- vector ("numeric", p)
    for(i in 1:p) {
      dists[i] <- sqrt (sum ((coordsY[,,i] - msX)^2)) # this is the formula for PD to a mean
    }
    dists
    
    
    # Create data frames
    CACZClass <- data.frame(Strain = c("CAST", "CAST", "CAST", "CZECHI", "CZECHI", "CZECHI", "F1", "F1", "F1", "F2", "F2", "F2", "B1", "B1", "B1"), 
                            Analysis = c("CAST", "CZECHI", "F1", "CAST", "CZECHI", "F1", "CAST", "CZECHI", "F1", "CAST", "CZECHI", "F1", "CAST", "CZECHI", "F1"), 
                            Centroid = c(28, 2, 0, 2, 26, 2, 0, 7, 23, 2, 23, 5, 1, 5, 0), Procrustes = c(30, 0, 0, 0, 30, 0, 0, 0, 30, 1, 0, 29, 0, 4, 2))
    
    #Create bar plot (stacked)
    
    CACZ_SBP <- ggplot(CACZClass, 
                       aes(x=Strain, y = Centroid, fill = Analysis)) +
      geom_bar(colour = "grey40", 
               stat = "identity") +
      labs(y = "", 
           x = "") + 
      theme_minimal() +
      scale_fill_manual(values = c("CAST" = "deepskyblue2",
                                   "CZECHI" = "seagreen3",
                                   "F1" = "firebrick3",
                                   "B1" = "plum",
                                   "F2" = "salmon1")) + 
      scale_x_discrete(limits= c("CAST",
                                 "CZECHI", 
                                 "F1", 
                                 "B1", 
                                 "F2")) +
      theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = 'none',
            panel.background = element_rect(color = "grey40", size = 0.5),
            axis.text = element_text(face = "italic", size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 14, face = "bold", vjust = 2)
      )
    
#Mixed models
    
    library("dplyr", lib.loc="C:/Program Files/R/R-3.3.0/library")
    
    #split up tables into parents and F1 hybrids
    CWP <- filter(CACZ, Genotype_Treatment %in% c("Cast-EiJ", "CZECHI-EiJ"))
    CWF1 <- filter(CACZ, Genotype_Treatment == "CastxCzech")
    
    CW.mean=rep(NA,999) # create table for distribution of means
    CW.variance=rep(NA,999) # create table for distribution of variances
    
    for (i in 1:999) {
      sample.parent<- sample(CWP$Centroid.Size, size,  replace=TRUE) #resampling parents
      sample.F1 <- sample(CWF1$Centroid.Size, size, replace=TRUE) #resampling F1 hybrids- add 0 to 30/50
      CW <- c(sample.parent, sample.F1) # combined numbers
      CW.mean[i] <- mean(CW) # mean of combined numbers
      CW.variance[i] <- var(CW) #variance of combined numbers
    }
    CW1_var_low <- quantile(CW.variance,0.025)
    CW1_var_high <- quantile(CW.variance, 0.975)
    CW1_var <- mean(CW.variance)
    CW1_mean_low <- quantile(CW.mean, 0.025)
    CW1_mean_high <- quantile(CW.mean, 0.975)
    CW1_mean <- mean(CW.mean)
    
    # Creating vectors and tables
    CW1_vector <- c(CW1_var_low, CW1_var_high, CW1_var, CW1_mean_low, CW1_mean_high, CW1_mean)
    CW2_vector <- c(CW2_var_low, CW2_var_high, CW2_var, CW2_mean_low, CW2_mean_high, CW2_mean)
    CW3_vector <- c(CW3_var_low, CW3_var_high, CW3_var, CW3_mean_low, CW3_mean_high, CW3_mean)
    CW4_vector <- c(CW4_var_low, CW4_var_high, CW4_var, CW4_mean_low, CW4_mean_high, CW4_mean)
    CW5_vector <- c(CW5_var_low, CW5_var_high, CW5_var, CW5_mean_low, CW5_mean_high, CW5_mean)
    CW6_vector <- c(CW6_var_low, CW6_var_high, CW6_var, CW6_mean_low, CW6_mean_high, CW6_mean)
    CW7_vector <- c(CW7_var_low, CW7_var_high, CW7_var, CW7_mean_low, CW7_mean_high, CW7_mean)
    CW8_vector <- c(CW8_var_low, CW8_var_high, CW8_var, CW8_mean_low, CW8_mean_high, CW8_mean)
    CW9_vector <- c(CW9_var_low, CW9_var_high, CW9_var, CW9_mean_low, CW9_mean_high, CW9_mean)
    CW10_vector <- c(CW10_var_low, CW10_var_high, CW10_var, CW10_mean_low, CW10_mean_high, CW10_mean)
    CW11_vector <- c(CW11_var_low, CW11_var_high, CW11_var, CW11_mean_low, CW11_mean_high, CW11_mean)
    
    levels <- factor(c("variance low", "variance high", "variance", "mean low", "mean high", "mean"))
    MMCACZ <- data.frame(CW1_vector, CW2_vector,CW3_vector,CW4_vector,CW5_vector,CW6_vector,CW7_vector,CW8_vector,CW9_vector,CW10_vector,CW11_vector, f = levels)
    MMCACZ <- t(MMCACZ)
    MMCACZ <- MMCACZ[-12,]
    MMCAWS <- as.data.frame(MMCAWS)
    MMCACZ$percentage <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    MMCACZ$group <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) 
    MMCACZ <- transform(MMCACZ, V5 = as.numeric(as.character(V5)))
    
    
    # Plot
    plot <- ggplot(MMCACZ, 
                   aes(x = percentage, y = V6, group = 1)) +
      geom_point(size = 5, shape = 18) + 
      geom_line(size = 0.7)  + 
      geom_ribbon(mapping = aes(ymin = V4, ymax =V5), 
                  alpha = 0.4, 
                  color = "tomato", 
                  fill = "tomato")
    plot1 <- plot + 
      theme_minimal() + labs(x = "Percentage (%)", y = "Variance" ) +
      theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = 'none',
            panel.background = element_rect(color = "grey40", size = 0.5),
            axis.text = element_text(face = "italic", size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 14, face = "bold", vjust = 1))
  
# Integration and modularity
    
  # Covariance analyses : https://www.mail-archive.com/morphmet@morphometrics.org/msg02388/RSfunctionR.R
    
    library("geomorph", lib.loc="C:/Program Files/R/R-3.3.0/library")
    library("dplyr", lib.loc="C:/Program Files/R/R-3.3.0/library")
    library("evolqg", lib.loc="C:/Program Files/R/R-3.3.0/library")
    
    # load or calculate procrustes coordinates, convert to data frame
    
    # separate out strains: use split instead?

    CACZF1 <- filter(CACZ, Genotype_Treatment == "CastxCzech")
    CAST <- filter(CACZ, Genotype_Treatment == "CAST-EiJ")
    CZECHI <- filter(CACZ, Genotype_Treatment == "CZECHI-EiJ")
    CACZF2 <- filter(CACZ, Genotype_Treatment == "CastxCzech*F2")
    CACZB1 <- filter(CACZ, Genotype_Treatment == "CastxCzech*Czech")
    
    # Calculate appropriate data tables (landmarks only)
    CACZF1 <- as.matrix(CACZF1[,17:ncol(CACZF1)])
    CAST <- as.matrix(CAST[,17:ncol(CAST)]) 
    CZECHI <- as.matrix(CZECHI[,17:ncol(CZECHI)]) 
    CACZF2 <- as.matrix(CACZF2[,17:ncol(CACZF2)])
    CACZB1 <- as.matrix(CACZB1[,17:ncol(CACZB1)])
    
    # Calculate V/CV matrix for each strain
    CACZF1cov <- cov(CACZF1, use = "complete.obs") 
    CASTcov <- cov(CAST, use = "complete.obs")
    CZECHIcov <- cov(CZECHI, use = "complete.obs")
    CACZF2cov <-cov(CACZF2, use = "complete.obs")
    CACZB1cov <- cov(CACZB1, use = "complete.obs")
    
    # Matrix repeatability for each strain
    #either:
    BootstrapRep(CAST, ComparisonFunc = RandomSkewers, iterations = 10)
    
    #Compare Matrices
    
    reps <- unlist(lapply(list(CASTcov, CZECHIcov, CACZF1cov, CACZF2cov, CACZB1cov), MonteCarloRep, sample.size = 30,
                          RandomSkewers, num.vectors = 999, 
                          iterations = 999))
    
    RandomSkewers(list(CASTcov, CZECHIcov, CACZF1cov, CACZF2cov, CACZB1cov), repeat.vector = reps)
    
    MatrixCompare(CACZF1cov, CASTcov) # EVOLQ package: Random Skewers, Mantel Cor, KrzCor and PCA similarity
    
    cor.list <- llply(list(CASTcov, CZECHIcov, CACZF1cov, CACZF2cov, CACZB1cov), cov2cor)
    MantelCor(cor.list)
    