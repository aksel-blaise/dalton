Analysis of Dalton point morphology
================
Robert Z. Selden, Jr.
17 September, 2020

## Generalised Procrustes Analysis

``` r
# library(devtools)
#devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
library(geomorph)
```

    ## Loading required package: RRPP

    ## Loading required package: rgl

``` r
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ---------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(wesanderson)
setwd(getwd())

# read GM data
source('readmulti.csv.R')

# read .csv files
setwd("./data")
filelist <- list.files(pattern = ".csv")
coords <- readmulti.csv(filelist)
setwd("../")

# read qualitative data
qdata <- read.csv("qdata.csv", header = TRUE, row.names = 1)
qdata <- qdata[match(dimnames(coords)[[3]],rownames(qdata)),]
qdata
```

    ##               heart.out heart.reg
    ## 11AI225               H         H
    ## 11HE445               N         P
    ## HK49_1462             N         I
    ## HK49_2                N         I
    ## HK49_3022             N         I
    ## HK49_3823             N         I
    ## HK49_4                N         I
    ## HK49_7                N         I
    ## KeilMangold13         H         H
    ## KeilMangold3          H         H
    ## Kinzer46              N         P
    ## Kinzer49              N         P

``` r
# gpa
Y.gpa <- gpagen(coords, PrinAxes = TRUE, ProcD = TRUE, Proj = TRUE, print.progress = FALSE)
# plot consensus configuration
par(mfrow=c(1,3))
plot(Y.gpa$consensus[,c("Y", "X")], pch=20)
plot(Y.gpa$consensus[,c("Z", "X")], pch=20)
plot(Y.gpa$consensus[,c("Z", "Y")], pch=20)
```

<img src="selden-baseanalysis_files/figure-gfm/gpa-1.png" width="100%" />

``` r
# gpa plot
# knitr::include_graphics('images/gpa3d.png')
# fig.cap="Results of generalized Procrustes analysis."
```

``` r
# geomorph data frame
gdf <- geomorph.data.frame(shape = Y.gpa$coords, size = Y.gpa$Csize, heart = qdata$heart.out, hreg = qdata$heart.reg)

# attributes for boxplots
csz <- Y.gpa$Csize # centroid size
heart <- qdata$heart.out # heartland in/out
hreg <- qdata$heart.reg # heartland region

# boxplot of Dalton point centroid size by in/out heartland
boxplot(csz~heart,
        names = c("H","N"), # heartland (H), and not heartland (N)
        xlab = "Heartland",
        ylab = "Centroid Size",
        col = wes_palette("Moonrise2"),
        )
```

<img src="selden-baseanalysis_files/figure-gfm/data.frame-1.png" width="100%" />

``` r
fig.cap = "Boxplot of centroid size by Heartland (in/out)."

# boxplot of Dalton point centroid size by heartland + regions
boxplot(csz~hreg,
        names = c("H","I","P"), # heartland (H), interior (I), and northern periphery (P)
        xlab = "Heartland Region",
        ylab = "Centroid Size",
        col = wes_palette("Moonrise2"),
        )
```

<img src="selden-baseanalysis_files/figure-gfm/data.frame-2.png" width="100%" />

``` r
fig.cap = "Boxplot of centroid size by Heartland region."
```

## Principal Components Analysis

``` r
# principal components analysis
pca<-gm.prcomp(Y.gpa$coords)
summary(pca)
```

    ## 
    ## Ordination type: Principal Component Analysis 
    ## Centering and projection: OLS 
    ## Number of observations 12 
    ## Number of vectors 12 
    ## 
    ## Importance of Components:
    ##                              Comp1       Comp2        Comp3       Comp4
    ## Eigenvalues            0.008924199 0.001495728 0.0008987617 0.000579829
    ## Proportion of Variance 0.671963425 0.112623520 0.0676738607 0.043659251
    ## Cumulative Proportion  0.671963425 0.784586944 0.8522608052 0.895920056
    ##                              Comp5        Comp6       Comp7        Comp8
    ## Eigenvalues            0.000404076 0.0003139634 0.000291867 0.0001797937
    ## Proportion of Variance 0.030425620 0.0236404304 0.021976646 0.0135378841
    ## Cumulative Proportion  0.926345676 0.9499861068 0.971962753 0.9855006375
    ##                               Comp9       Comp10       Comp11       Comp12
    ## Eigenvalues            8.108541e-05 6.602839e-05 4.544906e-05 2.006363e-33
    ## Proportion of Variance 6.105470e-03 4.971725e-03 3.422168e-03 1.510727e-31
    ## Cumulative Proportion  9.916061e-01 9.965778e-01 1.000000e+00 1.000000e+00

``` r
# set plot parameters to plot by heartland in (H) and out (N)
pch.gps.heart <- c(15,17)[as.factor(heart)]
col.gps.heart <- wes_palette("Moonrise2")[as.factor(heart)]
col.hull <- c("#798E87","#C27D38")

# plot pca by heartland in (H) and out (N)
pc.plot1 <- plot(pca, 
                 asp = 1,
                 pch = pch.gps.heart,
                 col = col.gps.heart)
                    shapeHulls(pc.plot1, 
                             groups = heart,
                             group.cols = col.hull)
```

<img src="selden-baseanalysis_files/figure-gfm/pca-1.png" width="100%" />

``` r
# set plot parameters to plot by heartland + regions
pch.gps.hreg <- c(15,17,18)[as.factor(hreg)]
col.gps.hreg <- wes_palette("Moonrise2")[as.factor(hreg)]
col.hull.2 <- c("#798E87","#CCC591","#C27D38")

# plot pca by heartland + regions
pc.plot2 <- plot(pca, 
                 asp = 1,
                 pch = pch.gps.hreg,
                 col = col.gps.hreg)
                    shapeHulls(pc.plot2, 
                             groups = hreg,
                             group.cols = col.hull.2)
```

<img src="selden-baseanalysis_files/figure-gfm/pca-2.png" width="100%" />

## Define models

``` r
## Define models
# size as a function of heart
fit.size.heart <- procD.lm(size ~ heart, 
                           data = gdf, 
                           print.progress = FALSE, 
                           iter = 9999)
# size as a function of hreg
fit.size.hreg <- procD.lm(size ~ hreg, 
                          data = gdf, 
                          print.progress = FALSE, 
                          iter = 9999)

# shape as a function of heart
fit.shape.heart <- procD.lm(shape ~ heart, 
                            data = gdf, 
                            print.progress = FALSE, 
                            iter = 9999)
# shape as a function of hreg
fit.shape.hreg <- procD.lm(shape ~ hreg, 
                           data = gdf, 
                           print.progress = FALSE, 
                           iter = 9999)
```

## Size

``` r
# ANOVA: do dalton projectile point sizes differ by heart?
anova(fit.size.heart)
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 10000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##           Df     SS      MS     Rsq      F      Z  Pr(>F)  
    ## heart      1 1348.6 1348.63 0.36056 5.6386 1.3779 0.04505 *
    ## Residuals 10 2391.8  239.18 0.63944                        
    ## Total     11 3740.4                                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Call: procD.lm(f1 = size ~ heart, iter = 9999, data = gdf, print.progress = FALSE)

``` r
# ANOVA: do dalton projectile point sizes differ by hreg?
anova(fit.size.hreg)
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 10000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##           Df     SS     MS     Rsq      F      Z Pr(>F)
    ## hreg       2 1353.8 676.90 0.36194 2.5526 1.0291 0.1384
    ## Residuals  9 2386.6 265.18 0.63806                     
    ## Total     11 3740.4                                    
    ## 
    ## Call: procD.lm(f1 = size ~ hreg, iter = 9999, data = gdf, print.progress = FALSE)

``` r
# pairwise comparison of LS means = which differ?
sz.hreg <- pairwise(fit.size.hreg, 
                    groups = qdata$heart.reg)
summary(sz.hreg, 
        confidence = 0.95, 
        test.type = "dist")
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: H I P 
    ## 
    ## RRPP: 10000 permutations
    ## 
    ## LS means:
    ## Vectors hidden (use show.vectors = TRUE to view)
    ## 
    ## Pairwise distances between means, plus statistics
    ##             d UCL (95%)         Z Pr > d
    ## H:I 23.946604  24.88171  1.769518 0.0623
    ## H:P 25.554212  28.39577  1.526889 0.0910
    ## I:P  1.607608  24.99907 -1.192242 0.9057

``` r
# pairwise distance between variances = standardization?
summary(sz.hreg, 
        confidence = 0.95, 
        test.type = "var")
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: H I P 
    ## 
    ## RRPP: 10000 permutations
    ## 
    ## 
    ## Observed variances by group
    ## 
    ##        H        I        P 
    ## 173.8491 162.8750 295.9358 
    ## 
    ## Pairwise distances between variances, plus statistics
    ##             d UCL (95%)           Z  Pr > d
    ## H:I  10.97416  291.3134 -1.34286166 0.95165
    ## H:P 122.08670  342.0968 -0.24363463 0.53445
    ## I:P 133.06086  296.4747  0.03641097 0.43275

## Shape

``` r
# ANOVA: do dalton projectile point shapes differ by heart?
anova(fit.shape.heart)
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 10000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##           Df       SS       MS     Rsq      F      Z Pr(>F)  
    ## heart      1 0.032292 0.032292 0.22104 2.8376 1.6911 0.0654 .
    ## Residuals 10 0.113797 0.011380 0.77896                       
    ## Total     11 0.146089                                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Call: procD.lm(f1 = shape ~ heart, iter = 9999, data = gdf, print.progress = FALSE)

``` r
# ANOVA: do dalton projectile point shapes differ by hreg?
anova(fit.shape.hreg)
```

    ## 
    ## Analysis of Variance, using Residual Randomization
    ## Permutation procedure: Randomization of null model residuals 
    ## Number of permutations: 10000 
    ## Estimation method: Ordinary Least Squares 
    ## Sums of Squares and Cross-products: Type I 
    ## Effect sizes (Z) based on F distributions
    ## 
    ##           Df       SS       MS     Rsq      F      Z Pr(>F)
    ## hreg       2 0.040505 0.020253 0.27727 1.7264 1.0623 0.1623
    ## Residuals  9 0.105583 0.011732 0.72273                     
    ## Total     11 0.146089                                      
    ## 
    ## Call: procD.lm(f1 = shape ~ hreg, iter = 9999, data = gdf, print.progress = FALSE)

``` r
# pairwise comparison of LS means = which differ?
sh.hreg <- pairwise(fit.shape.hreg, groups = qdata$heart.reg)
summary(sz.hreg, confidence = 0.95, test.type = "dist")
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: H I P 
    ## 
    ## RRPP: 10000 permutations
    ## 
    ## LS means:
    ## Vectors hidden (use show.vectors = TRUE to view)
    ## 
    ## Pairwise distances between means, plus statistics
    ##             d UCL (95%)         Z Pr > d
    ## H:I 23.946604  24.88171  1.769518 0.0623
    ## H:P 25.554212  28.39577  1.526889 0.0910
    ## I:P  1.607608  24.99907 -1.192242 0.9057

``` r
# pairwise distance between variances = standardization?
summary(sz.hreg, confidence = 0.95, test.type = "var")
```

    ## 
    ## Pairwise comparisons
    ## 
    ## Groups: H I P 
    ## 
    ## RRPP: 10000 permutations
    ## 
    ## 
    ## Observed variances by group
    ## 
    ##        H        I        P 
    ## 173.8491 162.8750 295.9358 
    ## 
    ## Pairwise distances between variances, plus statistics
    ##             d UCL (95%)           Z  Pr > d
    ## H:I  10.97416  291.3134 -1.34286166 0.95165
    ## H:P 122.08670  342.0968 -0.24363463 0.53445
    ## I:P 133.06086  296.4747  0.03641097 0.43275
