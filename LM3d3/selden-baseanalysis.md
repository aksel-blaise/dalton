Analysis of Dalton point morphology
================
Robert Z. Selden, Jr.
18 September, 2020

## Load packages + data

``` r
# load packages

# devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
library(geomorph)
```

    ## Loading required package: RRPP

    ## Loading required package: rgl

``` r
library(tidyverse)
```

    ## -- Attaching packages ---------------------------------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ------------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(ggfortify)
library(ggExtra)
library(wesanderson)

# setwd
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

# print qdata
knitr::kable(qdata, align = "cc", caption = "Attributes included in qdata.")
```

|               | heart.out | heart.reg |
| :------------ | :-------: | :-------: |
| 11AI225       |     H     |     H     |
| 11HE445       |     N     |     P     |
| HK49\_1462    |     N     |     I     |
| HK49\_2       |     N     |     I     |
| HK49\_3022    |     N     |     I     |
| HK49\_3823    |     N     |     I     |
| HK49\_4       |     N     |     I     |
| HK49\_5928    |     N     |     I     |
| HK49\_5961    |     N     |     I     |
| HK49\_7       |     N     |     I     |
| KeilMangold10 |     H     |     H     |
| KeilMangold13 |     H     |     H     |
| KeilMangold14 |     H     |     H     |
| KeilMangold17 |     H     |     H     |
| KeilMangold3  |     H     |     H     |
| KeilMangold8  |     H     |     H     |
| Kinzer46      |     N     |     P     |
| Kinzer49      |     N     |     P     |

Attributes included in qdata.

## Generalized Procrustes Analysis

``` r
# gpa
Y.gpa <- gpagen(coords, PrinAxes = TRUE, ProcD = TRUE, Proj = TRUE, print.progress = FALSE)

# output + consensus configuration coords
Y.gpa
```

    ## 
    ## Call:
    ## gpagen(A = coords, PrinAxes = TRUE, ProcD = TRUE, Proj = TRUE,  
    ##     print.progress = FALSE) 
    ## 
    ## 
    ## 
    ## Generalized Procrustes Analysis
    ## with Partial Procrustes Superimposition
    ## 
    ## 34 fixed landmarks
    ## 0 semilandmarks (sliders)
    ## 3-dimensional landmarks
    ## 2 GPA iterations to converge
    ## 
    ## 
    ## Consensus (mean) Configuration
    ## 
    ##               X             Y             Z
    ## 1  -0.325675032 -0.0005504061 -1.111179e-03
    ## 2   0.213810018  0.1202370207  2.413484e-04
    ## 3   0.208592114 -0.1268614549  1.292567e-03
    ## 4   0.085616980  0.1255693154 -1.582209e-03
    ## 5   0.080617786 -0.1227413320  5.934124e-04
    ## 6   0.193204163 -0.0027172747  6.393959e-04
    ## 7  -0.244087056  0.0483180025  1.953045e-03
    ## 8  -0.162326598  0.0742208861 -4.353625e-05
    ## 9  -0.080502652  0.0879158521  6.725165e-04
    ## 10  0.001346185  0.1015735485 -1.109478e-04
    ## 11  0.142848111  0.1179897042 -1.563014e-03
    ## 12  0.214037803  0.0602732283 -1.075357e-03
    ## 13  0.213443011 -0.0675825998  4.317139e-04
    ## 14  0.133888236 -0.1192393039 -1.590847e-04
    ## 15  0.002428225 -0.0977183027 -2.038832e-03
    ## 16 -0.079563101 -0.0855690877 -1.460510e-03
    ## 17 -0.161573124 -0.0729909602 -1.907247e-03
    ## 18 -0.243601414 -0.0509738329 -1.578314e-03
    ## 19 -0.243844471 -0.0002416158  2.134015e-02
    ## 20 -0.161945592  0.0002587260  2.900944e-02
    ## 21 -0.080058846  0.0006277119  3.124951e-02
    ## 22  0.001884748  0.0011342853  3.016402e-02
    ## 23  0.083092773  0.0015190651  2.661976e-02
    ## 24  0.138008721 -0.0005826870  1.920308e-02
    ## 25  0.138450265 -0.0006650425 -1.942717e-02
    ## 26  0.083123393  0.0013799618 -2.628140e-02
    ## 27  0.002002263  0.0009627902 -2.970941e-02
    ## 28 -0.079926357  0.0005687189 -2.957626e-02
    ## 29 -0.161872411  0.0001077377 -2.656883e-02
    ## 30 -0.243791192 -0.0002436710 -2.060505e-02
    ## 31  0.084427201  0.0670319779  2.337641e-02
    ## 32  0.081753723 -0.0635162822  2.429332e-02
    ## 33  0.081795311 -0.0638104774 -2.313601e-02
    ## 34  0.084396819  0.0663157982 -2.314534e-02

``` r
# plot consensus configuration
par(mfrow=c(1, 3))
plot(Y.gpa$consensus[,c("Y", "X")], pch=20)
plot(Y.gpa$consensus[,c("Z", "X")], pch=20)
plot(Y.gpa$consensus[,c("Z", "Y")], pch=20)
```

<div class="figure">

<img src="selden-baseanalysis_files/figure-gfm/consensus-1.png" alt="2D plot of 3D consensus configuration." width="100%" />

<p class="caption">

2D plot of 3D consensus configuration.

</p>

</div>

``` r
# render 3d gpa plot
#plot(Y.gpa)

# gpa plot
# knitr::include_graphics('images/gpa3d.png')
```

## Set gdf + boxplots by centroid size

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
    ## Number of observations 18 
    ## Number of vectors 18 
    ## 
    ## Importance of Components:
    ##                             Comp1      Comp2        Comp3        Comp4
    ## Eigenvalues            0.01101874 0.00143503 0.0009587457 0.0004452499
    ## Proportion of Variance 0.72284625 0.09414020 0.0628951867 0.0292090783
    ## Cumulative Proportion  0.72284625 0.81698645 0.8798816353 0.9090907137
    ##                               Comp5        Comp6        Comp7        Comp8
    ## Eigenvalues            0.0003533586 0.0002416372 0.0001990909 0.0001465823
    ## Proportion of Variance 0.0231808685 0.0158517698 0.0130606709 0.0096160249
    ## Cumulative Proportion  0.9322715822 0.9481233521 0.9611840230 0.9708000478
    ##                               Comp9       Comp10       Comp11       Comp12
    ## Eigenvalues            0.0000996901 8.139922e-05 7.017007e-05 5.066795e-05
    ## Proportion of Variance 0.0065398237 5.339914e-03 4.603264e-03 3.323895e-03
    ## Cumulative Proportion  0.9773398715 9.826798e-01 9.872830e-01 9.906069e-01
    ##                              Comp13       Comp14       Comp15       Comp16
    ## Eigenvalues            4.155661e-05 3.586756e-05 2.518719e-05 2.179588e-05
    ## Proportion of Variance 2.726177e-03 2.352967e-03 1.652318e-03 1.429843e-03
    ## Cumulative Proportion  9.933331e-01 9.956861e-01 9.973384e-01 9.987683e-01
    ##                              Comp17       Comp18
    ## Eigenvalues            1.877623e-05 2.012194e-33
    ## Proportion of Variance 1.231750e-03 1.320030e-31
    ## Cumulative Proportion  1.000000e+00 1.000000e+00

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
    ##           Df     SS      MS     Rsq      F      Z Pr(>F)  
    ## heart      1 1574.6 1574.65 0.28078 6.2464 1.3532 0.0236 *
    ## Residuals 16 4033.4  252.09 0.71922                       
    ## Total     17 5608.1                                       
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
    ## hreg       2 1606.6 803.29 0.28648 3.0112 1.2005  0.078 .
    ## Residuals 15 4001.5 266.77 0.71352                       
    ## Total     17 5608.1                                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##             d UCL (95%)          Z Pr > d
    ## H:I 18.142628  18.15233  1.9298847 0.0503
    ## H:P 21.968027  24.07614  1.6170807 0.0760
    ## I:P  3.825399  23.66874 -0.8365537 0.7575

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
    ## 263.9717 158.2345 295.9358 
    ## 
    ## Pairwise distances between variances, plus statistics
    ##             d UCL (95%)           Z Pr > d
    ## H:I 105.73722  234.3821  0.10137213 0.3983
    ## H:P  31.96414  310.5803 -1.05245785 0.8554
    ## I:P 137.70136  301.5930  0.08975482 0.4094

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
    ##           Df      SS       MS     Rsq      F      Z Pr(>F)   
    ## heart      1 0.07787 0.077870 0.30049 6.8733 2.6075 0.0054 **
    ## Residuals 16 0.18127 0.011329 0.69951                        
    ## Total     17 0.25914                                         
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
    ## hreg       2 0.090468 0.045234 0.34911 4.0226 2.2609  0.011 *
    ## Residuals 15 0.168672 0.011245 0.65089                       
    ## Total     17 0.259140                                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
    ##             d UCL (95%)          Z Pr > d
    ## H:I 18.142628  18.15233  1.9298847 0.0503
    ## H:P 21.968027  24.07614  1.6170807 0.0760
    ## I:P  3.825399  23.66874 -0.8365537 0.7575

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
    ## 263.9717 158.2345 295.9358 
    ## 
    ## Pairwise distances between variances, plus statistics
    ##             d UCL (95%)           Z Pr > d
    ## H:I 105.73722  234.3821  0.10137213 0.3983
    ## H:P  31.96414  310.5803 -1.05245785 0.8554
    ## I:P 137.70136  301.5930  0.08975482 0.4094
