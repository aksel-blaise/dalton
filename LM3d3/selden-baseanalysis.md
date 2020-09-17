Analysis of Dalton point morphology
================
Robert Z. Selden, Jr.
17 September, 2020

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

    ## -- Attaching packages --------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
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
knitr::kable(qdata, caption = "Attributes included in qdata.")
```

|               | heart.out | heart.reg |
| :------------ | :-------- | :-------- |
| 11AI225       | H         | H         |
| 11HE445       | N         | P         |
| HK49\_1462    | N         | I         |
| HK49\_2       | N         | I         |
| HK49\_3022    | N         | I         |
| HK49\_3823    | N         | I         |
| HK49\_4       | N         | I         |
| HK49\_5928    | N         | I         |
| HK49\_7       | N         | I         |
| KeilMangold13 | H         | H         |
| KeilMangold3  | H         | H         |
| Kinzer46      | N         | P         |
| Kinzer49      | N         | P         |

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
    ##                X             Y             Z
    ## 1  -0.3210737107 -1.204034e-04 -0.0008305918
    ## 2   0.2112947896  1.292704e-01 -0.0005340252
    ## 3   0.2071159485 -1.350098e-01  0.0016317622
    ## 4   0.0841677848  1.314787e-01 -0.0021968117
    ## 5   0.0785854204 -1.291301e-01  0.0013684320
    ## 6   0.1916859452 -2.534104e-03  0.0001599942
    ## 7  -0.2407227573  5.116942e-02  0.0030167704
    ## 8  -0.1602405966  7.874551e-02  0.0005143010
    ## 9  -0.0796906615  9.211458e-02  0.0016416122
    ## 10  0.0008834463  1.048382e-01  0.0001825901
    ## 11  0.1409639696  1.246187e-01 -0.0014128103
    ## 12  0.2136621313  6.503498e-02 -0.0010568461
    ## 13  0.2134795801 -7.171211e-02 -0.0001156529
    ## 14  0.1324044620 -1.258878e-01  0.0005509174
    ## 15  0.0017247893 -1.015790e-01 -0.0013896936
    ## 16 -0.0789309051 -9.001188e-02 -0.0022791408
    ## 17 -0.1596308556 -7.712807e-02 -0.0041430449
    ## 18 -0.2403705138 -5.420849e-02 -0.0029763169
    ## 19 -0.2405357398  7.184236e-05  0.0220410725
    ## 20 -0.1599331669  4.710396e-04  0.0294696030
    ## 21 -0.0793313725  7.178825e-04  0.0319430155
    ## 22  0.0013105903  1.100532e-03  0.0317543878
    ## 23  0.0813827170  1.334253e-03  0.0268172007
    ## 24  0.1364237320 -5.711318e-04  0.0188623315
    ## 25  0.1367435832 -7.316770e-04 -0.0188474304
    ## 26  0.0813389144  1.080900e-03 -0.0262606392
    ## 27  0.0014381682  7.899276e-04 -0.0302753065
    ## 28 -0.0792025637  5.317840e-04 -0.0298185936
    ## 29 -0.1598505579  2.159484e-04 -0.0272494484
    ## 30 -0.2404901120  1.221137e-05 -0.0210724777
    ## 31  0.0828785991  6.975180e-02  0.0237862766
    ## 32  0.0798787244 -6.673047e-02  0.0250819936
    ## 33  0.0798675097 -6.733899e-02 -0.0234053711
    ## 34  0.0827727083  6.934535e-02 -0.0249580596

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
    ## Number of observations 13 
    ## Number of vectors 13 
    ## 
    ## Importance of Components:
    ##                              Comp1       Comp2       Comp3        Comp4
    ## Eigenvalues            0.009725919 0.001890297 0.000982884 0.0005426515
    ## Proportion of Variance 0.666057230 0.129452676 0.067310554 0.0371622407
    ## Cumulative Proportion  0.666057230 0.795509906 0.862820460 0.8999827010
    ##                               Comp5        Comp6        Comp7        Comp8
    ## Eigenvalues            0.0004408853 0.0003052945 0.0002688807 0.0001665939
    ## Proportion of Variance 0.0301930155 0.0209073957 0.0184136754 0.0114088005
    ## Cumulative Proportion  0.9301757165 0.9510831122 0.9694967876 0.9809055881
    ##                               Comp9       Comp10       Comp11       Comp12
    ## Eigenvalues            0.0001122693 0.0000684058 5.661008e-05 4.153578e-05
    ## Proportion of Variance 0.0076885040 0.0046846141 3.876811e-03 2.844483e-03
    ## Cumulative Proportion  0.9885940921 0.9932787062 9.971555e-01 1.000000e+00
    ##                              Comp13
    ## Eigenvalues            1.503623e-33
    ## Proportion of Variance 1.029721e-31
    ## Cumulative Proportion  1.000000e+00

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
    ## heart      1 1180.6 1180.58 0.30387 4.8017 1.2302 0.0603 .
    ## Residuals 11 2704.6  245.87 0.69613                       
    ## Total     12 3885.1                                       
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
    ##           Df     SS     MS     Rsq     F       Z Pr(>F)
    ## hreg       2 1217.5 608.76 0.31338 2.282 0.95063 0.1665
    ## Residuals 10 2667.6 266.76 0.68662                     
    ## Total     12 3885.1                                    
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
    ## H:I 21.359973  23.95783  1.5385882 0.0887
    ## H:P 25.554212  28.43570  1.5647865 0.0833
    ## I:P  4.194239  23.91796 -0.8271761 0.7550

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
    ## 173.8491 179.7511 295.9358 
    ## 
    ## Pairwise distances between variances, plus statistics
    ##              d UCL (95%)          Z Pr > d
    ## H:I   5.901927  304.5800 -1.4632414 0.9778
    ## H:P 122.086700  366.2114 -0.3304380 0.5614
    ## I:P 116.184772  300.8819 -0.2282699 0.5327

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
    ## heart      1 0.025824 0.025824 0.14737 1.9013 1.1442 0.1518
    ## Residuals 11 0.149403 0.013582 0.85263                     
    ## Total     12 0.175227                                      
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
    ##           Df       SS       MS     Rsq      F       Z Pr(>F)
    ## hreg       2 0.039761 0.019880 0.22691 1.4676 0.79419 0.2262
    ## Residuals 10 0.135466 0.013547 0.77309                      
    ## Total     12 0.175227                                       
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
    ## H:I 21.359973  23.95783  1.5385882 0.0887
    ## H:P 25.554212  28.43570  1.5647865 0.0833
    ## I:P  4.194239  23.91796 -0.8271761 0.7550

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
    ## 173.8491 179.7511 295.9358 
    ## 
    ## Pairwise distances between variances, plus statistics
    ##              d UCL (95%)          Z Pr > d
    ## H:I   5.901927  304.5800 -1.4632414 0.9778
    ## H:P 122.086700  366.2114 -0.3304380 0.5614
    ## I:P 116.184772  300.8819 -0.2282699 0.5327
