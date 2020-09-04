# library(devtools)
devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
library(geomorph)
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
qdata <- read.csv("qdata.csv",header=TRUE,row.names=1)
qdata <- qdata[match(dimnames(coords)[[3]],rownames(qdata)),]
qdata

# gpa
Y.gpa <- gpagen(coords, PrinAxes = TRUE, ProcD = TRUE, Proj = TRUE, print.progress = FALSE)
plot(Y.gpa)

# gpa plot
# knitr::include_graphics('images/gpa3d.png')
# fig.cap="Results of generalized Procrustes analysis."

# geomorph data frame
gdf <- geomorph.data.frame(shape = Y.gpa$coords, size = Y.gpa$Csize, heart = qdata$heart.out, hreg = qdata$heart.reg)

# attributes for boxplots
csz <- Y.gpa$Csize # centroid size
heart <- qdata$heart.out # heartland in/out
hreg <- qdata$heart.reg # heartland region

# boxplot of Dalton point centroid size by in/out Heartland
boxplot(csz~heart,
        names = c("H","N"),
        xlab = "Heartland (H) or not (N)",
        ylab = "Centroid Size",
        col = wes_palette("Moonrise2"),
)
fig.cap = "Boxplot of centroid size by Heartland (in/out)."

# boxplot of Dalton point centroid size by Heartland Region
boxplot(csz~hreg,
        names = c("H","I","P"),
        xlab = "Heartland Region",
        ylab = "Centroid Size",
        col = wes_palette("Moonrise2"),
)
fig.cap = "Boxplot of centroid size by Heartland region."

# principal components analysis
pca<-gm.prcomp(Y.gpa$coords)
summary(pca)
# set plot parameters to plot by Heartland in/out
pch.gps.heart <- c(15,17)[as.factor(heart)]
col.gps.heart <- wes_palette("Moonrise2")[as.factor(heart)]
col.hull <- c("#798E87","#C27D38")
# plot pca by incision profile
pc.plot1 <- plot(pca, asp = 1,
                 pch = pch.gps.heart,
                 col = col.gps.heart)
shapeHulls(pc.plot1, 
           groups = heart,
           group.cols = col.hull)

# set plot parameters to plot by Heartland region
pch.gps.hreg <- c(15,17,18)[as.factor(hreg)]
col.gps.hreg <- wes_palette("Moonrise2")[as.factor(hreg)]
col.hull.2 <- c("#798E87","#C27D38","#CCC591")
# plot pca by incision profile (inc2)
pc.plot2 <- plot(pca, asp = 1,
                 pch = pch.gps.hreg,
                 col = col.gps.hreg)
shapeHulls(pc.plot2, 
           groups = hreg,
           group.cols = col.hull.2)
