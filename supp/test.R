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
