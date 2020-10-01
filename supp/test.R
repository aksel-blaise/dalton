# load packages

# devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
library(geomorph)
library(tidyverse)
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
qdata <- qdata[match(dimnames(coords)[[3]], rownames(qdata)),]

###############################################    NEED HELP HERE    ##########

# select landmarks from coords for cross-sections
a <- c(1,6,19:30) # midline landmarks
b <- c(4,5,31:34) # base/body intersection

# new coords midline
mid.cross <- coords[a,,]
arrayspecs(mid.cross, 14, 3)

# new coords blade/base
bb.cross <- coords[b,,]
arrayspecs(bb.cross, 6, 3)

#############################################################################

# gpa
Y.gpa <- gpagen(mid.cross, 
                PrinAxes = TRUE, 
                ProcD = TRUE, 
                Proj = TRUE, 
                print.progress = FALSE)

# output + consensus configuration coords
Y.gpa

# symmetric shape/design intent
symm.shape <- Y.gpa$mid.cross[,1:2,]
symm.shape2 <- res.bilat$symm.shape[,2:3,]

# plot all specimens/2D (symmetric component)
plotAllSpecimens(symm.shape)
plotAllSpecimens(symm.shape2)

# geomorph data frame
gdf <- geomorph.data.frame(shape = Y.gpa$mid.cross, 
                           size = Y.gpa$Csize, 
                           heart = qdata$heart.out, 
                           hreg = qdata$heart.reg,
                           bev.1 = qdata$bev, 
                           bev.2 = qdata$bev.type)


