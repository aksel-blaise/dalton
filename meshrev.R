#mesh processing
#uniform remesh and decimation
#output as PLY and OFF for analysis
# devtools::install_github("zarquon42b/Rvcg")
library(Rvcg)
library(rgl)
#set working directory
setwd(getwd())
###################

dir.create("ply")
dir.create("data3d/lowres", recursive = T)

processmesh <- function(x) {
  mesh <- vcgImport(x, updateNormals = TRUE, clean = TRUE, silent = FALSE)
  meshname <- gsub(".wrl", "", basename(x))
  #export as OFF to data3d file (auto3dgm)
  vcgOffWrite(remesh, filename = paste0("data3d/", meshname))
  #export as PLY to ply folder (Design X)
  vcgPlyWrite(remesh, filename = "ply/", meshname, ascii = TRUE)
  #decimate mesh
  decimated<-vcgQEdecim(remesh, percent = 0.25)
  #export as OFF to lowres file (auto3dgm)
  vcgOffWrite(decimated, filename = paste0("data3d/lowres/", meshname))
}

#batch loop assuming all wrl files are stored in a folder called wrl
mywrls <- list.files("wrl", pattern=".wrl", full.names = T)
runall <- lapply(mywrls, processmesh)
#end of script
