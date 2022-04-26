#mesh processing
#uniform remesh and decimation
#output as PLY for landmarking and OFF for auto3dgm
#devtools::install_github("zarquon42b/Rvcg")
library(Rvcg)
library(rgl)

#set working directory
setwd(getwd())
###################

dir.create("ply")
dir.create("data3d/lowres", recursive = T)

processmesh <- function(x) {
  mesh <- vcgImport(x, updateNormals = TRUE, clean = TRUE, silent = FALSE)
  meshname <- gsub(".stl", "", basename(x))
  #check for existance and validity of vertices, faces and vertex normals
  meshintegrity(mesh, facecheck = TRUE, normcheck = TRUE)
  #identify resolution of mesh for use as voxelSize in remesh
  mres<-vcgMeshres(mesh)
  #uniform remesh
  remesh<-vcgUniformRemesh(mesh, voxelSize = mres$res, offset = 0)
  #export as OFF to data3d file (auto3dgm)
  vcgOffWrite(remesh, filename = paste0("data3d/", meshname))
  #export as PLY to ply folder (landmarking in Design X)
  vcgPlyWrite(remesh, filename = "ply/", meshname, ascii = TRUE)
  #decimate mesh
  decimated<-vcgQEdecim(remesh, percent = 0.25)
  #export as OFF to lowres file (auto3dgm)
  vcgOffWrite(decimated, filename = paste0("data3d/lowres/", meshname))
}

#batch loop assuming all stl files are stored in a folder called stl
mystls <- list.files("stl", pattern=".stl", full.names = T)
runall <- lapply(mystls, processmesh)
#end of script
