# auto3dgm
# load required packages
library(Matrix)
library(clue)
library(linprog)
library(igraph)
library(MASS)

# load auto3dgm
rm(list=ls())
library(auto3dgm)

# set working directory
setwd(getwd())

#set data and output directories
Data_dir = "data3d"
Output_dir = "auto3dout"

Levels=c(150,500)
Ids = c('11AI225FSICAA119-107','11HE445-WISD-Bag800-31','HK49_2','HK49_4','HK49_7','HK49_1462','HK49_3022','HK49_3823','HK49_5928',
        'HK49_5961','Keil-Mangold3','KeilMangold8','KeilMangold10','KeilMangold13','KeilMangold14','KeilMangold17','KeilMangold24',
        'KeilMangold42','KeilMangold64','KeilMangold90','Kinzer46','Kinzer49','Kinzer50','Leprechaun11MS1983FSI800-5','MA1699p1',
        'MA1699p3','McL4_158','McL4_614','McL7_101','McL7_168','Nochta11MS12BFIS800-168','Nochta11MS128BFSI267-1','Nochta11MS128bFSIPP1001',
        'Nochta11MS128BFSIPP2495','Nochta11MS128BFSIPP3305','OH1_176','ReedVoss2','TR10_4','TR10_19','Welton1','Welton2','Welton3',
        'Welton4','Welton5','Welton6')  
Names = c('11AI225FSICAA119-107','11HE445-WISD-Bag800-31','HK49_2','HK49_4','HK49_7','HK49_1462','HK49_3022','HK49_3823','HK49_5928',
          'HK49_5961','Keil-Mangold3','KeilMangold8','KeilMangold10','KeilMangold13','KeilMangold14','KeilMangold17','KeilMangold24',
          'KeilMangold42','KeilMangold64','KeilMangold90','Kinzer46','Kinzer49','Kinzer50','Leprechaun11MS1983FSI800-5','MA1699p1',
          'MA1699p3','McL4_158','McL4_614','McL7_101','McL7_168','Nochta11MS12BFIS800-168','Nochta11MS128BFSI267-1','Nochta11MS128bFSIPP1001',
          'Nochta11MS128BFSIPP2495','Nochta11MS128BFSIPP3305','OH1_176','ReedVoss2','TR10_4','TR10_19','Welton1','Welton2','Welton3',
          'Welton4','Welton5','Welton6') 

#FULL is a list of 3 returned elements.  The user gets to specify what is returned.
FULL = align_shapes(Data_dir, Output_dir, Levels, Ids, Names)

ds = FULL[[1]]  #the whole shape dataset.
ga_full=FULL[[2]]  #the global alignment 
pa=FULL[[3]] #the pairwise alignments.  

# Each of these three returned items is useful in identifying and fixing misaligned shapes.  