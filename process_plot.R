library(ITSMe)
library(lidR)
library(dplyr)
library(data.table)
library(terra)

plot="GER02"
user="amob2"
path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
path_folder=file.path(path_plot,
                      paste0(plot,"-processing"))
path_trees=file.path(path_folder,"trees")
subfolders=c("adults","regen","hesitation","dead-incomplete","tree-parts")

### Create one point cloud ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txt_list=c()
for (subf in subfolders){
  if(subf!="dead-incomplete"){
    txt_list=c(txt_list,
               file.path(file.path(path_trees,subf),
                         list.files(file.path(path_trees,subf), pattern = ".txt")))
  }
}

cloud_list <- lapply(txt_list, function(f) {
  df <- fread(f)[,1:3]
  colnames(df)=c("X","Y","Z")  # adjust sep/header to match your file format
  df
})

# Merge into one data frame
las <- LAS(do.call(rbind, cloud_list))
rm(cloud_list)


### DTM LAS ####
#%%%%%%%%%%%%%%%
dtm_rast <- rast(file.path(path_folder,
                           list.files(path_folder,pattern="dtm.tif")))
las <- clip_rectangle(las,
                      xleft   = ext(dtm_rast)$xmin,
                      ybottom = ext(dtm_rast)$ymin,
                      xright  = ext(dtm_rast)$xmax,
                      ytop    = ext(dtm_rast)$ymax)
las_norm <- normalize_height(las, dtm_rast)

las_norm@data$X <- floor((las_norm@data$X - min(las_norm@data$X))/0.25)
las_norm@data$Y <- floor((las_norm@data$Y - min(las_norm@data$Y))/0.25)
# las@data$Z <- las@data$Z - min(las@data$Z)
las_norm <- las_update(las_norm)

# save(las_norm,file=paste0(plot_name,"_pc.rdata"))

# Slice
load("GER02_rb.rdata")
plot(rast(arr_norm[,,3/0.25]))


slice_0_1 <- filter_poi(las_norm, Z >= 0 & Z < 3)
plot(slice_0_1)
slice_0_1_chm <- rasterize_canopy(slice_0_1, res = 0.0625, algorithm = p2r())
plot(slice_0_1_chm)
as.data.frame(rast(arr_norm[,,3/0.25]),xy=TRUE)
as.data.frame(slice_0_1_chm,xy=TRUE)


slice_1_2 <- filter_poi(las, Z >= 1 & Z < 2)
plot(slice_1_2)
slice_2_3 <- filter_poi(las, Z >= 2 & Z < 3)
plot(slice_2_3)
slice_3_4 <- filter_poi(las, Z >= 3 & Z < 4)
plot(slice_3_4)
slice_3_4_chm <- rasterize_canopy(slice_3_4, res = 0.25, algorithm = p2r())
plot(slice_3_4_chm)

subplot1_csv=fread(file=file.path(path_trees,
                                  "subplot_000000.txt"))[,1:3]
colnames(subplot1_csv)=c("X","Y","Z")  

subplot1=LAS(subplot1_csv)
plot(subplot1)

# Read the point cloud file from its' specified path
tree_pc <- read_tree_pc(path = file.path(path_trees,
                                         "subplot_000000.txt"))
XYZ_pos <- tree_position_pc(pc = tree_pc)


plot(filter_poi(las, Z >= min(las@data$Z) & Z <min(las@data$Z)+3),axis=TRUE)
plot(las,axis=TRUE)
plot(rasterize_canopy(filter_poi(las_norm, Z >= min(las_norm@data$Z) & Z <min(las_norm@data$Z)+2),res=0.25))

plot(filter_poi(las_norm, Z >= min(las_norm@data$Z) & Z <min(las_norm@data$Z)+5),axis=TRUE)
