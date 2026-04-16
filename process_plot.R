library(ITSMe)
library(lidR)
library(dplyr)
library(data.table)
library(terra)

plot="GER02"
path_plot="C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project"
path_folder=file.path(path_plot,
                      paste0(plot,"-processing"))
path_trees=file.path(path_folder,"trees")
subfolders=c("adults")#,"regen","hesitation","dead-incomplete","tree-parts")

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

#normalize cp
# Use a low percentile instead of the absolute minimum

# plot(dtm_rast)
las_norm <- normalize_height(las, dtm_rast)
# plot(las_norm)
# plot(las)

### DTM LAS ####
#%%%%%%%%%%%%%%%
dtm_rast <- rast(file.path(path_folder,
                           list.files(path_folder,pattern="dtm.tif")))
las <- clip_rectangle(las,
                      xleft   = ext(dtm_rast)$xmin,
                      ybottom = ext(dtm_rast)$ymin,
                      xright  = ext(dtm_rast)$xmax,
                      ytop    = ext(dtm_rast)$ymax)
las@data$X <- las@data$X - min(las@data$X)
las@data$Y <- las@data$Y - min(las@data$Y)
las@data$Z <- las@data$Z - min(las@data$Z)
las <- las_update(las)

las_norm <- normalize_height(las, dtm_rast)
plot(las_norm)



# Slice

slice_0_1 <- filter_poi(las_norm, Z >= 0 & Z < 1)
plot(slice_0_1)
slice_0_1_chm <- rasterize_canopy(slice_0_1, res = 0.05, algorithm = p2r())
plot(slice_0_1_chm)
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
