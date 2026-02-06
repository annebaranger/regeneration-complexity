# ---- Packages ----
req <- c("lidR", "terra", "dplyr", "data.table", "spdep", "gstat", "sf")
lapply(req, library, character.only = TRUE)

path_data="C:/Users/amob2/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-regeneration-seg/"
path_regen=file.path(path_data,"understory_12.las")
path_dtm=file.path(path_data,"GER10_dtm.las")


LASfile=readLAS(path_regen)

### canopy height model & rugosity
dtm_las=readLAS(path_dtm)
understory=rbindlist(list(as.data.table(dtm_las@data), as.data.table(LASfile@data)), use.names = TRUE, fill = TRUE)
las_understory=LAS(understory)
plot(las_understory)

chm <- rasterize_canopy(las_understory, res = 0.2,p2r(na.fill = tin()))
col <- height.colors(25)
plot(chm, col = col)
