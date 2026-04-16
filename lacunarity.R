# ---- Packages ----
req <- c("lidR", "terra", "dplyr", "data.table", "lacunr")
to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(req, library, character.only = TRUE)



data_understory="C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-fullprocess/GER10_sample_understory.las"


las <- readLAS(data_understory)
# voxelize the LAS point cloud, taking care to input the correct S4 slot
vox <- voxelize(las@data, edge_length = c(0.25, 0.25, 0.25))

# generate binary map
box <- bounding_box(vox)
# calculate lacunarity curve
lac_curve <- lacunarity(box)

plot(lac_curve$box_size,lac_curve$lacunarity)
