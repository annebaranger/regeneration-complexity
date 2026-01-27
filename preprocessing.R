### Package imports ####
#%%%%%%%%%%%%%%%%%%%%%%%
library(lidR)
library(terra)
library(ggplot2)

### Path ####
#%%%%%%%%%%%%
data_raw = "C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-fullprocess/GER10_raw_sample.las"
data_understory="C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-fullprocess/GER10_sample_understory.las"

# 1) Read point cloud
raw_las <- readTLSLAS(data_raw)
stopifnot(!is.empty(las))
summary(understory_las)


# 2) denoise
# las_denoise= classify_noise(las,sor()) # removes to many points from upper canopy

ggplot(las@data, aes(x = Z)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  labs(
    x = "Z coordinate (height)",
    y = "Number of points",
    title = "Distribution of Z coordinates"
  ) +
  theme_minimal()+
  geom_vline(xintercept=c(quantile(las@data$Z,probs=0.001),quantile(las@data$Z,probs=0.999)))

z_min= -5
z_max= 40

las_filter <- filter_poi(
  las,
  Z > z_min & Z < z_max#,
  # Class = 7L
)

table(las_filter@data$Classification)

plot(las_filter,color="Classification")


# 2) Classify ground (choose one algorithm)
# CSF is popular and robust
las_grclass <- classify_ground(las_filter, csf(class_threshold = 0.2,cloth_resolution = 0.2))
table(las_grclass@data$Classification)

plot(las_grclass,color="Classification")

# 3) Build a DTM raster from ground points
dtm <- rasterize_terrain(las_grclass, res = 0.1, algorithm = knnidw())  # res in your CRS units (often meters)
plot(dtm)
# 4) Normalize heights using the DTM
las_n <- normalize_height(las_grclass, dtm)
plot(las_n)
# Optional: drop negative heights (can happen due to interpolation / noise)
las_n <- filter_poi(las_n, Z >= 0)

# 5) Save
writeLAS(las_n, "C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-fullprocess/GER10_normalized.laz")




tls = lidR::readTLSLAS(file)

# find tree positions as starting points for segmentation
map <- CspStandSegmentation::find_base_coordinates_raster(las_n)
# segment trees
segmented <- filter_poi(
  las_n,
  Classification == 0L
)
segmented<- segmented|>
  CspStandSegmentation::add_geometry(n_cores = parallel::detectCores()/2) |>
  CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = parallel::detectCores()/2)
writeLAS(segmented, "C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-fullprocess/GER10_segmented.laz")

# show results
lidR::plot(segmented, color = "TreeID")


### Inventory ####
understory_las <- readTLSLAS(data_understory)
map_init <- CspStandSegmentation::find_base_coordinates_raster(understory_las)

understory_segmented <- understory_las |>
  CspStandSegmentation::add_geometry(n_cores = 2) |>
  CspStandSegmentation::csp_cost_segmentation(map_init, 1, N_cores = 2)
plot(understory_segmented,color="TreeID")
# perform inventory
plot(understory_las,color="Original cloud index")
table(understory_las@data$`Original cloud index`)
inventory <- CspStandSegmentation::forest_inventory(understory_las,
                                                    slice_min = 0.2,
                                                    slice_max = 4,
                                                    increment = 0.2,
                                                    width = 0.1,
                                                    max_dbh = 1,
                                                    tree_id_col = "Original cloud index",
                                                    n_cores = 2)
