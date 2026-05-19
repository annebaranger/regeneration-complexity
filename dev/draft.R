rb=loadRData(tar_read(rb_GER01))[["12"]]
pc_norm=loadRData(tar_read(pc_norm_GER01))
rb_norm=loadRData(tar_read(rb_norm_GER01))[["12"]]
bbox=tar_read(bbox_GER01)
subplot_extent=tar_read(subplot_extent_GER01)
plot_name="GER01"
nstrat=4
tar_load(c("user","path_plot"))
vox_path=file.path(path_plot,paste0(plot_name,"-processing"),"vox_merge.txt")

library(readr)
df <- read_table(
  vox_path,
  col_names = TRUE
)

names(df) <- gsub('"', '', names(df))
dim(df)
# Optional: choose one z / k layer
df_k <- df %>%
  filter(j == 11)

k_lay=15
cowplot::plot_grid(
  # df %>% filter(k==k_lay) %>% ggplot(aes(x=i,y=j,fill=pad_transmittance))+geom_raster()+
  #   coord_equal()+
  #   scale_x_continuous(breaks = seq(0, 130, by = 10)) +
  #   scale_y_continuous(breaks = seq(0, 130, by = 10))+
  #   theme(legend.position = "none"),
  as.data.frame(rast(arr[,,k_lay]),xy=TRUE) %>% 
    ggplot(aes(x=x,y=y,fill=lyr.1))+geom_raster()+
    coord_equal()+
    scale_x_continuous(breaks = seq(0, 130, by = 10)) +
    scale_y_continuous(breaks = seq(0, 130, by = 10))+
    theme(legend.position = "none"),
  as.data.frame(rast(list_rb[['12']][,,k_lay]),xy=TRUE) %>% 
    ggplot(aes(x=x,y=y,fill=lyr.1))+geom_raster()+
    coord_equal()+
    scale_x_continuous(breaks = seq(0, 130, by = 10)) +
    scale_y_continuous(breaks = seq(0, 130, by = 10))+
    theme(legend.position = "none"),
  as.data.frame(rasterize_canopy(filter_poi(pc_norm, Z == k_lay/4), res = 0.1, algorithm = p2r()),xy=TRUE) %>%
    ggplot(aes(x=x,y=y,fill=Z))+geom_raster()+
    coord_equal()+
    scale_x_continuous(breaks = seq(0, 130, by = 10)) +
    scale_y_continuous(breaks = seq(0, 130, by = 10))+
    theme(legend.position = "none"),
  nrow = 1
)


# Turn into raster using x and y as coordinates
r <- rast(
  df_k[, c("i", "k", "pad_transmittance")],
  type = "xyz",
  crs = NA
)

plot(r)

x_vals <- sort(unique(df$x))
y_vals <- sort(unique(df$y))
z_vals <- sort(unique(df$z))

# Get all indices at once (vectorised)
ix <- match(df$x, x_vals)
iy <- match(df$y, y_vals)
iz <- match(df$z, z_vals)

# Create empty array
arr <- array(
  NA_real_,
  dim = c( length(y_vals), length(x_vals),length(z_vals)),
  dimnames = list( y = y_vals, x = x_vals, z = z_vals)
)

# Fill in one shot using a matrix of indices
arr[cbind( iy,ix, iz)] <- df$transmittance

# Dimensions
dim(arr)
dim(list_rb$`8`)
# Example access
arr[1, 1, 1]

# Slice at one z layer
plot(rast(arr[ , ,5 ]))

plot(rast(list_rb$`12`[ , ,5 ]))
plot(rast(tar_read(dtm_GER01)))



pc_norm=tar_read(pc_norm_GER01)
rb_norm=tar_read(rb_norm_GER01)
bbox=tar_read(bbox_GER01)
subplot_extent=tar_read(subplot_extent_GER01)
plot_name="GER01"
nstrat=4
pad=tar_read(voxnorm_GER01)

voxpad=loadRData(pad)



pad_slice=rast(apply(voxpad[,,(slice_low/0.25):(slice_up/0.25)], c(1, 2), sum))
ext(pad_slice)  <- ext(pc_slice_chm)
crs(pad_slice)  <- crs(pc_slice_chm)
pad_slice_matched <- resample(pad_slice, pc_slice_chm, method = "bilinear")
names(pad_slice_matched)="pad_sum"

pad_tot=sum(values(pad_slice_matched),na.rm=TRUE)
pad_mean=mean(values(pad_slice_matched),na.rm=TRUE)
pad_sd=sd(values(pad_slice_matched),na.rm=TRUE)
pad_cv=pad_sd/pad_mean