library(ncdf4)
library(dplyr)
library(terra)
plot_name="GER02"
user="Anne"
path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
path_folder=file.path(path_plot,
                      paste0(plot_name,"-processing"))
path_trees=file.path(path_folder,"trees")
dart_folder= paste0("C:/Users/",user,"/DART-1/user_data/simulations")


hours=c(8,10,12,14,16,18)



### Load netcdf ####
#%%%%%%%%%%%%%%%%%%%%%%%%
list_rb <- vector("list", length(hours))
names(list_rb) <- hours

for (i in seq_along(hours)) {
  h <- hours[i]
  
  path_ncdf <- file.path(
    dart_folder,
    paste0(plot_name, "_h", h),
    "output/netcdf/radiativeBudget/radiativeBudget_3D.nc"
  )
  
  nc <- nc_open(path_ncdf)
  
  list_rb[[i]] <- ncvar_get(nc, "3D/All/ITERX/Band_0/Total_Entry")
  
  nc_close(nc)
}
arr4d <- simplify2array(list_rb)
dims=dim(arr4d)

y_idx <- seq_len(dim(arr4d)[1])
arr4d_flipped <- arr4d[rev(y_idx), , ,]


sum_arr <- apply(arr4d_flipped, c(1, 2, 3), sum)
sd_arr <- apply(arr4d_flipped, c(1, 2, 3), sd)


# slice=apply(sd_arr, c(1, 3), sum)
# plot(rast(sd_arr[,,50])/rast(sum_arr[,,50]))
plot(rast(sum_arr[,,10]))

plot(c(rast(arr4d[,,30,1]),
       rast(arr4d[,,30,2]),
       rast(arr4d[,,30,3]),
       rast(arr4d[,,30,4]),
       rast(arr4d[,,30,5]),
       rast(arr4d[,,30,6])))

time_dense <- seq(min(hours), max(hours), by = 0.1)  # every 6 minutes

# Reshape to matrix [n_voxels × 6] for vectorised spline fitting
n_voxels <- prod(dims[1:3])
mat      <- matrix(arr4d, nrow = n_voxels, ncol = 6)  # each row = one voxel time series

# Fit a spline per voxel and evaluate on dense grid, then integrate
# (uses pracma::trapz or base diff/sum — no extra package needed)
daily_integral <- apply(mat, 1, function(y) {
  
  # Skip voxels that are all NA or all zero
  if (all(is.na(y)) || all(y == 0)) return(0)
  
  sp    <- splinefun(hours, y, method = "natural")
  y_hat <- sp(time_dense)
  y_hat[y_hat < 0] <- 0          # radiative values can't be negative
  
  # Trapezoid rule on dense grid (interval = 0.1 h)
  sum(diff(time_dense) * (head(y_hat, -1) + tail(y_hat, -1)) / 2)
})

daily_arr <- array(daily_integral, dim = dims[1:3])


plot(c(rast(daily_arr[,,8]),
     rast(sum_arr[,,8])))

### Normalize netcdf ####
#%%%%%%%%%%%%%%%%%%%%%%%%

dim(sum_arr)
nx <- dim(sum_arr)[2]
ny <- dim(sum_arr)[1]
nz <- dim(sum_arr)[3]

dtm_rast <- rast(file.path(path_folder,
                           list.files(path_folder,pattern="dtm.tif")))
# dtm_prj=project(dtm,crs(netcdf2rast))

xmin_vox <- ext(dtm_rast)[1]
ymin_vox <- ext(dtm_rast)[3]
zmin_vox <- min(dtm_rast)

# taille d'un voxel
dx <- 0.25
dy <- 0.25
dz <- 0.25

x_centers <- xmin_vox + (0:(nx - 1)) * dx + dx / 2
y_centers <- ymin_vox + (0:(ny - 1)) * dy + dy / 2

xy <- expand.grid(x = x_centers, y = y_centers)
# xy=cbind(xy,terra::extract(dtm_rast, xy)[, 2])
dtm_vals <- terra::extract(dtm_rast, xy)[, 2]
dtm_mat <- matrix(dtm_vals, nrow = nx, ncol = ny)

k_ground <- floor((dtm_mat - min(dtm_vals)) / dz) + 1

k_ground[k_ground < 1] <- 1
k_ground[k_ground > nz] <- nz

maxE=quantile(sum_arr,probs=0.99)[[1]]
arr_norm <- array(maxE, dim = dim(sum_arr))

for (i in 1:nx) {
  for (j in 1:ny) {
    kg <- k_ground[i, j]
    
    if (!is.na(kg) && kg < nz) {
      new_col <- sum_arr[j, i, (kg + 1):nz]
      arr_norm[j, i, 1:length(new_col)] <- new_col
    }
  }
}

plot(rast(sum_arr[6,,]))
plot(rast(arr_norm[6,,]))

plot(rast(arr_norm[,,1/dz]))


# save(arr_norm,file=paste0(plot_name,"_rb.rdata"))
