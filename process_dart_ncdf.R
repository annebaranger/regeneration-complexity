library(ncdf4)
library(dplyr)
library(terra)
plot_name="GER02"
user="amob2"
path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
path_folder=file.path(path_plot,
                      paste0(plot_name,"-processing"))
path_trees=file.path(path_folder,"trees")
dart_folder= paste0("C:/Users/",user,"/DART/user_data/simulations")


hours=c(8,10,12,14,16,18)


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
dim(arr4d)

y_idx <- seq_len(dim(arr4d)[1])
arr4d_flipped <- arr4d[rev(y_idx), , ,]


sum_arr <- apply(arr4d_flipped, c(1, 2, 3), sum)
sd_arr <- apply(arr4d_flipped, c(1, 2, 3), sd)


slice=apply(sd_arr, c(1, 3), sum)
plot(rast(sd_arr[,,50])/rast(sum_arr[,,50]))
plot(rast(sum_arr[,,10]))


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


#save(arr_norm,file=paste0(plot_name,"_rb.rdata"))
