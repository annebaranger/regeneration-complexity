library(ncdf4)
library(dplyr)
library(terra)
plot_name="GER02"
hours=c(8,10,12,14,16,18)

dart_folder="C:/Users/Anne/DART-1/user_data/simulations"

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
netcdf2rast=rast(sum_arr)
dtm_rast <- rast(file.path(path_folder,
                           list.files(path_folder,pattern="dtm.tif")))
dtm_prj=project(dtm,crs(netcdf2rast))
