library(ncdf4)
library(dplyr)
library(terra)
plot_name="GER01"
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
sum_arr <- apply(arr4d, c(1, 2, 3), sum)
sd_arr <- apply(arr4d, c(1, 2, 3), sd)


slice=apply(sd_arr, c(1, 3), sum)
plot(rast(sd_arr[,,50])/rast(sum_arr[,,50]))
plot(rast(sum_arr[,,50]))
