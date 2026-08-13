#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions.R  
#' @description R script containing all functions relative to data
#               processing
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#### ESSENTIAL FUNCTIONS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load Rdata file into a defined object
#' @param fileName n
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#### PLOT DESCRIPTORS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a bounding box of plot of interest, based on measures used in AMAPVox
#' @param plot_name
#' @param user user name for path 
get_bbox<-function(plot_name,
                   user){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  bbox_df=read.table(file.path(path_folder,
                               list.files(path_folder,pattern="bbox")),
                     sep="",header=F)
  bbox=list(xmin=bbox_df[1,2],xmax=bbox_df[1,3],
            ymin=bbox_df[2,2],ymax=bbox_df[2,3])
  return(bbox)
}

#' Load DTM and crop it to the accurate size
#' @param name description
get_dtm <-function(plot_name,
                   user,
                   bbox){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  
  # look for dtm file
  dtm_file=if_else(file.exists(file.path(path_folder,
                                         paste0(plot_name,"_dtm_bis.tif"))), # check if dtm was refined
                   file.path(path_folder,
                             paste0(plot_name,"_dtm_bis.tif")), # if so, chose refined file
                   file.path(path_folder,
                             paste0(plot_name,"_dtm.tif")) # if not, chose original dtm file
                   
  )
  # laod dtm file
  dtm_rast <- rast(dtm_file)
  
  dtm_rast<- crop(dtm_rast,
                  ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax),
                  ext=TRUE)
  
  dtm_path=paste0("output/",plot_name,"_dtm_cropped.tif")
  writeRaster(dtm_rast, filename = dtm_path, overwrite = TRUE)
  
  return(dtm_path)
}

#' Get the list of simulated times
#' @param dart_folder path of dart
#' @param plot_name plot name
get_simulation_time<-function(dart_folder,
                              plot_name){
  sim_tot=list.files(dart_folder)
  sim_plot=sim_tot[grepl(plot_name,sim_tot)]
  sim_plot_complete=c()
  for(s in sim_plot){
    if(file.exists(file.path(dart_folder,
                             s,
                             "output","netcdf","radiativeBudget","radiativeBudget_3D.nc"))){
      sim_plot_complete=c(sim_plot_complete,
                          substr(s,start = 8,stop = nchar(s)))
    }
  }
  sim_plot_complete=sort(as.numeric(sim_plot_complete))
  
  return(sim_plot_complete)
}

#' Get the extent of the subplots
#' @param plot_name plot name
#' @param user user
get_subplot_extent<-function(plot_name,
                             user){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  path_trees=file.path(path_folder,"trees")
  
  subplot_files=list.files(path_trees,pattern="subplot")
  
  list_hull <- vector("list", length(subplot_files))
  names(list_hull) <- lapply(subplot_files,function(x) substr(x,1,nchar(x)-4))
  
  for(f in  seq_along(subplot_files)){
    subplot=fread(file=file.path(path_trees,
                                 subplot_files[f]))[,1:3]
    colnames(subplot)=c("X","Y","Z")  
    subplot=LAS(subplot)
    
    list_hull[f]=st_convex_hull(subplot)
  }
  
  return(list_hull)
}

#### RADIATIVE BUDGET ####
#%%%%%%%%%%%%%%%%%%%%%%%%%

get_rb<-function(dart_folder,
                 plot_name,
                 hours){
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
    arr <- ncvar_get(nc, "3D/All/ITERX/Band_0/Total_Entry")
    nc_close(nc)
    
    y_idx <- seq_len(dim(arr)[1])
    list_rb[[i]] <- arr[rev(y_idx), , , drop = FALSE]
  }

  
  path_rb=paste0("output/",plot_name,"_rb_raw.rdata")
  save(list_rb,file=path_rb)
  
  return(path_rb)
}



normalize_rb<-function(plot_name,
                       user,
                       rb_path,
                       bbox,
                       dtm){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  
  list_rb=loadRData(rb_path)
  dtm_rast <- rast(dtm)
  

  for(r in seq_along(list_rb)){
    rb=list_rb[[r]]
    
    #get number of voxels
    nx <- dim(rb)[2]
    ny <- dim(rb)[1]
    nz <- dim(rb)[3]
    
    
    xmin_vox <- ext(dtm_rast)[1]
    ymin_vox <- ext(dtm_rast)[3]
    # zmin_vox <- min(values(dtm_rast))
    
    dx <- 0.25
    dy <- 0.25
    dz <- 0.25
    
    x_centers <- xmin_vox + (0:(nx - 1)) * dx + dx / 2
    y_centers <- ymin_vox + (0:(ny - 1)) * dy + dy / 2
    
    # check if dtm extent is too small, and replace highest value with max extent
    x_max <- ext(dtm_rast)[2][[1]]
    y_max <- ext(dtm_rast)[4][[1]]
    for (i in length(x_centers):1) {
      if (x_centers[i] > x_max) x_centers[i] <- x_max else break
    }
    for (i in length(y_centers):1) {
      if (y_centers[i] > y_max) y_centers[i] <- y_max else break
    }
    
    xy <- expand.grid(x = x_centers, y = y_centers)
    
    dtm_vals <- terra::extract(dtm_rast, xy)[, 2]
    dtm_mat <- matrix(dtm_vals, nrow = nx, ncol = ny)
    
    k_ground <- floor((dtm_mat - min(dtm_vals)) / dz) + 1
    
    k_ground[k_ground < 1] <- 1
    k_ground[k_ground > nz] <- nz
    
    maxE=quantile(rb,probs=0.99)[[1]]
    arr_norm <- array(maxE, dim = dim(rb))
    
    for (i in 1:nx) {
      for (j in 1:ny) {
        kg <- k_ground[i, j]
        
        if (!is.na(kg) && kg < nz) {
          new_col <- rb[j, i, (kg + 1):nz]
          arr_norm[j, i, 1:length(new_col)] <- new_col
        }
      }
    }
    list_rb[[r]]=arr_norm
  }
  
  path_rast_norm=paste0("output/",plot_name,"_norm.rdata")
  save(list_rb,file=path_rast_norm)
  return(path_rast_norm)
}


#### PLANT AREA DENSITY ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

get_pad<-function(user,
                  path_plot,
                  plot_name){
  vox_path=file.path(path_plot,paste0(plot_name,"-processing"),"vox_merge.txt")
  vox_df <- read_table(
    vox_path,
    col_names = TRUE
  )
  names(vox_df) <- gsub('"', '', names(vox_df))
  
  ### transform into array
  x_vals <- sort(unique(vox_df$x))
  y_vals <- sort(unique(vox_df$y))
  z_vals <- sort(unique(vox_df$z))
  
  # Get all indices at once (vectorised)
  ix <- match(vox_df$x, x_vals)
  iy <- match(vox_df$y, y_vals)
  iy_flipped <- length(y_vals) + 1 - iy
  iz <- match(vox_df$z, z_vals)
  
  # Create empty array
  arr <- array(
    NA_real_,
    dim = c(length(y_vals), length(x_vals), length(z_vals)),
    # dimnames = list( y = y_vals, x = x_vals, z = z_vals)
  )
  
  # Fill in one shot using a matrix of indices
  arr[cbind(iy_flipped, ix, iz)] <- vox_df$pad_transmittance
  
  ### normalize
  arr_norm <- apply(arr, c(1, 2), function(col) {
    first_valid <- which(!is.na(col))[1]
    if (is.na(first_valid)) return(rep(NA_real_, length(col)))  # all NA column
    new_col <- col[first_valid:length(col)]
    c(new_col, rep(NA_real_, length(col) - length(new_col)))    # pad end with NA
  })
  arr_norm <- aperm(arr_norm, c(2, 3, 1))
  
  path_vox=paste0("output/",plot_name,"_voxnorm.rdata")
  save(arr_norm,file=path_vox)
  return(path_vox)
}



#### POINT CLOUD ####
#%%%%%%%%%%%%%%%%%%%%


get_pc_norm_trees<-function(plot_name,
                      user,
                      subfolders=c("adults","regen","hesitation","dead-incomplete","tree-parts")
){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  path_trees=file.path(path_folder,"trees")

  ### Create one point cloud 
  #%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
  ### DTM LAS 
  #%%%%%%%%%%%
  dtm_rast <- rast(file.path(path_folder,
                             list.files(path_folder,pattern="dtm.tif")))
  las <- clip_rectangle(las,
                        xleft   = ext(dtm_rast)$xmin,
                        ybottom = ext(dtm_rast)$ymin,
                        xright  = ext(dtm_rast)$xmax,
                        ytop    = ext(dtm_rast)$ymax)
  las_norm <- normalize_height(las, dtm_rast)
  
  las_norm@data$X <- las_norm@data$X - min(las_norm@data$X)
  las_norm@data$Y <- las_norm@data$Y - min(las_norm@data$Y)
  
  las_norm <- las_update(las_norm)
  
  path_las=paste0("output/",plot_name,"_pcnorm.rdata")
  save(las_norm,file=path_las)
  return(path_las)
}


get_pc_norm_clean<-function(plot_name,
                            user,
                            bbox,
                            dtm){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))

  # load cleaned and denoised point cloud
  f=file.path(path_folder,
              list.files(path_folder,pattern="clean.las"))
  pc=readLAS(f)
  las <- clip_rectangle(pc,
                        xleft   = bbox$xmin,
                        ybottom = bbox$ymin,
                        xright  = bbox$xmax,
                        ytop    = bbox$ymax)
  
  # load dtm 
  dtm_rast <- rast(dtm)
  
  # normalize pc height
  las_norm <- normalize_height(las, dtm_rast)
  
  las_norm@data$X <- las_norm@data$X - min(las_norm@data$X)
  las_norm@data$Y <- las_norm@data$Y - min(las_norm@data$Y)
  
  las_norm <- las_update(las_norm)
  
  las_clip <- filter_poi(las_norm, Z < 5)
  
  path_las=paste0("output/",plot_name,"_pcnorm.rdata")
  save(las_clip,file=path_las)
  return(path_las)
}

get_gap_filling<-function(pc,voxel_size=0.25){
  pc_vox=voxelize_points(pc,res=0.25)
  occupied_voxels <- nrow(pc_vox)
  
  xrange <- range(pc@data$X)
  yrange <- range(pc@data$Y)
  zrange <- range(pc@data$Z)
  
  nx <- ceiling(diff(xrange) / voxel_size)
  ny <- ceiling(diff(yrange) / voxel_size)
  nz <- ceiling(diff(zrange) / voxel_size)
  
  total_voxels <- nx * ny * nz
  empty_voxels <- total_voxels - occupied_voxels
  
  empty_prop <- empty_voxels / total_voxels
  occupied_prop <- occupied_voxels / total_voxels
  return(list(empty=empty_prop,occupied=occupied_prop))

}


get_gini<-function(x){
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (all(x == 0)) return(0)
  x <- sort(x)
  n <- length(x)
  (2 * sum(seq_len(n) * x)) / (n * sum(x)) - (n + 1) / n
}
# pc_norm=tar_read(pc_norm_GER13)
# rb_norm=tar_read(rb_norm_GER13)
# voxnorm=tar_read(voxnorm_GER13)
# bbox=tar_read(bbox_GER13)
# subplot_extent=tar_read(subplot_extent_GER13)
# plot_name="GER13"
# nstrat=4


# Compute every per-slice metric from the already-prepared rasters.
# `chm`      : canopy height raster (layer named "Z"), NA-excluded stats
# `chm_num`  : same but NA -> 0
# `rb`       : resampled, quantile-filtered relative-brightness raster (named "rb")
# `rb_mask`  : rb masked by chm
# `pad`      : resampled plant-area-density raster
# `pc_gap`   : point cloud (slice or clipped slice) for gap filling
.slice_metrics <- function(chm, chm_num, rb, rb_mask, pad, pad_cube,
                           pc_gap, slice_up,
                           voxel_size = 0.25, min_pts_lm = 3) {
  
  v_chm <- values(chm)
  h_mean_na <- mean(v_chm, na.rm = TRUE)
  h_med_na  <- median(v_chm, na.rm = TRUE)
  h_sd_na   <- sd(v_chm, na.rm = TRUE)
  
  v_chm_num <- values(chm_num)
  h_mean <- mean(v_chm_num, na.rm = TRUE)
  h_med  <- median(v_chm_num, na.rm = TRUE)
  h_sd   <- sd(v_chm_num, na.rm = TRUE)
  
  v_rb <- values(rb)
  rb_mean <- mean(v_rb, na.rm = TRUE)
  rb_med  <- median(v_rb, na.rm = TRUE)
  rb_sd   <- sd(v_rb, na.rm = TRUE)
  
  v_rbm <- values(rb_mask)
  rb_mean_mask <- mean(v_rbm, na.rm = TRUE)
  rb_med_mask  <- median(v_rbm, na.rm = TRUE)
  rb_sd_mask   <- sd(v_rbm, na.rm = TRUE)
  
  # linear trend rb ~ Z over vegetated pixels below the slice top
  df <- as.data.frame(c(chm, rb), xy = TRUE, na.rm = TRUE) %>%
    filter(Z < (slice_up - 0.01))
  coef <- if (nrow(df) > min_pts_lm) {
    lm(rb ~ Z, data = df)$coefficients[["Z"]]
  } else {
    NA_real_
  }
  
  gf <- get_gap_filling(pc_gap, voxel_size = voxel_size)
  
  v_pad <- values(pad)
  pad_tot  <- sum(v_pad, na.rm = TRUE)
  pad_mean <- mean(v_pad, na.rm = TRUE)
  pad_sd   <- sd(v_pad, na.rm = TRUE)
  
  profile <- as.vector(global(pad_cube, "sum", na.rm = TRUE)[, 1])
  fhd     <- .fhd(profile)
  
  
  ri_na <- rumple_index(chm)
  ri <- rumple_index(chm_num)
  
  list(
    h_mean_na = h_mean_na, h_med_na = h_med_na, h_sd_na = h_sd_na,
    h_cv_na = h_sd_na / h_mean_na,
    h_mean = h_mean, h_med = h_med, h_sd = h_sd, h_cv = h_sd / h_mean,
    empty = gf$empty, occupied = gf$occupied,
    rb_mean = rb_mean, rb_med = rb_med, rb_sd = rb_sd, rb_cv = rb_sd / rb_mean,
    rb_mean_mask = rb_mean_mask, rb_med_mask = rb_med_mask, rb_sd_mask = rb_sd_mask,
    rb_cv_mask = rb_sd_mask / rb_mean_mask,
    coef = coef,
    pad_tot = pad_tot, pad_mean = pad_mean, pad_sd = pad_sd,
    pad_cv = pad_sd / pad_mean,
    ri=ri,
    ri_na=ri_na,
    fhd=fhd
  )
}

# Assemble one output row in the exact column order of the original data frame.
.make_row <- function(plot_name, subplot, time, slice_num, slice_low, slice_up, m) {
  data.frame(
    plot = plot_name, subplot = subplot, time = time,
    slice_num = slice_num, slice_low = slice_low, slice_up = slice_up,
    h_mean_na = m$h_mean_na, h_med_na = m$h_med_na, h_sd_na = m$h_sd_na, h_cv_na = m$h_cv_na,
    h_mean = m$h_mean, h_med = m$h_med, h_sd = m$h_sd, h_cv = m$h_cv,
    empty_gf = m$empty, occupied_gf = m$occupied,
    rb_mean = m$rb_mean, rb_med = m$rb_med, rb_sd = m$rb_sd, rb_cv = m$rb_cv,
    rb_mean_mask = m$rb_mean_mask, rb_med_mask = m$rb_med_mask,
    rb_sd_mask = m$rb_sd_mask, rb_cv_mask = m$rb_cv_mask,
    coef_lin = m$coef,
    pad_tot = m$pad_tot, pad_mean = m$pad_mean, pad_sd = m$pad_sd, pad_cv = m$pad_cv,
    ri=m$ri,ri_na=m$ri_na,fhd=m$fhd,
    stringsAsFactors = FALSE
  )
}



# Shannon diversity of a vertical PAD profile.
#   profile : numeric vector, one value per height layer
# Zero / non-finite layers are dropped (0*log(0) is NaN, log(0) is -Inf).
# Returns NA if fewer than 2 occupied layers, where FHD is undefined/degenerate.
.fhd <- function(profile) {
  p <- profile[is.finite(profile) & profile > 0]
  if (length(p) < 2L) return(NA_real_)
  p <- p / sum(p)
  -sum(p * log(p))
}



get_complexity <- function(pc_norm,
                           rb_norm,
                           voxnorm,
                           bbox,
                           subplot_extent,
                           plot_name,
                           sliced    = TRUE,   # TRUE = strata, FALSE = single band
                           nstrat    = 4,      # used only when sliced = TRUE
                           slice_low = 0,      # used only when sliced = FALSE
                           slice_up  = 4) {    # used only when sliced = FALSE
  
  pc_norm <- loadRData(pc_norm)
  list_rb <- loadRData(rb_norm)
  voxpad  <- loadRData(voxnorm)
  
  # --- define the height bands to process (the only thing the toggle changes) ---
  if (sliced) {
    edges  <- seq(from = 0.5, by = 1, length.out = nstrat + 1)
    slices <- data.frame(low = edges[-(nstrat + 1)],
                         up  = edges[-1],
                         num = seq_len(nstrat))
  } else {
    slices <- data.frame(low = slice_low, up = slice_up, num = 0)
  }
  
  # subplot geometries are constant -> compute once
  sub_geo <- lapply(subplot_extent, function(pe) {
    shift(vect(pe), dx = -bbox$xmin, dy = -bbox$ymin)
  })
  
  n_rows <- length(list_rb) * nrow(slices) * (1L + length(subplot_extent))
  rows <- vector("list", n_rows)
  k <- 0L
  
  for (t in seq_along(list_rb)) {
    rb_arr <- list_rb[[t]]
    time   <- as.numeric(names(list_rb)[t])
    
    for (s in seq_len(nrow(slices))) {
      s_low <- slices$low[s]
      s_up  <- slices$up[s]
      s_num <- slices$num[s]
      z_idx <- (s_low / 0.25):(s_up / 0.25)
      
      ## ---- height CHM ----
      pc_slice <- filter_poi(pc_norm, Z >= s_low & Z < s_up)
      chm <- rasterize_canopy(pc_slice, res = 0.1, algorithm = p2r())
      chm[chm > (s_up - 0.01)] <- NA
      
      chm_num <- chm
      chm_num[is.na(chm_num)] <- 0
      
      ## ---- relative brightness ----
      rb_slice <- rast(rowSums(rb_arr[, , z_idx, drop = FALSE], dims = 2))
      ext(rb_slice) <- ext(chm)
      crs(rb_slice) <- crs(chm)
      rb_matched <- resample(rb_slice, chm, method = "bilinear")
      names(rb_matched) <- "rb"
      thr <- quantile(values(rb_matched), na.rm = TRUE, probs = 0.005)
      rb_matched[rb_matched < thr] <- NA
      rb_matched_mask <- mask(rb_matched, chm)
      
      
      ## ---- plant area density ----
      pad_slice <- rast(rowSums(voxpad[, , z_idx, drop = FALSE], dims = 2))
      ext(pad_slice) <- ext(chm)
      crs(pad_slice) <- crs(chm)
      pad_matched <- resample(pad_slice, chm, method = "bilinear")
      names(pad_matched) <- "pad_sum"
      
      pad_cube <- rast(voxpad[, , z_idx, drop = FALSE])
      ext(pad_cube) <- ext(chm)
      crs(pad_cube) <- crs(chm)
      pad_cube_matched <- resample(pad_cube, chm, method = "bilinear")
      
      
      ## ---- plot-scale metrics ----
      m <- .slice_metrics(chm, chm_num, rb_matched, rb_matched_mask,
                          pad_matched, pad_cube, pc_slice, s_up)
      k <- k + 1L
      rows[[k]] <- .make_row(plot_name, "plot", time, s_num, s_low, s_up, m)
      
      ## ---- subplot-scale metrics ----
      for (sp in seq_along(subplot_extent)) {
        geo    <- sub_geo[[sp]]
        chm_s  <- crop(chm,             geo)
        chmn_s <- crop(chm_num,         geo)
        rb_s   <- crop(rb_matched,      geo)
        rbm_s  <- crop(rb_matched_mask, geo)
        pad_s  <- crop(pad_matched,     geo)
        cube_s <- crop(pad_cube_matched, geo)
        pc_s   <- clip_roi(pc_slice, sf::st_as_sf(geo))
        
        m_s <- .slice_metrics(chm_s, chmn_s, rb_s, rbm_s, pad_s, cube_s,
                              pc_s, s_up)
        k <- k + 1L
        rows[[k]] <- .make_row(plot_name, names(subplot_extent[sp]),
                               time, s_num, s_low, s_up, m_s)
      }
    }
  }
  
  dplyr::bind_rows(rows)
}


#' # get extinction characteristics
#' #' @param plot_name plot name
#' #' @param targets characteristic height to compute
#' get_Hext<-function(rb_norm,
#'                    pc_norm,
#'                    subplot_extent,
#'                    bbox,
#'                    plot_name,
#'                    targets=c(0.5,0.88)){
#'   list_rb<-loadRData(rb_norm)         # load radiative budget
#'   pc_norm<-loadRData(pc_norm)
#'   
#'   summary_height=data.frame(plot=character(),
#'                             subplot=character(),
#'                             time=numeric(),
#'                             extinction=numeric(),
#'                             mean=numeric(),
#'                             sd=numeric(),
#'                             mean_slope=numeric(),
#'                             sd_slope=numeric())
#'   
#'   
#'   for(t in seq_along(list_rb)){      # loop over simulation time
#'     rb_norm=list_rb[[t]]
#'     time=as.numeric(names(list_rb)[t])
#'     
#'     
#'     prof <- reshape2::melt(rb_norm, 
#'                            varnames = c("i", "j", "z"), 
#'                            value.name = "value") %>%
#'       group_by(i, j) %>%
#'       mutate(extinction = 1 - value / max(value),
#'              z   = z * 0.25)
#'     
#'     # z prédit pour chaque cible de ext, par courbe (i,j)
#'     half <- prof %>%
#'       group_modify(~{
#'         m   <- mgcv::gam(z ~ s(extinction, k = 5), data = .x)
#'         eps <- 1e-2
#'         z_pred <- as.numeric(predict(m, newdata = tibble(extinction = targets)))
#'         z_hi   <- as.numeric(predict(m, newdata = tibble(extinction = targets + eps)))
#'         z_lo   <- as.numeric(predict(m, newdata = tibble(extinction = targets - eps)))
#'         tibble(extinction = targets,
#'                z_pred = z_pred,
#'                slope  = (z_hi - z_lo) / (2 * eps))   # dz/d(extinction)
#'       }) %>%
#'       ungroup()
#'     
#'     summary_height <- summary_height %>%
#'       bind_rows(half %>%
#'                   group_by(extinction) %>%
#'                   summarise(mean       = mean(z_pred),
#'                             sd         = sd(z_pred),
#'                             mean_slope = mean(slope),
#'                             sd_slope   = sd(slope),
#'                             .groups = "drop") %>%
#'                   mutate(time = time, plot = plot_name,subplot="plot"))
#'     
#'     # get a raster 
#'     wide <- half %>%
#'       pivot_wider(id_cols = c(i, j),
#'                   names_from  = extinction,
#'                   values_from = c(z_pred, slope))
#'     
#'     stk <- terra::rast(wide, type = "xyz")
#'     names(stk) 
#'     
#'     # prepare reference raster
#'     pc_chm <- rasterize_canopy(pc_norm, res = 0.15, algorithm = p2r())
#'     ext(stk)  <- ext(pc_chm)
#'     crs(stk)  <- crs(pc_chm)
#'     stk_matched <- resample(stk, pc_chm, method = "bilinear")
#'     
#'     for(sp in seq_along(subplot_extent)){
#'       plot_ext=subplot_extent[[sp]]
#'       plot_ext_trans <- shift(vect(plot_ext), 
#'                               dx = -bbox$xmin, 
#'                               dy = -bbox$ymin)
#'       
#'       stk_sub=crop(stk_matched,plot_ext_trans)
#'       
#'       
#'       summary_height <- summary_height %>%
#'         bind_rows(as.data.frame(stk_sub,xy=TRUE) %>% 
#'                     pivot_longer(
#'                       cols = -c(x, y),
#'                       names_to  = c(".value", "extinction"),
#'                       names_sep = "_(?=[0-9])"          # split at the underscore before the number
#'                     ) %>% 
#'                     group_by(extinction) %>%
#'                     summarise(mean       = mean(z_pred),
#'                               sd         = sd(z_pred),
#'                               mean_slope = mean(slope),
#'                               sd_slope   = sd(slope),
#'                               .groups = "drop") %>%
#'                     mutate(extinction=as.numeric(extinction),
#'                            time = time,
#'                            plot = plot_name,
#'                            subplot=names(subplot_extent[sp])))
#'       
#'     }
#'   }
#'   return(summary_height)
#' }


## linear interpolation of a profile at arbitrary heights
.interp <- function(z, v, zout) stats::approx(z, v, xout = zout, rule = 2)$y

## height at which an extinction profile first reaches `target`, coming down
## from the top of the canopy. `z` must increase with the index.
.z_at_target <- function(z, e, target) {
  k <- which(e >= target)
  if (!length(k)&e[1]>=target) return(NA_real_)
  if (!length(k)&e[1]<target) return(0)
  k <- max(k)                                   # highest voxel still >= target
  if (k == length(e)) return(z[k])              # target reached at the very top
  z[k] + (target - e[k]) * (z[k + 1] - z[k]) / (e[k + 1] - e[k])
}

## all the metrics of ONE vertical profile (a column, or an aggregated unit)
.profile_metrics <- function(z, e, rb, targets, z_ref, dz, smooth = FALSE) {
  
  ok <- is.finite(e)
  if (sum(ok) < 5) return(NULL)
  zo <- z[ok]; eo <- e[ok]
  
  f <- if (smooth) {
    fit <- try(mgcv::gam(e ~ s(z, k = 5), data = data.frame(z = zo, e = eo)),
               silent = TRUE)
    if (inherits(fit, "try-error")) function(x) .interp(zo, eo, x)
    else function(x) as.numeric(predict(fit, newdata = data.frame(z = x)))
  } else {
    function(x) .interp(zo, eo, x)
  }
  
  ## central finite difference -> d(extinction)/dz, in m-1 (negative)
  slope <- function(z0) {
    if(z0!=0) {
      (f(z0 + dz) - f(z0 - dz)) / (2 * dz)
      } else (f(2*dz) - f(dz)) / dz
  } 
  
  zc <- vapply(targets, function(tg) .z_at_target(zo, eo, tg), numeric(1))
  
  tibble(
    target      = targets,
    z_cross     = zc,
    slope_cross = vapply(zc, function(x) if (is.na(x)) NA_real_ else slope(x),
                         numeric(1)),
    ext_ref     = .interp(zo, eo, z_ref),       # raw value, never smoothed
    slope_ref   = slope(z_ref),
    rb_ref      = .interp(z, rb, z_ref)
  )
}

## PAD metrics of ONE vertical PAI profile
.pad_metrics <- function(z_top, pai, z_ref, voxel_size) {
  pai <- ifelse(is.na(pai), 0, pai)
  tot <- sum(pai)
  w   <- pmin(pmax((z_top - z_ref) / voxel_size, 0), 1)   # fraction above z_ref
  z50 <- NA_real_
  if (tot > 0) {
    cs   <- cumsum(pai)
    k    <- which(cs >= 0.5 * tot)[1]
    prev <- if (k == 1L) 0 else cs[k - 1L]
    z50  <- (z_top[k] - voxel_size) + (0.5 * tot - prev) / pai[k] * voxel_size
  }
  c(padAbove = sum(pai * w), padTot = tot, zPAD50 = z50)
}

## sum an [i, j, k] array over a set of columns -> one 1-D profile
.profile_sum <- function(arr, ij) {
  m <- cbind(ij$i, ij$j)
  vapply(seq_len(dim(arr)[3]),
         function(k) sum(arr[cbind(m, k)], na.rm = TRUE), numeric(1))
}

## mean / sd of a set of columns, for every unit
.summarise_columns <- function(d, vars, groups) {
  d %>%
    group_by(across(all_of(groups))) %>%
    summarise(across(all_of(vars),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd   = ~sd(.x,   na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              n_col = n(), .groups = "drop") %>%
    .add_cv()
}

.add_cv <- function(d) {
  vars <- sub("_mean$", "", grep("_mean$", names(d), value = TRUE))
  for (v in vars) d[[paste0(v, "_cv")]] <- d[[paste0(v, "_sd")]] / d[[paste0(v, "_mean")]]
  d
}

## turn an (i, j, ...) table back into a georeferenced raster stack
.ij_rast <- function(d, tmpl) {
  ni <- nrow(tmpl); nj <- ncol(tmpl)
  vars <- setdiff(names(d), c("i", "j"))
  out <- lapply(vars, function(v) {
    m <- matrix(NA_real_, ni, nj)
    m[cbind(d$i, d$j)] <- d[[v]]
    terra::setValues(tmpl, as.vector(t(m[ni:1, , drop = FALSE])))  # row 1 = north
  })
  out <- terra::rast(out)
  names(out) <- vars
  out
}


#' Extinction and PAD characteristics of a plot and its subplots
#'
#' @param rb_norm path to the voxelised, height-normalised radiative budget
#'   (a named list of (i, j, z) arrays, names = simulation time)
#' @param pc_norm path to the height-normalised point cloud (used only for the
#'   extent / crs of the voxel grid)
#' @param voxnorm path to the voxelised, height-normalised PAD array (i, j, z),
#'   same horizontal grid as rb_norm
#' @param subplot_extent named list of subplot polygons (sf / SpatVector), in
#'   the original coordinate system
#' @param bbox bounding box used to normalise the point cloud (list with
#'   $xmin, $ymin)
#' @param plot_name plot name
#' @param targets characteristic extinction levels at which the height of the
#'   profile is computed
#' @param z_ref reference height (m) at which extinction, its vertical slope and
#'   the PAD above are evaluated
#' @param voxel_size vertical voxel size (m); also the finite-difference step
#' @param pad_unit "density" if voxnorm stores PAD in m2/m3 (values are
#'   multiplied by voxel_size to get a plant area index), "index" if each voxel
#'   already holds m2/m2
#' @param smooth_columns fit a gam(extinction ~ s(z, k = 5)) before taking the
#'   per-column slopes (your original behaviour). Costs one gam per column and
#'   per time step; set to FALSE for raw finite differences.
#' @param touches passed to terra::extract(): TRUE also keeps the cells merely
#'   touched by a subplot polygon, FALSE (default) keeps only those whose
#'   centre falls inside
#' @param return_rasters also return the per-column metrics as raster stacks
#'
#' @return a tibble with one row per (subplot, time, target), or a list if
#'   return_rasters = TRUE
get_Hext <- function(rb_norm,
                     pc_norm,
                     voxnorm,
                     subplot_extent,
                     bbox,
                     plot_name,
                     targets        = c(0.5, 0.9),
                     z_ref          = 4.5,
                     voxel_size     = 0.25,
                     pad_unit       = c("density", "index"),
                     smooth_columns = TRUE,
                     touches        = FALSE,
                     return_rasters = FALSE) {
  
  list_rb <- loadRData(rb_norm)
  pc      <- loadRData(pc_norm)
  voxpad  <- loadRData(voxnorm)
  
  if (!identical(dim(voxpad)[1:2], dim(list_rb[[1]])[1:2])){
    if(abs(dim(voxpad)[1]-dim(list_rb[[1]])[1])<=2 &
       abs(dim(voxpad)[2]-dim(list_rb[[1]])[2])<=2){
      dim_x=min(dim(voxpad)[1],dim(list_rb[[1]])[1])
      dim_y=min(dim(voxpad)[2],dim(list_rb[[1]])[2])
      voxpad=voxpad[1:dim_x,1:dim_y,]
      for(l in seq_along(list_rb)){
        list_rb[[l]]=list_rb[[l]][1:dim_x,1:dim_y,]
        }
    }else{
      stop("voxpad and rb_norm do not share the same horizontal grid")
    }
    }  
  ni <- dim(voxpad)[1]
  nj <- dim(voxpad)[2]
  nz <- dim(voxpad)[3]
  
  # -------------------------------------------------------------------
  # 1. geolocate the voxel grid and find which columns belong to which unit
  #    (i indexes rows south -> north, j indexes columns west -> east:
  #     this is the orientation implied by your rast(wide, type = "xyz"))
  # -------------------------------------------------------------------
  chm  <- lidR::rasterize_canopy(pc, res = 0.15, algorithm = lidR::p2r())
  tmpl <- terra::rast(nrows = ni, ncols = nj,
                      extent = terra::ext(chm), crs = terra::crs(chm))
  
  r_ij <- c(terra::setValues(tmpl, rep(ni:1,        each  = nj)),   # i
            terra::setValues(tmpl, rep(seq_len(nj), times = ni)))   # j
  names(r_ij) <- c("i", "j")
  
  sub_lut <- purrr::imap_dfr(subplot_extent, function(e, nm) {
    if (inherits(e, "bbox")) e <- sf::st_as_sfc(e)
    p <- terra::shift(terra::vect(e), dx = -bbox$xmin, dy = -bbox$ymin)
    v <- terra::extract(r_ij, p, touches = touches)
    if (!nrow(v)) warning("no voxel column falls inside subplot ", nm)
    tibble(subplot = nm, i = v$i, j = v$j)
  }) %>% filter(!is.na(i), !is.na(j))
  
  ## the plot itself = every column, handled exactly like a subplot
  sub_lut <- bind_rows(
    expand_grid(subplot = "plot", i = seq_len(ni), j = seq_len(nj)),
    sub_lut
  )
  
  # -------------------------------------------------------------------
  # 2. PAD metrics (time invariant)
  # -------------------------------------------------------------------
  pai   <- voxpad * if (pad_unit == "density") voxel_size else 1
  z_pad <- seq_len(nz) * voxel_size                       # top of each voxel
  
  pad_arr <- apply(pai, 1:2, .pad_metrics,
                   z_top = z_pad, z_ref = z_ref, voxel_size = voxel_size)
  
  pad_col <- tibble(i        = rep(seq_len(ni), times = nj),
                    j        = rep(seq_len(nj), each  = ni),
                    padAbove = as.vector(pad_arr["padAbove", , ]),
                    padTot   = as.vector(pad_arr["padTot",   , ]),
                    zPAD50   = as.vector(pad_arr["zPAD50",   , ]))
  
  ## distribution of the per-column values, per unit
  pad_sum <- sub_lut %>%
    left_join(pad_col, by = c("i", "j")) %>%
    .summarise_columns(vars   = c("padAbove", "padTot", "zPAD50"),
                       groups = "subplot")
  
  ## same metrics on the aggregated profile of each unit
  pad_agg <- sub_lut %>%
    group_by(subplot) %>%
    group_modify(~ tibble::as_tibble_row(
      .pad_metrics(z_pad, .profile_sum(pai, .x), z_ref, voxel_size))) %>%
    ungroup() %>%
    rename_with(~paste0(.x, "_agg"), c(padAbove, padTot, zPAD50))
  
  # -------------------------------------------------------------------
  # 3. light metrics, one pass per simulation time
  # -------------------------------------------------------------------
  grid  <- expand_grid(i = seq_len(ni), j = seq_len(nj))
  l_ras <- list()
  
  light <- imap_dfr(list_rb, function(rb, nm) {
    
    time   <- as.numeric(nm)
    nz_rb  <- dim(rb)[3]
    z_rb   <- seq_len(nz_rb) * voxel_size
    rb_max <- max(rb, na.rm = TRUE)          # brightest voxel of the scene
    
    ## --- per column -------------------------------------------------
    col_met <- map2_dfr(grid$i, grid$j, function(ii, jj) {
      rb_prof <- rb[ii, jj, ]
      if (all(!is.finite(rb_prof))) return(NULL)
      m <- .profile_metrics(z = z_rb, e = 1 - rb_prof / rb_max, rb = rb_prof,
                            targets = targets, z_ref = z_ref, dz = voxel_size,
                            smooth = smooth_columns)
      if (is.null(m)) return(NULL)
      mutate(m, i = ii, j = jj, .before = 1)
    })
    
    if (return_rasters)
      l_ras[[nm]] <<- .ij_rast(
        col_met %>%
          pivot_wider(id_cols     = c(i, j),
                      names_from  = target,
                      values_from = c(z_cross, slope_cross)) %>%
          left_join(distinct(col_met, i, j, ext_ref, slope_ref, rb_ref),
                    by = c("i", "j")),
        tmpl)
    
    col_sum <- sub_lut %>%
      left_join(col_met, by = c("i", "j"), relationship = "many-to-many") %>%
      filter(!is.na(target)) %>%
      .summarise_columns(vars   = c("z_cross", "slope_cross",
                                    "ext_ref", "slope_ref", "rb_ref"),
                         groups = c("subplot", "target"))
    
    ## --- aggregated profile of each unit ----------------------------
    agg <- sub_lut %>%
      group_by(subplot) %>%
      group_modify(~{
        rb_prof <- .profile_sum(rb, .x)
        .profile_metrics(z = z_rb, e = 1 - rb_prof / max(rb_prof, na.rm = TRUE),
                         rb = rb_prof, targets = targets, z_ref = z_ref,
                         dz = voxel_size, smooth = FALSE)
      }) %>%
      ungroup() %>%
      rename_with(~paste0(.x, "_agg"),
                  c(z_cross, slope_cross, ext_ref, slope_ref, rb_ref))
    
    col_sum %>%
      left_join(agg, by = c("subplot", "target")) %>%
      mutate(time = time, .before = 1)
  })
  
  # -------------------------------------------------------------------
  # 4. assemble
  # -------------------------------------------------------------------
  out <- light %>%
    left_join(pad_sum %>% select(-n_col), by = "subplot") %>%
    left_join(pad_agg,                    by = "subplot") %>%
    mutate(plot = plot_name, .before = 1) %>%
    arrange(subplot, time, target)
  
  if (!return_rasters) return(out)
  
  list(summary     = out)#,
  # pad_rast    = .ij_rast(pad_col, tmpl),
  # light_rast  = l_ras,
  # subplot_lut = sub_lut)
}



### REGENERATION METRICS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get regeneration metrics at the subplot level
#' @param regen_df table of the data from the field
get_regen_metric_subplot<-function(regen_df,shade_df){
  ### ALL SPECIES TOGETHER 
  #%%%%%%%%%%%%%%%%%%%%%%%
  
  # Abundance & richness
  general_df <- regen_df %>%
    group_by(plot, subplot) %>%
    summarise(
      richness     = n_distinct(Species),
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    )
  
  var_df <- regen_df %>%
    group_by(plot, subplot, Quadrat) %>%
    summarise(
      richness  = n_distinct(Species),
      abundance = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(plot, subplot) %>%
    summarise(
      richness_sd  = sd(richness),
      abundance_sd = sd(abundance),
      abundance_cv = sd(abundance) / mean(abundance),
      .groups = "drop"
    )
  
  # Abundance & richness - seedlings and saplings
  seedling_df <- regen_df %>%
    filter(!is.na(Species), Hclass %in% c(1, 2)) %>%
    group_by(plot, subplot) %>%
    summarise(richness_seedling     = n_distinct(Species),
              n_seedling = sum(Count), .groups = "drop")
  
  sapling_df <- regen_df %>%
    filter(!is.na(Species), !Hclass %in% c(1, 2)) %>%
    group_by(plot, subplot) %>%
    summarise(richness_sapling   = n_distinct(Species),
              n_sapling = sum(Count), .groups = "drop")
  
  # Shannon index total and per size class
  shannon_df <- regen_df %>%
    filter(!is.na(as.numeric(Count))) %>%
    group_by(plot, subplot, Species) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, subplot) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H = -sum(p), .groups = "drop")
  
  shannon_size_df <- regen_df %>%
    filter(!is.na(as.numeric(Count))) %>%
    group_by(plot, subplot, Species, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, subplot, Hclass) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H = -sum(p), .groups = "drop") %>%
    pivot_wider(names_from = Hclass, values_from = H, names_prefix = "H_hclass")
  
  # Growth rate, total and per size class
  growth_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, subplot) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    )
  
  growth_size_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, subplot, Hclass) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = Hclass,
      values_from = c(Height_increment_mn, Height_increment_cv),
      names_sep   = "_hclass"
    )
  
  # Height distribution
  shannon_height_df <-  regen_df %>%
    filter(!is.na(Hclass)) %>% 
    mutate(Count=replace_na(Count,0)) %>% 
    group_by(plot, subplot, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, subplot) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H_height = -sum(p), .groups = "drop")
  
  # browsing intensity
  browsing_df <- regen_df %>%
    filter(Hclass > 1) %>%
    group_by(plot, subplot) %>%
    summarise(browsing = sum(Browsing == "Y") / n(), .groups = "drop")
  
  # community wighted mean for shade tolerange
  shade_tol_df<- regen_df %>%
    group_by(plot, subplot,Species) %>%
    summarise(
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    left_join(shade_df,by=c("Species"="species")) %>% 
    group_by(plot,subplot) %>% select(plot,,subplot,Species,abundance,shade_tolerance) %>% 
    summarise(cwm_shade=sum(abundance*shade_tolerance)/sum(abundance))
  
  ### SPECIES SPECIFIC 
  #%%%%%%%%%%%%%%%%%%%%%%%
  
  # Abundance mean and variability
  
  general_sp_df <- regen_df %>%
    group_by(plot, subplot,Species) %>%
    summarise(
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    )
  
  var_sp_df <- regen_df %>%
    group_by(plot, subplot,Species, Quadrat) %>%
    summarise(
      abundance = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(plot, subplot,Species) %>%
    summarise(
      abundance_sd = sd(abundance),
      abundance_cv = sd(abundance) / mean(abundance),
      .groups = "drop"
    )
  
  # Height distribution
  shannon_height_sp_df <-  regen_df %>%
    filter(!is.na(Hclass)) %>% 
    mutate(Count=replace_na(Count,0)) %>% 
    group_by(plot, subplot,Species, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, subplot,Species) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H_height = -sum(p), .groups = "drop")
  
  # Growth rate, total and per size class
  growth_sp_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, subplot, Species) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    )
  
  growth_sp_size_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, subplot, Species, Hclass) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = Hclass,
      values_from = c(Height_increment_mn, Height_increment_cv),
      names_sep   = "_hclass"
    )
  
  # browsing intensity
  browsing_sp_df <- regen_df %>%
    filter(Hclass > 1) %>%
    group_by(plot, subplot, Species) %>%
    summarise(browsing = sum(Browsing == "Y") / n(), .groups = "drop")
  
  # ── Final join at plot/subplot ────────────────────────────────────────────────
  final_df <- general_df %>%
    left_join(var_df,          by = c("plot", "subplot")) %>%
    left_join(shannon_df,      by = c("plot", "subplot")) %>%
    left_join(browsing_df,     by = c("plot", "subplot")) %>%
    left_join(seedling_df,     by = c("plot", "subplot")) %>%
    left_join(sapling_df,     by = c("plot", "subplot")) %>%
    left_join(shannon_size_df, by = c("plot", "subplot")) %>%
    left_join(growth_df,  by = c("plot", "subplot")) %>% 
    left_join(growth_size_df,  by = c("plot", "subplot")) %>% 
    left_join(shannon_height_df, by = c("plot", "subplot"))%>% 
    left_join(shade_tol_df, by = c("plot", "subplot"))
  final_sp_df<-general_sp_df %>%
    left_join(var_sp_df,             by = c("plot", "subplot","Species")) %>%
    left_join(shannon_height_sp_df,  by = c("plot", "subplot","Species")) %>%
    left_join(growth_sp_df,          by = c("plot", "subplot","Species")) %>%
    left_join(growth_sp_size_df,     by = c("plot", "subplot","Species")) %>%
    left_join(browsing_sp_df,        by = c("plot", "subplot","Species")) 
  return(list(all=final_df,species=final_sp_df))
}

# Get regeneration metrics at the plot level
#' @param regen_df table of the data from the field

get_regen_metric_plot<-function(regen_df,shade_df){
  ### ALL SPECIES TOGETHER 
  #%%%%%%%%%%%%%%%%%%%%%%%
  
  # Abundance & richness
  general_df <- regen_df %>%
    group_by(plot) %>%
    summarise(
      richness     = n_distinct(Species),
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    )
  
  var_df <- regen_df %>%
    group_by(plot, Quadrat) %>%
    summarise(
      richness  = n_distinct(Species),
      abundance = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(plot) %>%
    summarise(
      richness_sd  = sd(richness),
      abundance_sd = sd(abundance),
      abundance_cv = sd(abundance) / mean(abundance),
      .groups = "drop"
    )
  
  # Abundance & richness - seedlings and saplings
  seedling_df <- regen_df %>%
    filter(!is.na(Species), Hclass %in% c(1, 2)) %>%
    group_by(plot) %>%
    summarise(richness_seedling     = n_distinct(Species),
              n_seedling = sum(Count), .groups = "drop")
  
  sapling_df <- regen_df %>%
    filter(!is.na(Species), !Hclass %in% c(1, 2)) %>%
    group_by(plot) %>%
    summarise(richness_sapling   = n_distinct(Species),
              n_sapling = sum(Count), .groups = "drop")
  
  # Shannon index total and per size class
  shannon_df <- regen_df %>%
    filter(!is.na(as.numeric(Count))) %>%
    group_by(plot, Species) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H = -sum(p), .groups = "drop")
  
  shannon_size_df <- regen_df %>%
    filter(!is.na(as.numeric(Count))) %>%
    group_by(plot, Species, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, Hclass) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H = -sum(p), .groups = "drop") %>%
    pivot_wider(names_from = Hclass, values_from = H, names_prefix = "H_hclass")
  
  # Growth rate, total and per size class
  growth_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    )
  
  growth_size_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, Hclass) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = Hclass,
      values_from = c(Height_increment_mn, Height_increment_cv),
      names_sep   = "_hclass"
    )
  
  # Height distribution
  shannon_height_df <-  regen_df %>%
    filter(!is.na(Hclass)) %>% 
    mutate(Count=replace_na(Count,0)) %>% 
    group_by(plot, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H_height = -sum(p), .groups = "drop")
  
  # community wighted mean for shade tolerange
  shade_tol_df<- regen_df %>%
    group_by(plot,Species) %>%
    summarise(
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    left_join(shade_df,by=c("Species"="species")) %>% 
    group_by(plot) %>% select(plot,Species,abundance,shade_tolerance) %>% 
    summarise(cwm_shade=sum(abundance*shade_tolerance)/sum(abundance))
  
  
  
  # browsing intensity
  browsing_df <- regen_df %>%
    filter(Hclass > 1) %>%
    group_by(plot) %>%
    summarise(browsing = sum(Browsing == "Y") / n(), .groups = "drop")
  
  ### SPECIES SPECIFIC 
  #%%%%%%%%%%%%%%%%%%%%%%%
  
  # Abundance mean and variability
  
  general_sp_df <- regen_df %>%
    group_by(plot,Species) %>%
    summarise(
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    )
  
  var_sp_df <- regen_df %>%
    group_by(plot,Species, Quadrat) %>%
    summarise(
      abundance = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(plot,Species) %>%
    summarise(
      abundance_sd = sd(abundance),
      abundance_cv = sd(abundance) / mean(abundance),
      .groups = "drop"
    )
  
  # Height distribution
  shannon_height_sp_df <-  regen_df %>%
    filter(!is.na(Hclass)) %>% 
    mutate(Count=replace_na(Count,0)) %>% 
    group_by(plot,Species, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(plot,Species) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H_height = -sum(p), .groups = "drop")
  
  # Growth rate, total and per size class
  growth_sp_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, Species) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    )
  
  growth_sp_size_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(plot, Species, Hclass) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = Hclass,
      values_from = c(Height_increment_mn, Height_increment_cv),
      names_sep   = "_hclass"
    )
  
  # browsing intensity
  browsing_sp_df <- regen_df %>%
    filter(Hclass > 1) %>%
    group_by(plot, Species) %>%
    summarise(browsing = sum(Browsing == "Y") / n(), .groups = "drop")
  
  # ── Final join at plot/subplot ────────────────────────────────────────────────
  final_df <- general_df %>%
    left_join(var_df,          by = "plot") %>%
    left_join(shannon_df,      by = "plot") %>%
    left_join(browsing_df,     by = "plot") %>%
    left_join(seedling_df,     by = "plot") %>%
    left_join(sapling_df,     by = "plot") %>%
    left_join(shannon_size_df, by = "plot") %>%
    left_join(growth_df,  by = "plot") %>% 
    left_join(growth_size_df,  by = "plot") %>% 
    left_join(shannon_height_df, by = "plot")%>% 
    left_join(shade_tol_df, by = "plot")
  final_sp_df<-general_sp_df %>%
    left_join(var_sp_df,             by = c("plot","Species")) %>%
    left_join(shannon_height_sp_df,  by = c("plot","Species")) %>%
    left_join(growth_sp_df,          by = c("plot","Species")) %>%
    left_join(growth_sp_size_df,     by = c("plot","Species")) %>%
    left_join(browsing_sp_df,        by = c("plot","Species")) 
  return(list(all=final_df,species=final_sp_df))
}


### ADULT METRICS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

get_adults <- function(inventory_df){
  
  plot_stand_metrics <- inventory_df %>%
    mutate(
      alive = Dead17 == 0 & Cut17 == 0,
      basal_area_m2 = pi * (DBH2017 / 200)^2,  # DBH in cm -> radius in m
      alive_big = Dead17 == 0 & Cut17 == 0 & DBH2017 > 20
    ) %>%
    group_by(plot) %>%
    summarise(
      n_trees_total = n(),
      n_trees_alive = sum(alive, na.rm = TRUE),
      n_dead = sum(Dead17 == 1, na.rm = TRUE),
      n_cut = sum(Cut17 == 1, na.rm = TRUE),
      n_recruits = sum(Recrut17 == 1, na.rm = TRUE),
      
      species_richness_total = n_distinct(SpeciesCode[!is.na(SpeciesCode)]),
      species_richness_alive = n_distinct(SpeciesCode[alive & !is.na(SpeciesCode)]),
      species_richness_dominant = n_distinct(SpeciesCode[alive_big & !is.na(SpeciesCode)]),
      
      mean_dbh_alive_cm = mean(DBH2017[alive], na.rm = TRUE),
      median_dbh_alive_cm = median(DBH2017[alive], na.rm = TRUE),
      max_dbh_alive_cm = max(DBH2017[alive], na.rm = TRUE),
      sd_dbh_alive_cm = sd(DBH2017[alive], na.rm = TRUE),
      
      mean_height_alive_m = mean(Htot2017[alive], na.rm = TRUE),
      max_height_alive_m = max(Htot2017[alive], na.rm = TRUE),
      
      basal_area_alive_m2 = sum(basal_area_m2[alive], na.rm = TRUE),
      mean_basal_area_alive_m2 = mean(basal_area_m2[alive], na.rm = TRUE),
      
      biomass_alive = sum(Biomass2017[alive], na.rm = TRUE),
      biomass_total = sum(Biomass2017, na.rm = TRUE),
      virtual_biomass_total = sum(VirtualBiomass2017, na.rm = TRUE),
      
      mortality_rate = n_dead / n_trees_total,
      cutting_rate = n_cut / n_trees_total,
      recruitment_rate = n_recruits / n_trees_total,
      
      .groups = "drop"
    )
}


### SEGMENTED PLOTS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%



file.segmented<-function(path_plot,
                         plot_name,
                         folders=c("regen","hesitation")){
  path_trees=file.path(path_plot,
                       paste0(plot_name,"-processing"),
                       "trees",folders)
  is.segmented=sum(file.exists(path_trees))
  
  if(is.segmented>0){
    file=paste0("output/dbh_seg_",plot_name,".csv")
  }else file=NULL
  return(file)
}

get_pc_norm_regen<-function(path_plot,
                            plot_name,
                            file.segmented,
                            dtm,
                            bbox){
  if(!is.null(file.segmented)){
    
    ## gather trees txt files -----------------------
    dbh_seg=read.csv(file.segmented) %>% 
      select(-any_of("X")) %>%
      mutate(across(-any_of("File"), as.numeric)) %>% 
      filter(DBH<0.12|is.na(DBH)) %>% 
      mutate(file_long=case_when(grepl("hesit",File)~file.path(path_plot,
                                                               paste0(plot_name,"-processing"),
                                                               "trees","hesitation",File),
                                 grepl("regen",File)~file.path(path_plot,
                                                               paste0(plot_name,"-processing"),
                                                               "trees","regen",File)))
    
    files=dbh_seg$file_long
    
    ## load dtm ------------------------------------
    dtm_rast <- rast(dtm)
    
    
    ## load trees files and merge in LAS -----------
    read_xyz <- function(fn) {
      dt <- read_tree_pc(fn)
      setDT(dt)
      dt <- dt[
        is.finite(X) &
          is.finite(Y) &
          is.finite(Z)
      ]
      dt[]
    }
    
    tree_list <- lapply(files, read_xyz)
    tree_list <- Filter(Negate(is.null), tree_list)
    merged <- rbindlist( tree_list,
                         use.names = TRUE,
                         fill = TRUE)
    las_merged=LAS(merged[, c("X", "Y", "Z")])
    
    las_clipped <- clip_rectangle(las_merged,
                                  xleft   = bbox$xmin,
                                  ybottom = bbox$ymin,
                                  xright  = bbox$xmax,
                                  ytop    = bbox$ymax)
    
    
    # normalize pc height ---------------------------
    las_norm <- normalize_height(las_clipped, dtm_rast)
    
    las_norm@data$X <- las_norm@data$X - min(las_norm@data$X)
    las_norm@data$Y <- las_norm@data$Y - min(las_norm@data$Y)
    
    las_norm <- las_update(las_norm)
    
    # las_clip <- filter_poi(las_norm, Z < 5)
    
    path_las=paste0("output/",plot_name,"_pcnorm_regen.rdata")
    save(las_norm,file=path_las)
  }else{path_las=NULL}
  return(path_las)
}

# 
# pc_norm_seg=tar_read(pc_norm_seg_GER02)
# rb_norm=tar_read(rb_norm_GER02)
# voxnorm=tar_read(voxnorm_GER02)
# bbox=tar_read(bbox_GER02)
# subplot_extent=tar_read(subplot_extent_GER02)
# plot_name="GER02"
# sliced=FALSE
get_below_canopy_complexity <- function(pc_norm,
                                        rb_norm,
                                        voxnorm,
                                        bbox,
                                        subplot_extent,
                                        plot_name,
                                        voxel_size = 0.25,
                                        chm_res    = 0.1) {
  
  pc_norm <- loadRData(pc_norm)
  list_rb <- loadRData(rb_norm)
  voxpad  <- loadRData(voxnorm)
  
  chm <- rasterize_canopy(pc_norm, res = chm_res, algorithm = p2r())   # full canopy
  
  # ---- below-canopy keep-mask for ONE voxel array, on its own grid -----------
  # keep[i, j, a] = TRUE when layer a's top (a * voxel_size) is at/below canopy.
  make_keep <- function(arr) {
    tmpl <- rast(arr[, , 1]); ext(tmpl) <- ext(chm); crs(tmpl) <- crs(chm)
    cv   <- resample(chm, tmpl, method = "max")
    tl   <- as.matrix(cv, wide = TRUE) / voxel_size          # [nr, nc] for THIS array
    kp   <- outer(tl, seq_len(dim(arr)[3]), function(t, a) a <= t)
    kp[is.na(kp)] <- FALSE
    list(keep = kp, no_canopy = is.na(tl))
  }
  
  rb1    <- list_rb[[1]]
  km_rb  <- make_keep(rb1)       # mask shaped like the rb arrays
  km_pad <- make_keep(voxpad)    # mask shaped like voxpad
  
  # ---- below-canopy PAD (time-invariant): matched raster + cube for FHD ------
  padkeep   <- voxpad * km_pad$keep
  pad_below <- rowSums(padkeep, dims = 2, na.rm = TRUE)
  # pad_below[km_pad$no_canopy] <- NA
  pad_slice <- rast(pad_below); ext(pad_slice) <- ext(chm); crs(pad_slice) <- crs(chm)
  pad_matched <- resample(pad_slice, chm, method = "bilinear")
  names(pad_matched) <- "pad_sum"
  
  padkeep_cube <- rast(padkeep); ext(padkeep_cube) <- ext(chm); crs(padkeep_cube) <- crs(chm)
  
  # ---- spatial units: plot (geo = NULL) + subplots ---------------------------
  sub_geo <- lapply(subplot_extent, function(pe)
    shift(vect(pe), dx = -bbox$xmin, dy = -bbox$ymin))
  units <- c(list(list(name = "plot", geo = NULL)),
             lapply(seq_along(subplot_extent), function(sp)
               list(name = names(subplot_extent)[sp], geo = sub_geo[[sp]])))
  
  cropu <- function(r, geo) if (is.null(geo)) r else crop(r, geo)
  
  # ---- pre-crop the time-invariant inputs once per unit ----------------------
  unit_static <- lapply(units, function(u) {
    chm_u     <- cropu(chm, u$geo)
    chm_num_u <- chm_u; chm_num_u[is.na(chm_num_u)] <- 0
    list(
      name     = u$name,
      geo      = u$geo,
      chm      = chm_u,
      chm_num  = chm_num_u,
      pad      = cropu(pad_matched,  u$geo),
      pad_cube = cropu(padkeep_cube, u$geo),
      pc_gap   = if (is.null(u$geo)) pc_norm else clip_roi(pc_norm, sf::st_as_sf(u$geo))
    )
  })
  
  n_rows <- length(list_rb) * length(units)
  rows <- vector("list", n_rows)
  k <- 0L
  
  # ---- time loop: only rb (light) varies -------------------------------------
  for (t in seq_along(list_rb)) {
    rb_arr <- list_rb[[t]]
    time   <- as.numeric(names(list_rb)[t])
    
    rb_below <- rowSums(rb_arr * km_rb$keep, dims = 2, na.rm = TRUE)
    # rb_below[km_rb$no_canopy] <- NA
    rb_slice <- rast(rb_below); ext(rb_slice) <- ext(chm); crs(rb_slice) <- crs(chm)
    rb_matched <- resample(rb_slice, chm, method = "bilinear")
    names(rb_matched) <- "rb"
    thr <- quantile(values(rb_matched), na.rm = TRUE, probs = 0.005)
    rb_matched[rb_matched < thr] <- NA
    rb_matched_mask <- mask(rb_matched, chm)
    
    for (us in unit_static) {
      rb_u  <- cropu(rb_matched,      us$geo)
      rbm_u <- cropu(rb_matched_mask, us$geo)
      
      m <- .slice_metrics(us$chm, us$chm_num, rb_u, rbm_u,
                          us$pad, us$pad_cube, us$pc_gap,
                          slice_up = Inf, voxel_size = voxel_size)
      
      k <- k + 1L
      rows[[k]] <- .make_row(plot_name, us$name, time,
                             slice_num = 0, slice_low = 0, slice_up = NA_real_, m)
    }
  }
  
  dplyr::bind_rows(rows)
}



#### GINI AND BIG LEAF ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%

.gini <- function(x) {
  x <- x[is.finite(x)]
  x <- x[x >= 0]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  s <- sum(x)
  if (s <= 0) return(NA_real_)
  x <- sort(x)
  (2 * sum(seq_len(n) * x) / (n * s)) - (n + 1) / n
}

## --------------------------------------------------------------------------
## Array 3D [row, col, z] -> SpatRaster multicouche (1 couche = 1 tranche z)
## L'emprise est calée sur l'origine (0,0), comme les géométries décalées
## par -bbox$xmin / -bbox$ymin.
## --------------------------------------------------------------------------
.cube_rast <- function(a, vox, crs_ref = "") {
  r <- terra::rast(a)
  terra::ext(r) <- c(0, terra::ncol(r) * vox, 0, terra::nrow(r) * vox)
  if (nzchar(crs_ref)) terra::crs(r) <- crs_ref
  names(r) <- paste0("z", seq_len(terra::nlyr(r)))
  r
}

## --------------------------------------------------------------------------
## Métriques pour une zone (plot entier ou subplot déjà croppé)
##   rb_cube  : SpatRaster, énergie reçue par voxel
##   pad_cube : SpatRaster, PAD par voxel (mêmes dimensions)
##   z_band   : indices de couches pour le Gini horizontal
##   vox      : taille de voxel (m)
##   k        : coefficient d'extinction du big-leaf
##   top_idx  : indice de la couche sommitale (NULL = déduit du PAD)
## --------------------------------------------------------------------------
.light_metrics <- function(rb_cube, pad_cube, z_band, vox, k,
                           top_idx = NULL, keep_profiles = FALSE,
                           pad_is_density = TRUE) {
  
  rb_prof  <- as.numeric(terra::global(rb_cube,  "mean", na.rm = TRUE)[, 1])
  pad_prof <- as.numeric(terra::global(pad_cube, "mean", na.rm = TRUE)[, 1])
  nz <- length(rb_prof)
  
  ## ---- sommet de canopée -------------------------------------------------
  if (is.null(top_idx)) {
    w <- which(pad_prof > 0 & is.finite(pad_prof))
    top <- if (length(w)) max(w) else nz
  } else {
    top <- min(top_idx, nz)
  }
  zz     <- seq_len(top)
  rb_p   <- rb_prof[zz]
  pad_p  <- pad_prof[zz]
  z_mid  <- (zz - 0.5) * vox
  
  ## ---- 1. Gini horizontal ------------------------------------------------
  z_band <- z_band[z_band >= 1 & z_band <= nz]
  band   <- terra::app(terra::subset(rb_cube, z_band), fun = sum, na.rm = TRUE)
  gini_h <- .gini(terra::values(band, mat = FALSE))
  
  ## ---- 2. Gini vertical --------------------------------------------------
  gini_v <- .gini(rb_p)
  
  ## ---- 3. Écart au big-leaf ---------------------------------------------
  ## PAD stocké en densité (m2/m3) -> surface par couche = densité * épaisseur
  lai_layer <- if (pad_is_density) pad_p * vox else pad_p
  ## PAI total moyen de la zone (m2/m2), redistribué uniformément sur [0, top]
  pai    <- sum(lai_layer, na.rm = TRUE)
  pad_bl <- pai / top                     # surface foliaire par couche
  ## surface foliaire cumulée depuis le sommet, au milieu de chaque couche
  Lcum   <- pad_bl * (top - zz + 0.5)
  rb_bl  <- exp(-k * Lcum)
  
  ## normalisation : 1 dans la couche sommitale pour les deux profils
  ref <- rb_p[top]
  if (!is.finite(ref) || ref <= 0) {
    rb_obs <- rep(NA_real_, top)
  } else {
    rb_obs <- rb_p / ref
  }
  rb_bl <- rb_bl / rb_bl[top]
  
  d          <- rb_obs - rb_bl
  diff_abs   <- sum(abs(d), na.rm = TRUE)
  diff_mean  <- mean(abs(d), na.rm = TRUE)
  diff_sign  <- sum(d, na.rm = TRUE)
  diff_rmse  <- sqrt(mean(d^2, na.rm = TRUE))
  
  if (top >= 2L) {
    sl_obs      <- diff(rb_obs) / vox
    sl_bl       <- diff(rb_bl)  / vox
    ds          <- sl_obs - sl_bl
    slope_abs   <- sum(abs(ds), na.rm = TRUE)
    slope_mean  <- mean(abs(ds), na.rm = TRUE)
    slope_sign  <- sum(ds, na.rm = TRUE)
  } else {
    slope_abs <- slope_mean <- slope_sign <- NA_real_
  }
  
  out <- data.frame(
    gini_h_rb          = gini_h,
    gini_v_rb          = gini_v,
    bl_diff_abs        = diff_abs,
    bl_diff_abs_mean   = diff_mean,
    bl_diff_signed     = diff_sign,
    bl_diff_rmse       = diff_rmse,
    bl_slope_abs       = slope_abs,
    bl_slope_abs_mean  = slope_mean,
    bl_slope_signed    = slope_sign,
    pai                = pai,
    z_top              = top * vox,
    n_layers           = top,
    rb_top_mean        = ref
  )
  
  if (keep_profiles) {
    attr(out, "profile") <- data.frame(
      z         = z_mid,
      rb_mean   = rb_p,
      pad       = pad_p,
      lai_layer = lai_layer,
      rb_obs  = rb_obs,
      rb_bl   = rb_bl
    )
  }
  out
}

## --------------------------------------------------------------------------
## Fonction principale
## --------------------------------------------------------------------------
get_light_metrics <- function(rb_norm,
                              voxnorm,
                              bbox,
                              subplot_extent,
                              plot_name,
                              vox_size      = 0.25,
                              horiz_low     = 0,      # m
                              horiz_up      = 4.5,    # m
                              k             = 0.5,    # extinction big-leaf
                              z_top         = NULL,   # m, NULL = déduit du PAD
                              keep_profiles = FALSE,
                              pad_is_density = TRUE) { # PAD en m2/m3
  
  list_rb <- loadRData(rb_norm) 
  voxpad  <- loadRData(voxnorm)
  
  if (!identical(dim(voxpad)[1:2], dim(list_rb[[1]])[1:2])){
    if(abs(dim(voxpad)[1]-dim(list_rb[[1]])[1])<=2 &
       abs(dim(voxpad)[2]-dim(list_rb[[1]])[2])<=2){
      dim_x=min(dim(voxpad)[1],dim(list_rb[[1]])[1])
      dim_y=min(dim(voxpad)[2],dim(list_rb[[1]])[2])
      voxpad=voxpad[1:dim_x,1:dim_y,]
      for(l in seq_along(list_rb)){
        list_rb[[l]]=list_rb[[l]][1:dim_x,1:dim_y,]
      }
    }else{
      stop("voxpad and rb_norm do not share the same horizontal grid")
    }
  } 
  
  ## indices de couches pour le Gini horizontal (conversion 1-based)
  z_band <- seq(floor(horiz_low / vox_size) + 1L,
                ceiling(horiz_up / vox_size))
  top_idx <- if (is.null(z_top)) NULL else round(z_top / vox_size)
  
  ## géométries des subplots (constantes) ----------------------------------
  sub_geo <- NULL
  crs_ref <- ""
  sub_geo <- lapply(subplot_extent, function(pe) {
    terra::shift(terra::vect(pe), dx = -bbox$xmin, dy = -bbox$ymin)
  })
  crs_ref <- terra::crs(sub_geo[[1]])
  
  
  pad_cube <- .cube_rast(voxpad, vox_size, crs_ref)
  pad_sub  <- if (is.null(sub_geo)) NULL else
    lapply(sub_geo, function(g) terra::crop(pad_cube, g))
  
  n_rows <- length(list_rb) * (1L + length(sub_geo))
  rows   <- vector("list", n_rows)
  profs  <- vector("list", n_rows)
  k_row  <- 0L
  
  for (t in seq_along(list_rb)) {
    
    time    <- as.numeric(names(list_rb)[t])
    rb_cube <- .cube_rast(list_rb[[t]], vox_size, crs_ref)
    
    ## ---- échelle plot ----------------------------------------------------
    m <- .light_metrics(rb_cube, pad_cube, z_band, vox_size, k,
                        top_idx, keep_profiles, pad_is_density)
    k_row <- k_row + 1L
    rows[[k_row]] <- data.frame(plot = plot_name, scale = "plot",
                                subplot = "plot", time = time,
                                band_low = horiz_low, band_up = horiz_up,
                                k = k, m, row.names = NULL)
    if (keep_profiles)
      profs[[k_row]] <- transform(attr(m, "profile"),
                                  plot = plot_name, subplot = "plot",
                                  time = time)
    
    ## ---- échelle subplot -------------------------------------------------
    for (sp in seq_along(sub_geo)) {
      rb_s  <- terra::crop(rb_cube, sub_geo[[sp]])
      m_s   <- .light_metrics(rb_s, pad_sub[[sp]], z_band, vox_size, k,
                              top_idx, keep_profiles, pad_is_density)
      k_row <- k_row + 1L
      rows[[k_row]] <- data.frame(plot = plot_name, scale = "subplot",
                                  subplot = names(subplot_extent)[sp],
                                  time = time,
                                  band_low = horiz_low, band_up = horiz_up,
                                  k = k, m_s, row.names = NULL)
      if (keep_profiles)
        profs[[k_row]] <- transform(attr(m_s, "profile"),
                                    plot = plot_name,
                                    subplot = names(subplot_extent)[sp],
                                    time = time)
    }
  }
  
  res <- dplyr::bind_rows(rows)
  if (keep_profiles) attr(res, "profiles") <- dplyr::bind_rows(profs)
  res
}