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
# bbox=tar_read(bbox_GER13)
# subplot_extent=tar_read(subplot_extent_GER13)
# plot_name="GER13"
# nstrat=4
get_complexity_rb<-function(pc_norm,
                            rb_norm,
                            voxnorm,
                            bbox,
                            subplot_extent,
                            plot_name,
                            nstrat=4){
  
  pc_norm<-loadRData(pc_norm)
  list_rb<-loadRData(rb_norm)
  voxpad<-loadRData(voxnorm)
  
  slice_list = seq(from = 0.5, by = 1, length.out = nstrat+1)
  summary_plot=data.frame(plot=character(),
                          subplot=character(),
                          time=numeric(),
                          slice_num=numeric(),
                          slice_low=numeric(),
                          slice_up=numeric(),
                          h_mean_na=numeric(),
                          h_med_na=numeric(),
                          h_sd_na=numeric(),
                          h_cv_na=numeric(),
                          h_mean=numeric(),
                          h_med=numeric(),
                          h_sd=numeric(),
                          h_cv=numeric(),
                          empty_gf=numeric(),
                          occupied_gf=numeric(),
                          rb_mean=numeric(),
                          rb_med=numeric(),
                          rb_sd=numeric(),
                          rb_cv=numeric(),
                          rb_mean_mask=numeric(),
                          rb_med_mask=numeric(),
                          rb_sd_mask=numeric(),
                          rb_cv_mask=numeric(),
                          coef_lin=numeric(),
                          pad_tot=numeric(),
                          pad_mean=numeric(),
                          pad_sd=numeric(),
                          pad_cv=numeric())
  
  
  for(t in seq_along(list_rb)){
    rb_norm=list_rb[[t]]
    time=as.numeric(names(list_rb)[t])
    
    for(i in 1:nstrat){
      slice_low=slice_list[i]
      slice_up=slice_list[i+1]
      
      ### metric at the plot scale for height
      pc_slice=filter_poi(pc_norm, Z >= slice_low & Z < slice_up)
      pc_slice_chm <- rasterize_canopy(pc_slice, res = 0.1, algorithm = p2r())
      pc_slice_chm[pc_slice_chm > (slice_up-0.01)] <- NA
      
      # excluding zone without vegetation
      h_mean_na=mean(values(pc_slice_chm),na.rm=TRUE)
      h_med_na=median(values(pc_slice_chm),na.rm=TRUE)
      h_sd_na=sd(values(pc_slice_chm),na.rm=TRUE)
      h_cv_na=h_sd_na/h_mean_na
      
      # accounting for empty space
      pc_slice_chm_num=pc_slice_chm
      pc_slice_chm_num[is.na(pc_slice_chm_num)] <- 0
      
      h_mean=mean(values(pc_slice_chm_num),na.rm=TRUE)
      h_med=median(values(pc_slice_chm_num),na.rm=TRUE)
      h_sd=sd(values(pc_slice_chm_num),na.rm=TRUE)
      h_cv=h_sd/h_mean
      
      ### metric at the plot scale for light
      rb_slice=rast(apply(rb_norm[,,(slice_low/0.25):(slice_up/0.25)], c(1, 2), sum))
      ext(rb_slice)  <- ext(pc_slice_chm)
      crs(rb_slice)  <- crs(pc_slice_chm)
      rb_slice_matched <- resample(rb_slice, pc_slice_chm, method = "bilinear")
      names(rb_slice_matched)="rb"
      rb_slice_matched[rb_slice_matched<quantile(values(rb_slice_matched),
                                                 na.rm=TRUE,probs=0.005)] <-NA
      rb_slice_matched_mask = mask(rb_slice_matched,pc_slice_chm)
      
      
      rb_mean=mean(values(rb_slice_matched),na.rm=TRUE)
      rb_med=median(values(rb_slice_matched),na.rm=TRUE)
      rb_sd=sd(values(rb_slice_matched),na.rm=TRUE)
      rb_cv=rb_sd/rb_mean
      rb_mean_mask=mean(values(rb_slice_matched_mask),na.rm=TRUE)
      rb_med_mask=median(values(rb_slice_matched_mask),na.rm=TRUE)
      rb_sd_mask=sd(values(rb_slice_matched_mask),na.rm=TRUE)
      rb_cv_mask=rb_sd_mask/rb_mean_mask
      
      df=as.data.frame(c(pc_slice_chm,
                         rb_slice_matched),xy=TRUE,
                       na.rm=TRUE) %>% 
        filter(Z<(slice_up-0.01))
      reglin=lm(rb~Z,data=df)
      coef = reglin$coefficients["Z"][[1]]
      
      
      ## gap filling
      gap_filling=get_gap_filling(pc_slice,voxel_size = 0.25)
      
      
      ## pad
      pad_slice=rast(apply(voxpad[,,(slice_low/0.25):(slice_up/0.25)], c(1, 2), sum))
      ext(pad_slice)  <- ext(pc_slice_chm)
      crs(pad_slice)  <- crs(pc_slice_chm)
      pad_slice_matched <- resample(pad_slice, pc_slice_chm, method = "bilinear")
      names(pad_slice_matched)="pad_sum"
      
      pad_tot=sum(values(pad_slice_matched),na.rm=TRUE)
      pad_mean=mean(values(pad_slice_matched),na.rm=TRUE)
      pad_sd=sd(values(pad_slice_matched),na.rm=TRUE)
      pad_cv=pad_sd/pad_mean
      
      
      summary_plot[nrow(summary_plot) + 1, ] <- list(
        plot_name, "plot",time,i, slice_low, slice_up,
        h_mean_na, h_med_na, h_sd_na, h_cv_na,
        h_mean, h_med, h_sd, h_cv,
        gap_filling$empty,gap_filling$occupied,
        rb_mean, rb_med, rb_sd, rb_cv,
        rb_mean_mask, rb_med_mask, rb_sd_mask, rb_cv_mask,
        coef,
        pad_tot,pad_mean,pad_sd,pad_cv
      )
      
      ### do the same for subplots
      for(sp in seq_along(subplot_extent)){
        plot_ext=subplot_extent[[sp]]
        plot_ext_trans <- shift(vect(plot_ext), 
                                dx = -bbox$xmin, 
                                dy = -bbox$ymin)
        pc_sub_slice_chm=crop(pc_slice_chm,plot_ext_trans)
        pc_sub_slice_chm_num=crop(pc_slice_chm_num,plot_ext_trans)
        rb_sub_slice_matched=crop(rb_slice_matched,plot_ext_trans)
        rb_sub_slice_matched_mask=crop(rb_slice_matched_mask,plot_ext_trans)
        pad_sub_slice_matched=crop(pad_slice_matched,plot_ext_trans)
        
        
        ### height metrics
        h_mean_s_na=mean(values(pc_sub_slice_chm),na.rm=TRUE)
        h_med_s_na=median(values(pc_sub_slice_chm),na.rm=TRUE)
        h_sd_s_na=sd(values(pc_sub_slice_chm),na.rm=TRUE)
        h_cv_s_na=h_sd_s_na/h_mean_s_na
        
        ### height metrics, na set as 0
        h_mean_s=mean(values(pc_sub_slice_chm_num),na.rm=TRUE)
        h_med_s=median(values(pc_sub_slice_chm_num),na.rm=TRUE)
        h_sd_s=sd(values(pc_sub_slice_chm_num),na.rm=TRUE)
        h_cv_s=h_sd_s/h_mean_s
        
        ### rb metrics 
        rb_mean_s=mean(values(rb_sub_slice_matched),na.rm=TRUE)
        rb_med_s=median(values(rb_sub_slice_matched),na.rm=TRUE)
        rb_sd_s=sd(values(rb_sub_slice_matched),na.rm=TRUE)
        rb_cv_s=rb_sd_s/rb_mean_s
        
        ### rb metrics masked
        rb_mean_mask_s=mean(values(rb_sub_slice_matched_mask),na.rm=TRUE)
        rb_med_mask_s=median(values(rb_sub_slice_matched_mask),na.rm=TRUE)
        rb_sd_mask_s=sd(values(rb_sub_slice_matched_mask),na.rm=TRUE)
        rb_cv_mask_s=rb_sd_mask_s/rb_mean_mask_s
        
        names(rb_sub_slice_matched)="rb"
        df_s=as.data.frame(c(pc_sub_slice_chm,
                             rb_sub_slice_matched),xy=TRUE,
                           na.rm=TRUE) %>% 
          filter(Z<(slice_up-0.01))
        if(dim(df_s)[1]>3){
          reglin_s=lm(rb~Z,data=df_s)
          coef_s = reglin_s$coefficients["Z"][[1]]
        }else{coef_s=NA}
        
        ## gap filling
        gap_filling_s=get_gap_filling(clip_roi(pc_slice,
                                             sf::st_as_sf(plot_ext_trans)),
                                    voxel_size = 0.25)
        
        ## pad metrics
        pad_tot_s=sum(values(pad_sub_slice_matched),na.rm=TRUE)
        pad_mean_s=mean(values(pad_sub_slice_matched),na.rm=TRUE)
        pad_sd_s=sd(values(pad_sub_slice_matched),na.rm=TRUE)
        pad_cv_s=pad_sd_s/pad_mean_s
        
        
        summary_plot[nrow(summary_plot) + 1, ] <- list(
          plot_name, names(subplot_extent[sp]),time, i, slice_low, slice_up,
          h_mean_s_na, h_med_s_na, h_sd_s_na, h_cv_s_na,
          h_mean_s, h_med_s, h_sd_s, h_cv_s,
          gap_filling_s$empty,gap_filling_s$occupied,
          rb_mean_s, rb_med_s, rb_sd_s, rb_cv_s,
          rb_mean_mask_s, rb_med_mask_s, rb_sd_mask_s, rb_cv_mask_s,
          coef_s,
          pad_tot_s,pad_mean_s,pad_sd_s,pad_cv_s
        )
      }
    }
  }
  return(summary_plot)
}


### REGENERATION METRICS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get regeneration metrics at the subplot level
#' @param regen_df table of the data from the field
get_regen_metric_subplot<-function(regen_df){
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
    left_join(shannon_height_df, by = c("plot", "subplot"))
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

get_regen_metric_plot<-function(regen_df){
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
    left_join(shannon_height_df, by = "plot")
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