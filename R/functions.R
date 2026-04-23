loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


get_rb<-function(dart_folder,
                 plot_name,
                 hours=c(8,10,12,14,16,18)){
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
  
  path_sum=paste0("output/",plot_name,"_rb_sum.rdata")
  save(sum_arr,file=path_sum)
  path_sd=paste0("output/",plot_name,"_rb_sd.rdata")
  save(sd_arr,file=path_sd)
  
  return(c(path_sum,path_sd))
}


normalize_rb<-function(plot_name,
                       user,
                       rb_path,
                       rb_type){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  
  rb_rast=loadRData(rb_path)

  nx <- dim(rb_rast)[2]
  ny <- dim(rb_rast)[1]
  nz <- dim(rb_rast)[3]
  
  dtm_rast <- rast(file.path(path_folder,
                             list.files(path_folder,pattern="dtm.tif")))

  xmin_vox <- ext(dtm_rast)[1]
  ymin_vox <- ext(dtm_rast)[3]
  zmin_vox <- min(dtm_rast)
  
  dx <- 0.25
  dy <- 0.25
  dz <- 0.25
  
  x_centers <- xmin_vox + (0:(nx - 1)) * dx + dx / 2
  y_centers <- ymin_vox + (0:(ny - 1)) * dy + dy / 2
  
  xy <- expand.grid(x = x_centers, y = y_centers)

  dtm_vals <- terra::extract(dtm_rast, xy)[, 2]
  dtm_mat <- matrix(dtm_vals, nrow = nx, ncol = ny)
  
  k_ground <- floor((dtm_mat - min(dtm_vals)) / dz) + 1
  
  k_ground[k_ground < 1] <- 1
  k_ground[k_ground > nz] <- nz
  
  maxE=quantile(rb_rast,probs=0.99)[[1]]
  arr_norm <- array(maxE, dim = dim(rb_rast))
  
  for (i in 1:nx) {
    for (j in 1:ny) {
      kg <- k_ground[i, j]
      
      if (!is.na(kg) && kg < nz) {
        new_col <- rb_rast[j, i, (kg + 1):nz]
        arr_norm[j, i, 1:length(new_col)] <- new_col
      }
    }
  }
  
  path_rast_norm=paste0("output/",plot_name,"_",rb_type,".rdata")
  save(arr_norm,file=path_rast_norm)
  return(path_rast_norm)
}

get_pc_norm<-function(plot_name,
                      user,
                      subfolders=c("adults","regen","hesitation","dead-incomplete","tree-parts")
){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  path_trees=file.path(path_folder,"trees")

  ### Create one point cloud ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
  ### DTM LAS ####
  #%%%%%%%%%%%%%%%
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


get_plot_offset<-function(plot_name,
                          user){
  path_plot=paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")
  path_folder=file.path(path_plot,
                        paste0(plot_name,"-processing"))
  dtm_rast <- rast(file.path(path_folder,
                             list.files(path_folder,pattern="dtm.tif")))
  xoffset=ext(dtm_rast)$xmin
  yoffset=ext(dtm_rast)$ymin
  return(c(xoffset,yoffset))
}


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

get_complexity_rb<-function(pc_norm,
                            rb_norm,
                            xyoffset,
                            subplot_extent,
                            plot_name,
                            nstrat=4){
  pc_norm<-loadRData(pc_norm)
  rb_norm<-loadRData(rb_norm)

  slice_list = seq(from = 0.5, by = 1, length.out = nstrat+1)
  summary_plot=data.frame(plot=character(),
                          subplot=character(),
                          slice_num=numeric(),
                          slice_low=numeric(),
                          slice_up=numeric(),
                          h_mean=numeric(),
                          h_med=numeric(),
                          h_sd=numeric(),
                          h_cv=numeric(),
                          rb_mean=numeric(),
                          rb_med=numeric(),
                          rb_sd=numeric(),
                          rb_cv=numeric(),
                          rb_mean_mask=numeric(),
                          rb_med_mask=numeric(),
                          rb_sd_mask=numeric(),
                          rb_cv_mask=numeric(),
                          coef_lin=numeric())
  for(i in 1:nstrat){
    slice_low=slice_list[i]
    slice_up=slice_list[i+1]
    
    ### metric at the plot scale for height
    pc_slice=filter_poi(pc_norm, Z >= slice_low & Z < slice_up)
    pc_slice_chm <- rasterize_canopy(pc_slice, res = 0.1, algorithm = p2r())
    pc_slice_chm[pc_slice_chm > (slice_up-0.01)] <- NA
    
    
    h_mean=mean(values(pc_slice_chm),na.rm=TRUE)
    h_med=median(values(pc_slice_chm),na.rm=TRUE)
    h_sd=sd(values(pc_slice_chm),na.rm=TRUE)
    h_cv=h_sd/h_mean
    
    
    ### metric at the plot scale for light
    rb_slice=rast(apply(rb_norm[,,(slice_low/0.254):(slice_up/0.25)], c(1, 2), sum))
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
                       rb_slice_matched),xy=TRUE) %>% 
      filter(across(everything(),
                    ~!is.na(.))) %>% 
      filter(Z<(slice_up-0.01))
    reglin=lm(rb~Z,data=df)
    coef = reglin$coefficients["Z"][[1]]
    
    
    summary_plot[nrow(summary_plot) + 1, ] <- list(
      plot_name, "plot",i, slice_low, slice_up,
      h_mean, h_med, h_sd, h_cv,
      rb_mean, rb_med, rb_sd, rb_cv,
      rb_mean_mask, rb_med_mask, rb_sd_mask, rb_cv_mask,
      coef
    )
    
    ### do the same for subplots
    for(sp in seq_along(subplot_extent)){
      plot_ext=subplot_extent[[sp]]
      plot_ext_trans <- shift(vect(plot_ext), 
                              dx = -xyoffset[[1]], 
                              dy =-xyoffset[[2]])
      pc_sub_slice_chm=crop(pc_slice_chm,plot_ext_trans)
      rb_sub_slice_matched=crop(rb_slice_matched,plot_ext_trans)
      rb_sub_slice_matched_mask=crop(rb_slice_matched_mask,plot_ext_trans)
      
      ### height metrics
      h_mean_s=mean(values(pc_sub_slice_chm),na.rm=TRUE)
      h_med_s=median(values(pc_sub_slice_chm),na.rm=TRUE)
      h_sd_s=sd(values(pc_sub_slice_chm),na.rm=TRUE)
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
                           rb_sub_slice_matched),xy=TRUE) %>% 
        filter(across(everything(),
                      ~!is.na(.))) %>% 
        filter(Z<(slice_up-0.01))
      reglin_s=lm(rb~Z,data=df_s)
      coef_s = reglin_s$coefficients["Z"][[1]]
      
      
      summary_plot[nrow(summary_plot) + 1, ] <- list(
        plot_name, names(subplot_extent[sp]), i, slice_low, slice_up,
        h_mean_s, h_med_s, h_sd_s, h_cv_s,
        rb_mean_s, rb_med_s, rb_sd_s, rb_cv_s,
        rb_mean_mask_s, rb_med_mask_s, rb_sd_mask_s, rb_cv_mask_s,
        coef_s
      )
    }
  }
  return(summary_plot)
}

get_regen_metric<-function(regen_df){
  general_df <- regen_df %>%
    group_by(Plot, Subplot) %>%
    summarise(
      richness     = n_distinct(Species),
      abundance    = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    )
  
  var_df <- regen_df %>%
    group_by(Plot, Subplot, Quadrat) %>%
    summarise(
      richness  = n_distinct(Species),
      abundance = sum(as.numeric(Count), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Plot, Subplot) %>%
    summarise(
      richness_sd  = sd(richness),
      abundance_cv = sd(abundance) / mean(abundance),
      .groups = "drop"
    )
  
  shannon_df <- regen_df %>%
    filter(!is.na(as.numeric(Count))) %>%
    group_by(Plot, Subplot, Species) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(Plot, Subplot) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H = -sum(p), .groups = "drop")
  
  browsing_df <- regen_df %>%
    filter(Hclass > 1) %>%
    group_by(Plot, Subplot) %>%
    summarise(browsing = sum(Browsing == "Y") / n(), .groups = "drop")
  
  seedling_df <- regen_df %>%
    filter(!is.na(Species), Hclass %in% c(1, 2)) %>%
    group_by(Plot, Subplot) %>%
    summarise(n_seedling = sum(Count), .groups = "drop")
  
  growth_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(Plot, Subplot, Species) %>%
    summarise(
      Height_increment_mn = mean(Height_increment, na.rm = TRUE),
      Height_increment_cv = sd(Height_increment,   na.rm = TRUE) / Height_increment_mn,
      .groups = "drop"
    )
  
  # ── Tables at Plot/Subplot/Hclass level → pivot wide before joining ──────────
  shannon_size_df <- regen_df %>%
    filter(!is.na(as.numeric(Count))) %>%
    group_by(Plot, Subplot, Species, Hclass) %>%
    summarise(abundance = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
    group_by(Plot, Subplot, Hclass) %>%
    mutate(
      abundance_tot = sum(abundance),
      p = (abundance / abundance_tot) * log(abundance / abundance_tot)
    ) %>%
    summarise(H = -sum(p), .groups = "drop") %>%
    pivot_wider(names_from = Hclass, values_from = H, names_prefix = "H_hclass")
  
  growth_size_df <- regen_df %>%
    filter(!is.na(Height_increment_1)) %>%
    mutate(Height_increment = if_else(
      !is.na(Height_increment_2),
      (Height_increment_1 + Height_increment_2) / 2,
      Height_increment_1
    )) %>%
    group_by(Plot, Subplot, Species, Hclass) %>%
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
  
  # ── Final join at Plot/Subplot ────────────────────────────────────────────────
  final_df <- general_df %>%
    left_join(var_df,          by = c("Plot", "Subplot")) %>%
    left_join(shannon_df,      by = c("Plot", "Subplot")) %>%
    left_join(browsing_df,     by = c("Plot", "Subplot")) %>%
    left_join(seedling_df,     by = c("Plot", "Subplot")) %>%
    left_join(shannon_size_df, by = c("Plot", "Subplot")) %>%
    left_join(growth_size_df,  by = c("Plot", "Subplot"))
  return(final_df)
}
