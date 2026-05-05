tar_load(dart_folder)
plot_name = "GER11"
hours     = c(8, 10,11, 12,13, 14, 15, 16, 18)


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


# sum_arr <- apply(arr4d_flipped, c(1, 2, 3), sum)
# sd_arr <- apply(arr4d_flipped, c(1, 2, 3), sd)

path_rb=paste0("output/",plot_name,"_rb_raw.rdata")
save(list_rb,file=path_rb)


#get plot time 
plot_name
tar_load(dart_folder)
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

#--- NORMALIzE
plot_name = "GER11"
user="Anne"
path_rb=paste0("output/",plot_name,"_rb_raw.rdata")
bbox=tar_read(bbox_GER11)

load(path_rb)

dtm_rast <- get_dtm(plot_name,user,bbox)


for(r in seq_along(list_rb)){
  rb=list_rb[[r]]
  
  nx <- dim(rb)[2]
  ny <- dim(rb)[1]
  nz <- dim(rb)[3]
  
  
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


load(tar_read(rb_norm_GER11))
i <- 10  # slice index

df <- do.call(rbind, lapply(seq_along(list_rb), function(k) {
  z <- apply(list_rb[[k]][i:(i+10) , ,], c(2,3), sum)
  
  expand.grid(
    x = seq_len(nrow(z)),
    y = seq_len(ncol(z))
  ) |>
    transform(
      value = as.vector(z),
      band = factor(names(list_rb)[k], levels = names(list_rb))
    )
}))

ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  facet_wrap(~ band, ncol = 3) +
  coord_equal() +
  scale_fill_viridis_c(na.value = "transparent") +
  labs(
    title = paste("Slice", i),
    x = NULL,
    y = NULL,
    fill = "Value"
  ) +
  theme_minimal()


## get rb mean, sd for whole plot and subplot
pc_norm=tar_read(pc_norm_GER11)
rb_norm=tar_read(rb_norm_GER11)
subplot_extent=tar_read(subplot_extent_GER11)
bbox=tar_read(bbox_GER11)
plot_name="GER11"
nstrat=4
###

pc_norm<-loadRData(pc_norm)
list_rb<-loadRData(rb_norm)

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
                        rb_mean=numeric(),
                        rb_med=numeric(),
                        rb_sd=numeric(),
                        rb_cv=numeric(),
                        rb_mean_mask=numeric(),
                        rb_med_mask=numeric(),
                        rb_sd_mask=numeric(),
                        rb_cv_mask=numeric(),
                        coef_lin=numeric())


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
                       rb_slice_matched),xy=TRUE) %>% 
      filter(across(everything(),
                    ~!is.na(.))) %>% 
      filter(Z<(slice_up-0.01))
    reglin=lm(rb~Z,data=df)
    coef = reglin$coefficients["Z"][[1]]
    
    
    summary_plot[nrow(summary_plot) + 1, ] <- list(
      plot_name, "plot",time,i, slice_low, slice_up,
      h_mean_na, h_med_na, h_sd_na, h_cv_na,
      h_mean, h_med, h_sd, h_cv,
      rb_mean, rb_med, rb_sd, rb_cv,
      rb_mean_mask, rb_med_mask, rb_sd_mask, rb_cv_mask,
      coef
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
                           rb_sub_slice_matched),xy=TRUE) %>% 
        filter(across(everything(),
                      ~!is.na(.))) %>% 
        filter(Z<(slice_up-0.01))
      reglin_s=lm(rb~Z,data=df_s)
      coef_s = reglin_s$coefficients["Z"][[1]]
      
      
      summary_plot[nrow(summary_plot) + 1, ] <- list(
        plot_name, names(subplot_extent[sp]),time, i, slice_low, slice_up,
        h_mean_s_na, h_med_s_na, h_sd_s_na, h_cv_s_na,
        h_mean_s, h_med_s, h_sd_s, h_cv_s,
        rb_mean_s, rb_med_s, rb_sd_s, rb_cv_s,
        rb_mean_mask_s, rb_med_mask_s, rb_sd_mask_s, rb_cv_mask_s,
        coef_s
      )
    }
  }
}


pc_norm=tar_read(pc_norm_GER11)
rb_norm=tar_read(rb_norm_GER11)
bbox=tar_read(bbox_GER11)
subplot_extent=tar_read(subplot_extent_GER11)
plot_name="GER11"
nstrat=4
