lapply(c("dplyr", "ggplot2","data.table","tidyr","readxl","readr",
         "ncdf4", "terra","sf",
         "lidR","ITSMe"),require,character.only=TRUE)
library(targets)

pc_norm=tar_read(pc_norm_GER11)
rb_norm=tar_read(rb_norm_GER11)
voxnorm=tar_read(voxnorm_GER11)
bbox=tar_read(bbox_GER11)
subplot_extent=tar_read(subplot_extent_GER11)
plot_name="GER11"
targets    = c(0.5, 0.9)
z_ref      = 4.5
voxel_size = 0.25
pad_unit   = c("density", "index")
pad_unit="density"
# nstrat=4


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



z_cross=prof %>% 
  filter(i==100,j==100)  %>%
  group_by(i, j) %>%
  arrange(desc(z), .by_group = TRUE) %>%          # top of canopy -> down
  mutate(prev_ext = lag(extinction),
         prev_z   = lag(z)) %>%
  expand_grid(ext=c(.5,.9)) %>% 
  filter(extinction >= ext, prev_ext < ext) %>%   # the crossing pair
  group_by(i,j,ext) %>% 
  slice(1) %>%                                    # first crossing only
  mutate(z_cross = prev_z + (ext - prev_ext) *
           (z - prev_z) / (extinction - prev_ext)) %>%
  select(i, j,ext, z_cross) %>%
  ungroup()

extRef=prof %>%
  filter(i==100,j==100) %>% 
  filter(z==4.5) %>% 
  select(i,j,rb_ref=value,ext_ref=extinction)


prof %>%
  left_join(z_cross) %>% 
  group_by(i,j,ext) %>% 
  group_modify(~{
    m_z <- mgcv::gam(extinction ~ s(z, k = 5), data = .x)
    dz  <- voxel_size
    
    zc <- .x$z_cross[1]

    e <- as.numeric(predict(m_z, newdata = tibble(
      z = c(zc + dz, zc - dz, z_ref + dz, z_ref - dz))))
    
    tibble(slope_cross = (e[1] - e[2]) / (2 * dz),
           slope_ref   = (e[3] - e[4]) / (2 * dz))
  })%>% 
  left_join(extRef) %>% 
  left_join(z_cross)


data.frame(ext= 1 - rb[1,1,] / rb_max,
           x=1:length(rb[1,1,])) %>% 
  ggplot(aes(x,ext))+geom_point()+geom_hline(yintercept=c(0.9,0.5))


# get extinction characteristics
#' @param voxnorm path to the voxelised, height-normalised PAD array (i, j, z),
#'   same grid as rb_norm
#' @param plot_name plot name
#' @param targets characteristic extinction levels at which the height of the
#'   profile is computed
#' @param z_ref reference height (m) at which extinction, its vertical slope and
#'   the PAD above are evaluated
#' @param voxel_size vertical voxel size (m); also the finite-difference step for
#'   the vertical slope at z_ref
#' @param pad_unit "density" if voxpad stores PAD in m2/m3 (values are multiplied
#'   by voxel_size to get a plant area index), "index" if each voxel already
#'   holds m2/m2. Only affects the magnitude of padAbove / padTot, not zPAD50.
get_Hext <- function(rb_norm,
                     pc_norm,
                     voxnorm,
                     subplot_extent,
                     bbox,
                     plot_name,
                     targets    = c(0.5, 0.88),
                     z_ref      = 4.5,
                     voxel_size = 0.25,
                     pad_unit   = c("density", "index")) {
  
  pad_unit <- match.arg(pad_unit)
  
  list_rb <- loadRData(rb_norm)         # load radiative budget
  pc_norm <- loadRData(pc_norm)
  voxpad  <- loadRData(voxnorm)         # voxelised PAD
  
  if (!identical(dim(voxpad)[1:2], dim(list_rb[[1]])[1:2]))
    stop("voxpad and rb_norm do not share the same horizontal grid")
  
  ## ------------------------------------------------------------------------
  ## PAD metrics, computed once (voxpad does not depend on simulation time).
  ## Convention: voxel k spans [(k-1)*voxel_size, k*voxel_size], i.e. the same
  ## indexing as z = z * voxel_size used on the radiative budget below.
  ##   padAbove = plant area index located above z_ref
  ##   padTot   = total plant area index of the column
  ##   zPAD50   = height below which 50 % of the PAD is found (median height of
  ##              the vertical PAD distribution)
  ## ------------------------------------------------------------------------
  nz    <- dim(voxpad)[3]
  z_top <- seq_len(nz) * voxel_size
  z_bot <- z_top - voxel_size
  pai   <- voxpad * if (pad_unit == "density") voxel_size else 1
  
  # fraction of each voxel lying above z_ref (handles a z_ref inside a voxel)
  w_above <- pmin(pmax((z_top - z_ref) / voxel_size, 0), 1) # which voxels are above z_ref
  
  
  # qunatify variation
  pad_above <- apply(pai, 1:2, function(v) sum(v * w_above, na.rm = TRUE))
  pad_tot   <- apply(pai, 1:2, sum, na.rm = TRUE)
  z_pad50   <- apply(pai, 1:2, function(v) {
    v   <- ifelse(is.na(v), 0, v)
    tot <- sum(v)
    if (tot <= 0) return(NA_real_)
    cs   <- cumsum(v)
    k    <- which(cs >= 0.5 * tot)[1]
    prev <- if (k == 1L) 0 else cs[k - 1L]
    z_bot[k] + (0.5 * tot - prev) / v[k] * voxel_size   # linear within the voxel
  })
  
  pad_tab <- tibble(i        = as.vector(row(pad_above)),
                    j        = as.vector(col(pad_above)),
                    padAbove = as.vector(pad_above),
                    padTot   = as.vector(pad_tot),
                    zPAD50   = as.vector(z_pad50)) %>% 
    summarize(padAbove_mean=mean(padAbove),
            padTot_mean=mean(padTot),
            zPAD50_mean=mean(zPAD50,na.rm=TRUE),
            padAbove_sd=sd(padAbove),
            padTot_sd=sd(padTot),
            zPAD50_sd=sd(zPAD50,na.rm=TRUE),
            padAbove_cv=padAbove_sd/padAbove_mean,
            padTot_cv=padTot_sd/padTot_mean,
            zPAD50_cv=zPAD50_sd/zPAD50_mean) %>% 
    mutate(plot=plot_name,subplot="plot")
  
  #quantify mean
  pad_1d<-apply(pai, 3, sum, na.rm = TRUE)
  cs_1d   <- cumsum(pad_1d)
  pad_tot_plot=sum(pad_1d)
  pad_above_plot <- sum(pad_1d * w_above, na.rm = TRUE)
  k    <- which(cs_1d >= 0.5 * pad_tot_plot)[1]
  prev <- if (k == 1L) 0 else cs_1d[k - 1L]
  z_pad50_plot=z_bot[k] + (0.5 * pad_tot_plot - prev) / pad_1d[k] * voxel_size
  
  pad_tab$z_pad50_plot=z_pad50_plot
  pad_tab$pad_tot_plot=pad_tot_plot
  pad_tab$pad_above_plot=pad_above_plot
  
  ## get a raster --------------------------------------------------------
  rast_pai=rast(pai)
  names(rast_pai)=1:nz
  
  summary <- data.frame()
  
  for (t in seq_along(list_rb)) {      # loop over simulation time
    rb_norm <- list_rb[[t]]
    time    <- as.numeric(names(list_rb)[t])
    
    prof <- reshape2::melt(rb_norm,
                           varnames   = c("i", "j", "z"),
                           value.name = "value") %>%
      mutate(max_val=max(value)) %>% 
      group_by(i, j) %>%
      mutate(extinction = 1 - value / max_val,
             z          = z * voxel_size)

    ## per-plot metrics --------------------------------------------
    nz_rb    <- dim(rb_norm)[3]
    z_top_rb <- seq_len(nz_rb) * voxel_size
    z_bot_rb <- z_top_rb - voxel_size
    
    rb_norm_1d<-apply(rb_norm, 3, sum, na.rm = TRUE)
    ext_1d<-1 - rb_norm_1d / max(rb_norm_1d)
    
    get_hx<-function(ext_1d,th){
      kx<- max(which(ext_1d >= th ))
      prev <- if (kx == 1L) 0 else ext_1d[kx - 1L]
      aft <- if (kx == 1L) 0 else ext_1d[kx + 1L]
      return(c(z_bot_rb[kx],(aft - prev)/2 ))
    }
    k_50    <- get_hx(ext_1d,th=0.5)[1]
    slope_50    <- get_hx(ext_1d,th=0.5)[2] 
    k_90    <- get_hx(ext_1d,th=0.9)[1]
    slope_90    <- get_hx(ext_1d,th=0.9)[2] 
    extRef_plot =ext_1d[4.5/0.25]
    rbref_plot=rb_norm_1d[4.5/0.25]
    slope_ref=(ext_1d[4.5/0.25 + 1L]-ext_1d[4.5/0.25 - 1L])/2
    
        
    ## per-column (i,j) metrics --------------------------------------------
    
    # get height at which target extinction is reached
    z_cross=prof %>% 
      group_by(i, j) %>%
      arrange(desc(z), .by_group = TRUE) %>%          # top of canopy -> down
      mutate(prev_ext = lag(extinction),
             prev_z   = lag(z)) %>%
      expand_grid(ext=c(.5,.9)) %>% 
      filter(extinction >= ext, prev_ext < ext) %>%   # the crossing pair
      group_by(i,j,ext) %>% 
      slice(1) %>%                                    # first crossing only
      mutate(z_cross = prev_z + (ext - prev_ext) *
               (z - prev_z) / (extinction - prev_ext)) %>%
      select(i, j,ext, z_cross) %>%
      ungroup()
    
    # get extinction at targetted height
    extRef=prof %>%
      filter(z==4.5) %>% 
      select(i,j,rb_ref=value,ext_ref=extinction)
    
    # get slope at those characteristic heights
    half<- prof %>%
      left_join(z_cross) %>% 
      group_by(i,j,ext) %>% 
      group_modify(~{
        m_z <- mgcv::gam(extinction ~ s(z, k = 5), data = .x)
        dz  <- voxel_size
        
        zc <- .x$z_cross[1]
        
        e <- as.numeric(predict(m_z, newdata = tibble(
          z = c(zc + dz, zc - dz, z_ref + dz, z_ref - dz))))
        
        tibble(slope_cross = (e[1] - e[2]) / (2 * dz),
               slope_ref   = (e[3] - e[4]) / (2 * dz))
      })%>% 
      left_join(extRef) %>% 
      left_join(z_cross) %>% 
      ungroup()
    
    
    

    
    ## one value per column: time-varying light metrics + static PAD metrics
    summary <- summary %>% 
      bind_rows(half %>%
                  group_by(ext) %>% 
                  summarise(across(c(slope_cross, slope_ref, z_cross, ext_ref, rb_ref),
                                   list(mean = ~mean(.x, na.rm = TRUE),
                                        sd   = ~sd(.x,   na.rm = TRUE)),
                                   .names = "{.col}_{.fn}")) %>%
                  mutate(plot=plot_name,subplot="plot") %>% 
                  ungroup() %>% 
                  left_join(data.frame(ext=c(0.5,0.9),
                                       z_cross_plot=c(k_50,k_90),
                                       slope_cross_plot=c(slope_50,slope_90),
                                       rb_ref_plot=rbref_plot,
                                       ext_ref_plot=extRef_plot,
                                       slope_ref_plot=slope_ref))
      )
    

    ## get a raster --------------------------------------------------------
    wide <- half %>%
      pivot_wider(id_cols     = c(i, j),
                  names_from  = ext,
                  values_from = c(slope_cross, z_cross)) %>%
      left_join(half %>% select(i,j,slope_ref, ext_ref, rb_ref), by = c("i", "j")) %>% 
      rename(x=j,y=i) %>% relocate(x, .before = y)# adds the 5 scalar layers
    
    stk <- terra::rast(wide, type = "xyz")
    
    # prepare reference raster
    pc_chm <- rasterize_canopy(pc_norm, res = 0.15, algorithm = p2r())
    ext(stk) <- ext(pc_chm)
    crs(stk) <- crs(pc_chm)
    stk_matched <- resample(stk, pc_chm, method = "bilinear")
    
    
    #### HERE ###
    for (sp in seq_along(subplot_extent)) {
      plot_ext <- subplot_extent[[sp]]
      plot_ext_trans <- shift(vect(plot_ext),
                              dx = -bbox$xmin,
                              dy = -bbox$ymin)
      
      stk_sub <- crop(stk_matched, plot_ext_trans)
      df_sub  <- as.data.frame(stk_sub, xy = TRUE)
      
      ## target-based metrics (long format, as before)
      targ_sub <- df_sub %>%
        # select(x, y, starts_with("z_pred_"), starts_with("slope_")) %>%
        pivot_longer(
          cols      = -c(x, y,slope_ref, ext_ref, rb_ref),
          names_to  = c(".value", "extinction"),
          names_sep = "_(?=[0-9])"        # split at the underscore before the number
        ) %>%
        group_by(extinction) %>%
        summarise(mean       = mean(z_pred),
                  sd         = sd(z_pred),
                  mean_slope = mean(slope),
                  sd_slope   = sd(slope),
                  .groups = "drop")
      
      ## height-based metrics (one value per pixel)
      scal_sub <- df_sub %>%
        summarise(across(c(extRef, slopeRef, padAbove, padTot, zPAD50),
                         list(mean = ~mean(.x, na.rm = TRUE),
                              sd   = ~sd(.x,   na.rm = TRUE)),
                         .names = "{.fn}_{.col}"))
      
      summary_height <- summary_height %>%
        bind_rows(targ_sub %>%
                    mutate(!!!as.list(scal_sub),
                           extinction = as.numeric(extinction),
                           time       = time,
                           plot       = plot_name,
                           subplot    = names(subplot_extent[sp])))
    }
  }
  return(summary_height)
}