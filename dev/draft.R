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


# ── 1. Données longues ────────────────────────────────────────────────────────
data <- regen_metrics_subplot$all %>%
  select(plot, subplot,
         richness, richness_sd, richness_seedling, richness_sapling,
         H, H_hclass1, H_hclass2, H_hclass3, H_hclass4, H_hclass5,
         H_height) %>%
  left_join(
    metric3d_all %>%
      select(plot, subplot, slice_num, time, rb_sd, rb_mean, rb_cv, pad_sd) %>%
      filter(time == 12, subplot != "plot"),
    by = c("plot", "subplot")
  ) %>%
  filter(plot %in% metric3d_all$plot) %>%
  mutate(
    slice_num = factor(slice_num, levels = sort(unique(slice_num))),
    plot      = as.factor(plot)
  ) %>%
  pivot_longer(
    cols = -c(plot, subplot, slice_num, time, rb_sd, rb_mean, rb_cv),
    names_to  = "y",
    values_to = "y_val"
  ) %>%
  group_by(y) %>%
  mutate(y_val = as.numeric(scale(y_val))) %>%   # as.numeric : évite la colonne-matrice
  ungroup()

vars <- c("Total", "richness_sd", "richness_seedling", "richness_sapling",
          "H", "H_hclass1", "H_hclass2", "H_hclass3", "H_hclass4", "H_hclass5",
          "H_height", "pad_sd")

# ── 2. Bootstrap (bootMer) par variable ──────────────────────────────────────
# Modèle mixte : y_val ~ slice_num + slice_num:rb_cv + (1 | plot)
#   - slice_num discret (facteur)
#   - chaque coef "slice_num<niveau>:rb_cv" = pente de rb_cv DANS ce niveau
#   - (1 | plot) = intercept aléatoire par plot

slope_fun <- function(m) {
  fx <- fixef(m)
  fx[grepl(":rb_cv$", names(fx))]   # une pente par niveau de slice_num
}

boot_results <- map_dfr(vars, function(v) {
  
  data_lm <- data %>%
    filter(y == v) %>%
    drop_na(y_val, rb_cv, slice_num) %>%
    droplevels()
  
  if (dplyr::n_distinct(data_lm$plot) < 5) return(NULL)   # assez de plots pour l'effet aléatoire
  
  fit <- suppressMessages(suppressWarnings(
    lmer(y_val ~ slice_num + slice_num:rb_cv + (1 | plot),
         data = data_lm, REML = FALSE,
         control = lmerControl(check.conv.singular = "ignore"))
  ))
  
  bm <- suppressWarnings(
    bootMer(fit, slope_fun, nsim = 1000,
            type = "parametric", use.u = FALSE, .progress = "none")
  )
  
  cn   <- colnames(bm$t)                                   # noms des pentes
  lvls <- sub("^slice_num", "", sub(":rb_cv$", "", cn))    # niveau de slice_num
  
  map_dfr(seq_along(cn), function(i) {
    tibble(
      var       = v,
      slice_num = lvls[i],
      slope     = bm$t[, i]
    )
  })
})

# ── 3. Plot ridgeline ─────────────────────────────────────────────────────────
slice_levels <- levels(data$slice_num)
n_vars       <- length(vars)
n_slices     <- length(slice_levels)

boot_results %>%
  mutate(
    slice_num  = factor(slice_num, levels = slice_levels),
    var_idx    = as.numeric(factor(var, levels = vars)),
    slice_code = as.numeric(slice_num),                    # 1..n_slices
    slice_pos  = var_idx +
      (slice_code - mean(seq_len(n_slices))) / max(n_slices - 1, 1) * 0.8
  ) %>%
  ggplot(aes(x = slope, y = slice_pos, group = interaction(var, slice_num),
             fill = slice_num)) +
  stat_density_ridges(
    quantile_lines = TRUE, quantiles = c(0.05, 0.5, 0.95),
    alpha = 0.6, color = "white", linewidth = 0.3,
    scale = 2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(breaks = 1:n_vars, labels = vars) +
  scale_fill_viridis_d(option = "mako", name = "Slice height") +   # _d : facteur
  labs(x = "Slope of rb_cv (bootstrap)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(face = "bold")
  )


boot_results %>%
  mutate(
    slice_num  = factor(slice_num, levels = slice_levels),
    var_idx    = as.numeric(factor(var, levels = vars)),
    slice_code = as.numeric(slice_num),                    # 1..n_slices
    slice_pos  = var_idx +
      (slice_code - mean(seq_len(n_slices))) / max(n_slices - 1, 1) * 0.8
  ) %>%
  filter(var %in% c("H", "H_hclass1", "H_hclass2", "H_hclass3", "H_hclass4", "H_hclass5")) %>%
  # filter(var %in% c("H_height", "pad_sd")) %>%
  # mutate(var=case_when(var=="H_height"~"Height heterogeneity",
  #                      var=="pad_sd"~"Biomass distribution\n heterogeneity"
  # )) %>% 
  mutate(var=case_when(var=="H_hclass3"~"H (Height class 3)",
                       var=="H_hclass1"~"H (Height class 1)",
                       var=="H_hclass2"~"H (Height class 2)",
                       var=="H_hclass4"~"H (Height class 4)",
                       var=="H_hclass5"~"H (Height class 5)",
                       var=="H"~"H"
                       )) %>%
  # mutate(var=case_when(var=="richness"~"Total",
  #                      var=="richness_seedling"~"Seedlings",
  #                      var=="richness_sapling"~"Saplings",
  #                      var=="richness_sd"~"SD"
  #                      ),
  #        var=factor(var,
  #                   levels=c("Total","Seedlings","Saplings","SD"))) %>% 
  
  ggplot(aes(x = slope, y = var, group = interaction(var, slice_num),
             fill = slice_num)) +
  stat_density_ridges(
    quantile_lines = TRUE, quantiles = c(0.05, 0.5, 0.95),
    alpha = 0.6, color = "white", linewidth = 0.3,
    scale = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_viridis_d(option = "mako", name = "Slice height") +   # _d : facteur
  labs(x = "Slope of Y~(light CV)", y = NULL,
       title="Sublot level") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(face = "bold")
  )


regen_metrics_subplot$all %>% 
  select(plot,subplot,richness,abundance,H_height,n_sapling) %>% 
  left_join(metric3d_all %>% 
              filter(subplot!="plot") %>% 
              filter(time==12) %>% 
              select(plot,subplot,slice_num,h_sd,occupied_gf,pad_sd)) ->t
  
t %>% 
  ggplot(aes(occupied_gf,richness))+
  geom_point()+
  facet_wrap(~slice_num)
  