library(targets)
library(tarchetypes)
library(crew)
tar_source(files = "R")
tar_option_set(packages = c("dplyr", "ggplot2","data.table","tidyr","readxl","readr",
                            "ncdf4", "terra","sf","purrr","tibble","mgcv",
                            "lidR","ITSMe"),
               controller = crew_controller_local(workers = 3),
               error = "null")
# lapply(c("targets",
#          "dplyr", "ggplot2","data.table","tidyr","readxl",
#          "ncdf4", "terra","sf",
#          "lidR","ITSMe"),require,character.only=TRUE)

plot_values <- tibble::tibble(
  plot_name = c("GER02", "GER01","GER06","GER09","GER11","GER12","GER13","GER14",
                "GER18","GER19","GER20","GER21","GER27","GER29","GER34","GER35",
                "GER37","GER38","GER32")   # <-- add all your plots here
)

# =============================================================================
# Variable sets + targets for the paper figures
# Add `paper_targets` to the list returned at the end of _targets.R
# (targets tracks `paper_vars`: editing it invalidates the dependent targets)
# =============================================================================

paper_vars <- list(
  pca_regen  = c("richness", "H", "richness_sapling", "cwm_shade", "abundance", "n_sapling",
                 "Height_increment_mn", "H_height", "empty_gf", "pad_tot", "pad_cv", "ri", "fhd"),
  lidar      = c("empty_gf", "pad_tot", "pad_cv", "ri", "fhd"),
  light      = c("lrb_mean", "ext_ref_agg", "rb_cv", "gini_h_rb", "bl_diff_abs", "slope_ref_agg"),
  light_box  = c("lrb_mean", "ext_ref_agg", "rb_cv", "gini_h_rb", "bl_diff_abs", "slope_ref_agg",
                 "padAbove_agg", "padAbove_cv"),
  regen      = c("richness", "richness_sapling", "H", "cwm_shade", "abundance", "n_sapling",
                 "Height_increment_mn", "empty_gf", "H_height", "fhd", "ri", "pad_cv"),
  light_expl = c("lrb_mean", "gini_h_rb", "slope_ref_agg", "bl_diff_abs"),
  adult_expl = c("n_trees_alive", "species_richness_alive", "mean_dbh_alive_cm",
                 "basal_area_alive_m2", "mean_height_alive_m", "browsing"),
  variability = c("ext_ref_agg", "slope_ref_agg", "gini_h_rb"),
  vp_adult   = c("species_richness_alive", "mean_dbh_alive_cm", "mean_height_alive_m"),
  vp_light   = c("lrb_mean", "gini_h_rb", "slope_ref_agg"),
  vp_regen   = c("richness", "abundance", "ri"),
  vp_log     = c("abundance"),
  vp_labels  = c(species_richness_alive = "Adult species richness",
                 mean_dbh_alive_cm = "Mean DBH", mean_height_alive_m = "Mean height",
                 lrb_mean = "Mean light", gini_h_rb = "Horizontal heterogeneity",
                 slope_ref_agg = "Vertical heterogeneity",
                 richness = "Richness", abundance = "Abundance", ri = "Complexity")
)


# Mapped targets (run once per plot) ────────────────────────────────────────
mapped<-tar_map(
    values = plot_values,      # iterates over plot_name
    names  = plot_name,        # used as suffix in target names e.g. rb_GER02

    tar_target(bbox,
               get_bbox(plot_name=plot_name,
                        user)),
    tar_target(dtm,
               get_dtm(plot_name=plot_name,
                       user,
                       bbox),
               format="file"),
    tar_target(simu_time,
               get_simulation_time(dart_folder,
                                   plot_name)),
    tar_target(rb,
               get_rb(dart_folder,
                      plot_name = plot_name,
                      hours     = simu_time),
               format = "file"),
    
    tar_target(rb_norm,
               normalize_rb(plot_name = plot_name,
                            user,
                            rb_path = rb,
                            bbox,
                            dtm),
               format = "file"),
    
    tar_target(voxnorm,
               get_pad(user,
                       path_plot,
                       plot_name),
               format = "file"),
    
    tar_target(pc_norm,
               get_pc_norm_clean(plot_name = plot_name,
                                 user,
                                 bbox,
                                 dtm),
               format = "file"),
    tar_target(subplot_extent,
               get_subplot_extent(plot_name = plot_name,
                                  user)),
    
    tar_target(metric3d_df,
               get_complexity(pc_norm,
                              rb_norm,
                              voxnorm,
                              bbox,
                              subplot_extent,
                              plot_name = plot_name,
                              sliced=TRUE,
                              nstrat    = 4)),
    tar_target(metric3d_df_unsliced,
               get_complexity(pc_norm,
                              rb_norm,
                              voxnorm,
                              bbox,
                              subplot_extent,
                              plot_name = plot_name,
                              sliced=FALSE,
                              slice_low=0,
                              slice_up=4)),
    tar_target(H_ext,
               get_Hext(rb_norm,
                        pc_norm,
                        voxnorm,
                        subplot_extent,
                        bbox,
                        plot_name,
                        targets        = c(0.5, 0.9),
                        z_ref          = 4.5,
                        voxel_size     = 0.25,
                        pad_unit       = c("density"),
                        smooth_columns = TRUE,
                        touches        = FALSE,
                        return_rasters = FALSE)),
    tar_target(gini,
               get_light_metrics(rb_norm,
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
                                 pad_is_density = TRUE)),
    tar_target(file_segmented,
               file.segmented(path_plot,
                              plot_name,
                              folders=c("regen","hesitation"))),
    tar_target(pc_norm_seg,
               get_pc_norm_regen(path_plot,
                                 plot_name,
                                 file_segmented,
                                 dtm,
                                 bbox)),
    tar_target(
      metric3d_df_reg, {
        targets::tar_cancel(is.null(pc_norm_seg))
        get_below_canopy_complexity(
          pc_norm_seg, rb_norm, voxnorm, bbox,
          subplot_extent, plot_name,
          voxel_size = 0.25, chm_res = 0.1
        )
      }
    )
    )


list(
  # Static targets (run once) ─────────────────────────────────────────────────
  tar_target(user, "Anne"),
  tar_target(path_plot,
             paste0("C:/Users/", user, "/OneDrive - University of Cambridge/2. FLF project")),
  tar_target(dart_folder,
             paste0("C:/Users/", user, "/DART-1/user_data/simulations")),
  # Combine all plots into one df at the end ─────────────────────────────────
  mapped, 
  tar_combine(
    metric3d_all,
    mapped[["metric3d_df"]],   # grabs metric3d_df_GER02, _GER20, etc.
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    metric3d_all_unsliced,
    mapped[["metric3d_df_unsliced"]],   # grabs metric3d_df_GER02, _GER20, etc.
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    metric3d_all_seg,
    mapped[["metric3d_df_reg"]],   # grabs metric3d_df_GER02, _GER20, etc.
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    H_ext_all,
    mapped[["H_ext"]],   # grabs metric3d_df_GER02, _GER20, etc.
    command = dplyr::bind_rows(!!!.x)
  ),
  tar_combine(
    gini_all,
    mapped[["gini"]],   # grabs metric3d_df_GER02, _GER20, etc.
    command = dplyr::bind_rows(!!!.x)
  ),
  # Analyse regeneration survey ─────────────────────────────────────────────── 
  tar_target(regen_path,
             "//ifs-prod-596-cifs.ifs.uis.private.cam.ac.uk/geog-forest/germany-2025/germany-regen/Regen_Fundiv_Germany_0825.xlsx",
             format="file"),
  tar_target(regen_df,
             readxl::read_excel(regen_path,sheet=3) %>%
               dplyr::mutate(dplyr::across(6:13, as.numeric)) %>% 
               rename(plot=Plot,
                      subplot=Subplot) %>% 
               mutate(subplot=paste0("subplot_",subplot))),
  tar_target(regen_metrics_subplot,
             get_regen_metric_subplot(regen_df,shade_df)),
  tar_target(regen_metrics_plot,
             get_regen_metric_plot(regen_df,shade_df)),
  tar_target(regen_cat,
             readxl::read_excel(regen_path,sheet=2) %>%
               dplyr::mutate(dplyr::across(c(2,4:12), as.numeric)) %>% 
               rename(plot=Plot,
                      subplot=Subplot) %>% 
               mutate(subplot=paste0("subplot_",subplot)) %>% 
               select(-c(14,13)) %>% 
               filter(!is.na(Regeneration)) %>% 
               group_by(plot) %>% 
               summarize(mean_reg=mean(Regeneration),
                         .groups = "drop"
               ) %>% 
               dplyr::mutate(
                 reg_category = dplyr::ntile(mean_reg, 3),
                 reg_category = factor(
                   reg_category,
                   levels = 1:3,
                   labels = c("Low", "Medium", "High")
                 )
               )
  ),
  tar_target(regen_cat_subplot,
             readxl::read_excel(regen_path,sheet=2) %>%
               dplyr::mutate(dplyr::across(c(2,4:12), as.numeric)) %>% 
               rename(plot=Plot,
                      subplot=Subplot) %>% 
               mutate(subplot=paste0("subplot_",subplot)) %>% 
               select(-c(14,13)) %>% 
               filter(!is.na(Regeneration)) %>% 
               dplyr::mutate(
                 reg_category = dplyr::ntile(Regeneration, 3),
                 reg_category = factor(
                   reg_category,
                   levels = 1:3,
                   labels = c("Low", "Medium", "High")
                 )
               )
  ),
  
  # Get inventory data ─────────────────────────────────────────────────────────
  tar_target(inventory_path,
             "C:/Users/Anne/OneDrive - University of Cambridge/2. FLF project/germany-2025/regeneration-germany/Plot_descriptors_tree_data_Year2017_Inventory2_Germany.xls",
             format="file"
  ),
  tar_target(inventory_df,
             readxl::read_excel(inventory_path,sheet = "Raw data") %>% 
               filter(PlotID %in% plot_values$plot_name) %>% 
               rename(plot=PlotID) %>% 
               dplyr::mutate(dplyr::across(c(2,3,6:15), as.numeric))
               
  ),
  tar_target(plot_df,
             get_adults(inventory_df)),
  
  
  # get shade tolerance ─────────────────────────────────────────────────────────
  tar_target(
    shade_df,
    read.csv("shadetolerance.csv") %>% 
      separate(Species,into=c("Genus","species")) %>% 
      filter(Data.set.1=="Europe") %>% 
      mutate(species=paste0(toupper(substr(Genus,1,2)),toupper(substr(species,1,2)))) %>% 
      filter(species %in% regen_df$Species) %>% 
      select(species,shade_tolerance=Shade.tolerance) %>% 
      bind_rows(data.frame(species="CRSP",shade_tolerance=1.93))
  ),
  
  # --- format data ------------------------------------------------------------
  tar_target(var_meta, make_var_meta()),
  tar_target(tls_metrics,
             build_tls_metrics(metric3d_all, metric3d_all_unsliced, H_ext_all, gini_all)),
  tar_target(forest_data_plot,
             build_forest_data_plot(regen_metrics_plot, regen_cat, plot_df)),
  tar_target(forest_data_subplot,
             build_forest_data_subplot(regen_metrics_subplot, regen_cat_subplot, plot_df)),
  
  # --- regeneration structure PCA ------------------------------------------
  tar_target(pca_regen,
             run_pca_regen(tls_metrics, forest_data_plot, forest_data_subplot, regen_cat,
                           paper_vars$pca_regen, paper_vars$lidar, var_meta)),
  tar_target(cor_regen, cor_table(pca_regen$data)),
  tar_target(fig_pca_regen, plot_pca_regen(pca_regen)),
  tar_target(file_fig1,
             save_fig(fig_pca_regen, "mainfig1_PCA_regen", 8.5, 11, scale = 1.2),
             format = "file"),
  
  # --- 2. light descriptors, plot level ---------------------------------------
  tar_target(pca_light_plot,
             run_pca_light(tls_metrics, forest_data_plot, paper_vars$light, "plot", var_meta)),
  tar_target(cor_light, cor_table(pca_light_plot$data)),
  tar_target(box_light_plot,
             light_box_data(tls_metrics, forest_data_plot, paper_vars$light_box, "plot", var_meta)),
  tar_target(fig_light_plot,
             plot_light_panel(plot_pca_light(pca_light_plot),
                              plot_light_boxes(box_light_plot, kw = kw_by_variable(box_light_plot)))),
  tar_target(file_fig2,
             save_fig(plot_pca_light(pca_light_plot), "mainfig2_PCA_light", 8.5, 11, scale = 1.1),
             format = "file"),
  tar_target(file_fig2.1,
             save_fig( plot_light_boxes(box_light_plot, kw = kw_by_variable(box_light_plot)), 
                       "si_lightdist", 16, 9, scale = 1.1),
             format = "file"),
  
  # --- 3. light effects, plot level -------------------------------------------
  tar_target(light_effect_plot,
             build_light_effect_data(tls_metrics, forest_data_plot,
                                     c(paper_vars$regen, paper_vars$light_expl), "plot")),
  tar_target(boot_light_plot,
             boot_effects(light_effect_plot, paper_vars$regen, paper_vars$light_expl, B = 1000)),
  tar_target(fig_light_eff,
             plot_effects(boot_light_plot, var_meta, paper_vars$regen,
                          paper_vars$light_expl, drop = "bl_diff_abs")),
  tar_target(file_fig3,
             save_fig(fig_light_eff, "mainfig3_light_eff", 16, 9, scale = 1.2),
             format = "file"),
  
  # --- 4. subplot level --------------------------------------------------------
  tar_target(pca_light_subplot,
             run_pca_light(tls_metrics, forest_data_subplot, paper_vars$light, "subplot", var_meta)),
  tar_target(box_light_subplot,
             light_box_data(tls_metrics, forest_data_subplot, paper_vars$light_box, "subplot", var_meta)),
  tar_target(fig_light_subplot,
             plot_light_panel(plot_pca_light(pca_light_subplot, point_size = 1.5),
                              plot_light_boxes(box_light_subplot,
                                               cld = cld_by_variable(box_light_subplot)))),
  # tar_target(light_effect_subplot,
  #            build_light_effect_data(tls_metrics, forest_data_subplot,
  #                                    c(paper_vars$regen, paper_vars$light_expl), "subplot")),
  # tar_target(boot_light_subplot,
  #            boot_effects(light_effect_subplot, paper_vars$regen, paper_vars$light_expl,
  #                         B = 1000, random = "plot", cluster_boot = FALSE)),
  # tar_target(fig_light_eff_subplot,
  #            plot_effects(boot_light_subplot, var_meta, paper_vars$regen,
  #                         paper_vars$light_expl, drop = "bl_diff_abs", lw = c(0.9, 1.5))),
  
  # --- 5. temporal vs spatial variability ---------------------------------------
  tar_target(light_variability, compute_light_variability(tls_metrics, paper_vars$variability)),
  tar_target(fig_variability, plot_variability(light_variability, var_meta)),
  tar_target(file_fig5.1,
             save_fig(fig_variability, 
                       "si_tempspatvar", 16, 9, scale = 1.1),
             format = "file"),
  
  # --- 6. adult structure effects -----------------------------------------------
  tar_target(adult_effect,
             prep_data(tls_metrics, forest_data_plot,
                       id_cols = c(paper_vars$regen, paper_vars$adult_expl), plot_level = "plot")),
  tar_target(boot_adult,
             boot_effects(adult_effect, paper_vars$regen, paper_vars$adult_expl, B = 1000)),
  tar_target(fig_adult_eff,
             plot_effects(boot_adult, var_meta, paper_vars$regen,
                          paper_vars$adult_expl, lw = c(0.9, 1.5))),
  tar_target(file_fig6,
             save_fig(fig_light_eff, "si_adult_eff", 16, 9, scale = 1.2),
             format = "file"),
  
  # --- 7. variance partitioning --------------------------------------------------
  tar_target(varpart,
             run_varpart(tls_metrics, forest_data_plot,
                         adult_vars = paper_vars$vp_adult, light_vars = paper_vars$vp_light,
                         regen_vars = paper_vars$vp_regen, log_vars = paper_vars$vp_log,
                         n_perm = 999, n_perm_hp = 999)),
  tar_target(fig_varpart_venn, plot_varpart_venn(varpart)),
  tar_target(file_fig7.1,
             save_fig(fig_varpart_venn, "si_venndiag", 8.5, 6, scale = 1.2),
             format = "file"),
  tar_target(fig_varpart_uni,  plot_varpart_univariate(varpart, paper_vars$vp_labels)),
  tar_target(file_fig7,
             save_fig(fig_varpart_uni, "mainfig7_varpart", 16, 9, scale = 1.2),
             format = "file"),
  tar_target(fig_rda_triplot,  plot_rda_triplot(varpart, paper_vars$vp_labels)),
  tar_target(file_fig7.2,
             save_fig(fig_rda_triplot, "si_triplot", 8.5, 9, scale = 1.2),
             format = "file")
  )

