library(targets)
library(tarchetypes)
library(crew)
tar_source(files = "R")
tar_option_set(packages = c("dplyr", "ggplot2","data.table","tidyr","readxl","readr",
                            "ncdf4", "terra","sf",
                            "lidR","ITSMe"),
               controller = crew_controller_local(workers = 6),
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
               get_complexity_rb(pc_norm,
                                 rb_norm,
                                 voxnorm,
                                 bbox,
                                 subplot_extent,
                                 plot_name = plot_name,
                                 nstrat    = 4))
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
             get_regen_metric_subplot(regen_df)),
  tar_target(regen_metrics_plot,
             get_regen_metric_plot(regen_df)),
  # Get inventory data
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
             get_adults(inventory_df))
)