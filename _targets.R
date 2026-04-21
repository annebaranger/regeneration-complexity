library(targets)
# library(tarchetypes)
tar_source(files = "R")
tar_option_set(packages = c("dplyr", "ggplot2","data.table",
                            "ncdf4", "terra",
                            "lidR","ITSMe"))

list(
  tar_target(user, "Anne"),
  tar_target(path_plot,
             paste0("C:/Users/",user,"/OneDrive - University of Cambridge/2. FLF project")),
  tar_target(dart_folder,
             paste0("C:/Users/",user,"/DART-1/user_data/simulations")),
  tar_target(rb,
             get_rb(dart_folder,
                    plot_name="GER02",
                    hours=c(8,10,12,14,16,18)),
             format="file"),
  tar_target(rb_norm,
             normalize_rb(plot_name="GER02",
                          user,
                          rb_path=rb[1],
                          rb_type="sum"))
)