#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_dataformat.R  
#' @description R script containing all functions relative to data
#               processing
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# =============================================================================
# Building the analysis tables
# =============================================================================

#' Bind all tls metrics
build_tls_metrics <- function(metric3d_all, metric3d_all_unsliced, H_ext_all, gini_all) {
  rbind(metric3d_all, metric3d_all_unsliced) %>%
    left_join(H_ext_all, by = c("plot", "subplot", "time")) %>%
    left_join(gini_all,  by = c("plot", "subplot", "time")) %>%
    mutate(lrb_mean = log(rb_mean))
}

#' Bind all forest metrics at the plot level
build_forest_data_plot <- function(regen_metrics_plot, regen_cat, plot_df) {
  regen_metrics_plot$all %>%
    left_join(regen_cat) %>%
    left_join(plot_df) %>%
    mutate(across(c(richness_sapling, n_sapling), ~ tidyr::replace_na(.x, 0)),
           subplot = "plot")
}


#' Bind all forest metrics at the subplot level
build_forest_data_subplot <- function(regen_metrics_subplot, regen_cat_subplot, plot_df) {
  regen_metrics_subplot$all %>%
    left_join(regen_cat_subplot) %>%
    left_join(plot_df) %>%
    unique() %>%
    mutate(across(c(richness_sapling, n_sapling), ~ tidyr::replace_na(.x, 0)))
}

# Table used by the light-effect models (slope sign flipped as in the original)
build_light_effect_data <- function(tls_metrics, forest_data, vars, plot_level) {
  prep_data(tls_metrics, forest_data, vars, plot_level = plot_level, t = 12, s = 0, tg = 0.9) %>%
    mutate(slope_ref_agg = -slope_ref_agg)
}

# Long, scaled table for the light boxplots
light_box_data <- function(tls_metrics, forest_data, vars, plot_level, var_meta) {
  prep_data(tls_metrics, forest_data, c(vars, "reg_category"),
            plot_level = plot_level, t = 12, s = 0) %>%
    tidyr::pivot_longer(-c(plot, subplot, reg_category)) %>%
    filter(!is.na(reg_category)) %>%
    mutate(name = factor(name, levels = vars),
           reg_category = factor(reg_category)) %>%
    group_by(name) %>%
    mutate(value = as.numeric(scale(value))) %>%
    ungroup() %>%
    left_join(var_meta, by = c("name" = "variable")) %>%
    arrange(name) %>%
    mutate(label = forcats::fct_inorder(fix_cross_labels(label, as.character(name))),
           type  = if_else(type == "PAD", "Plant area density\n distribution", type))
}
