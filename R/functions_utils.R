#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_utils.R  
#' @description R script containing helpers functions
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# --- Variable metadata (labels, abbreviations, types) -------------------------

#' Create a correspondance table between vars names and figures names

make_var_meta <- function() {
  tibble::tribble(
    ~variable,                ~type,                            ~label,                                 ~abrv,
    "richness",               "Composition",                    "Richness \ntotal",                     "Richn. tot.",
    "richness_sapling",       "Composition",                    "Richness \nsaplings",                  "Richn. sap.",
    "H",                      "Composition",                    "Shannon \nindex",                      "H",
    "cwm_shade",              "Composition",                    "Shade tolerance \n (CWM)",             "Shade tol.",
    "abundance",              "Quantity",                       "Abundance",                            "Ab.",
    "abundance_cv",           "Quantity",                       "Abundance CV",                         "Ab. (cv)",
    "n_sapling",              "Quantity",                       "Abundance of \nsaplings",              "Ab. sap.",
    "h_mean",                 "Quantity",                       "Mean \nheight",                        "h mean",
    "pad_tot",                "Quantity",                       "PAD \ntotal",                          "PAD tot.",
    "empty_gf",               "Quantity",                       "Gap fraction",                         "Gap frac.",
    "Height_increment_mn",    "Quantity",                       "Height \nincrement",                   "h incr.",
    "H_height",               "Complexity",                     "Vertical \ndistribution",              "Vert. dist.",
    "h_cv",                   "Complexity",                     "Height \nvariability",                 "h var.",
    "pad_cv",                 "Complexity",                     "PAD \nvariability",                    "PAD (cv)",
    "ri",                     "Complexity",                     "Rumple index",                         "RI",
    "fhd",                    "Complexity",                     "Foliage height \ndiversity",           "FHD",
    "browsing",               "Other",                          "Browsing",                             "Brows.",
    "lrb_mean",               "Light quantity",                 "Radiative budget \nmean (logged)",     "RB (log)",
    "rb_mean",                "Light quantity",                 "Radiative\nbudget mean",               "RB",
    "rb_cv",                  "Light horizontal heterogeneity", "Radiative\nbudget CV",                 "RB (cv)",
    "rb_sd",                  "Light horizontal heterogeneity", "Radiative\nbudget SD",                 "",
    "ext_ref_agg",            "Light quantity",                 "Mean extinction \n at 4.5m",           "Ext. mean",
    "ext_ref_cv",             "Light horizontal heterogeneity", "Mean extinction \n at 4.5m CV",        "",
    "ext_ref_sd",             "Light horizontal heterogeneity", "Mean extinction \n at 4.5m SD",        "",
    "padAbove_agg",           "PAD",                            "PAD above \nregeneration",             "",
    "padAbove_cv",            "PAD",                            "PAD above \nCV",                       "",
    "padAbove_sd",            "PAD",                            "PAD above \nSD",                       "",
    "padTot_agg",             "PAD",                            "PAD total",                            "",
    "padTot_cv",              "PAD",                            "PAD total CV",                         "",
    "padTot_sd",              "PAD",                            "PAD total sd",                         "",
    "rb_ref_agg",             "Light quantity",                 "RB. at 4.5m",                          "",
    "rb_ref_cv",              "Light horizontal heterogeneity", "RB. at 4.5m \nCV",                     "",
    "rb_ref_sd",              "Light horizontal heterogeneity", "RB. at 4.5m \nSD",                     "",
    "slope_cross_agg",        "Light vertical heterogeneity",   "Light profile slope \nat ext = ",      "",
    "slope_cross_cv",         "Light vertical heterogeneity",   "Light profile slope, CV \nat ext = ",  "",
    "slope_cross_sd",         "Light vertical heterogeneity",   "Light profile slope, SD \nat ext = ",  "",
    "slope_ref_agg",          "Light vertical heterogeneity",   "Light profile slope \nat 4.5m",        "Light profile \nslope",
    "slope_ref_cv",           "Light vertical heterogeneity",   "Light profile slope, CV\nat 4.5m",     "",
    "slope_ref_sd",           "Light vertical heterogeneity",   "Light profile slope, SD\nat 4.5m",     "",
    "z_cross_agg",            "Light",                          "Height\nat ext = ",                    "",
    "z_cross_cv",             "Light",                          "Height, CV\nat ext = ",                "",
    "z_cross_sd",             "Light",                          "Height, SD\nat ext = ",                "",
    "zPAD50_agg",             "PAD",                            "Height at\n 50% of max(PAD)",          "",
    "zPAD50_cv",              "PAD",                            "Height (CV) at\n 50% of max(PAD)",     "",
    "zPAD50_sd",              "PAD",                            "Height (SD) at\n 50% of max(PAD)",     "",
    "gini_h_rb",              "Light horizontal heterogeneity", "Gini index",                           "Gini",
    "gini_v_rb",              "Light vertical heterogeneity",   "Gini vert.",                           "",
    "bl_diff_abs",            "Light vertical heterogeneity",   "Big leaf dif.",                        "Big leaf dif.",
    "bl_slope_abs",           "Light vertical heterogeneity",   "Big leaf \nslope dif.",                "",
    "PC1",                    "Light",                          "Light PCA axis 1",                     "",
    "n_trees_alive",          "Adult",                          "Number of\n living trees",             "",
    "species_richness_alive", "Adult",                          "Adult richness",                       "",
    "mean_dbh_alive_cm",      "Adult",                          "Mean diameter",                        "",
    "mean_height_alive_m",    "Adult",                          "Mean height",                          "",
    "mortality_rate",         "Adult",                          "Mortality",                            "",
    "basal_area_alive_m2",    "Adult",                          "Total basal area",                     ""
  )
}



# --- Dataset builder (unchanged logic) ----------------------------------------

#' Create a analysis table with the targetted variables 
#' @param df first dataframe with lidar/light vars
#' @param df_2 second optional df with forest metrics
#' @param id_cols names of the variables
#' @param plot_level plot or subplot grouping
#' @param t time of the simul
#' @param s section of the simul
#' @param tg targetted extinction levels for some vars

prep_data <- function(df, df_2 = NULL, id_cols, plot_level, t = 12, s = 0, tg = 0.9) {
  df_f <- df %>%
    filter(if (plot_level == "plot") subplot == "plot" else subplot != "plot",
           time == t, slice_num == s,
           if ("target" %in% colnames(df)) target == .env$tg else TRUE) %>%
    dplyr::select(any_of(c("plot", "subplot", id_cols))) %>%
    unique()
  if (!is.null(df_2)) {
    df_2f <- df_2 %>%
      filter(
        if (plot_level == "plot") subplot == "plot" else subplot != "plot",
        if ("time" %in% colnames(df_2)) time == .env$t else TRUE,
        if ("slice_num" %in% colnames(df_2)) slice_num == .env$s else TRUE
      ) %>%
      dplyr::select(any_of(c("plot", "subplot", id_cols))) %>%
      unique()
    df_f <- left_join(df_f, df_2f, by = c("plot", "subplot"))
  }
  df_f
}

# --- Small helpers -------------------------------------------------------------
#' Edit the names of the variables that have several targette values
fix_cross_labels <- function(label, variable) {
  if_else(variable %in% c("slope_cross_agg", "z_cross_agg"), paste0(label, " 0.9"), label)
}

#' Write significance levels
p_stars <- function(p, ns = "ns") {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.1 ~ "(*)", TRUE ~ ns)
}

pretty_name <- function(v, labels) unname(ifelse(v %in% names(labels), labels[v], v))


p_class_colours <- c("< 0.001" = "#08306b", "< 0.01" = "#2171b5",
                     "< 0.05"  = "#6baed6", "< 0.1"  = "#c6dbef", "ns" = "grey78")

pca_axis_lab <- function(pca, k) {
  paste0("PCA axis ", k, " (",
         round(summary(pca)$importance["Proportion of Variance", k] * 100, 1),
         "% of variation)")
}

# common square limits / breaks for a biplot
pca_axes <- function(scores, var_coord, scale_f = 5) {
  v <- c(scores$PC1, scores$PC2, var_coord$PC1 * scale_f, var_coord$PC2 * scale_f)
  r <- range(v, na.rm = TRUE)
  list(breaks = pretty(r, n = 5), limits = c(floor(r[1]), ceiling(r[2])))
}

# labels at arrow tips + invisible points along arrows (so ggrepel avoids them)
arrow_label_df <- function(var_coord, scale_f = 5) {
  tips <- var_coord %>%
    transmute(x = PC1 * scale_f * 1.01, y = PC2 * scale_f * 1.01, abrv, type)
  along <- var_coord %>%
    tidyr::crossing(t = seq(0.05, 1, length.out = 7)) %>%
    transmute(x = PC1 * scale_f * t, y = PC2 * scale_f * t, abrv = "", type)
  bind_rows(tips, along)
}

# --- Themes --------------------------------------------------------------------
theme_pca_paper <- function() {
  theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_line(linewidth = 0.2),
      legend.position = "bottom",
      legend.location = "panel",
      legend.direction = "vertical",
      legend.justification = "left",
      legend.key.spacing = unit(0, "pt"),
      legend.text = element_text(margin = margin(l = 0), size = 8),
      aspect.ratio = 1,
      text = element_text(size = 8)
    )
}

theme_varpart <- function() {
  theme_minimal(base_size = 9) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_line(linewidth = 0.2),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey35", size = 10),
      plot.caption     = element_text(colour = "grey45", size = 9, hjust = 0),
      axis.title       = element_text(colour = "grey25"),
      legend.position  = "bottom"
    )
}

# --- Saving (returns file paths -> usable with format = "file") ----------------
save_fig <- function(plot, name, width, height, dir = "figures",
                     units = "cm", dpi = 600, scale = 1) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  png <- file.path(dir, paste0(name, ".png"))
  pdf <- file.path(dir, paste0(name, ".pdf"))
  ggsave(png, plot, width = width, height = height, units = units,
         dpi = dpi, bg = "white", scale = scale)
  ggsave(pdf, plot, width = width, height = height, units = units,
         device = cairo_pdf, scale = scale)
  c(png, pdf)
}
