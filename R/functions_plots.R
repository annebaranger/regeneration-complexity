#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_figures.R  
#' @description R script containing helpers functions
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# --- Fig 1: regeneration PCA with side boxplots --------------------------------
plot_pca_regen <- function(res, scale_f = 5, kw_x = 2.7, kw_y = 2.5, mid_level = "Medium") {
  pts     <- filter(res$scores, subplot == "plot")
  pts_cat <- filter(pts, !is.na(reg_category))
  vc      <- res$var_coord
  ax      <- pca_axes(pts, vc, scale_f)
  kw1 <- filter(res$kw, axis == "PC1") %>% mutate(x = kw_x, y = mid_level)
  kw2 <- filter(res$kw, axis == "PC2") %>% mutate(x = mid_level, y = kw_y)
  
  ggplot() +
    geom_vline(xintercept = 0, colour = "darkgrey", linewidth = 0.2) +
    geom_hline(yintercept = 0, colour = "darkgrey", linewidth = 0.2) +
    geom_point(data = pts, aes(PC1, PC2, color = reg_category), size = 0.5, shape = 19) +
    ggside::geom_xsideboxplot(
      data = pts_cat, aes(x = PC1, y = reg_category, fill = reg_category),
      whisker.linewidth = 0.2, median.linewidth = 0.2, box.linewidth = 0.2,
      orientation = "y", outliers = FALSE) +
    ggside::geom_ysideboxplot(
      data = pts_cat, aes(y = PC2, x = reg_category, fill = reg_category),
      whisker.linewidth = 0.2, median.linewidth = 0.2, box.linewidth = 0.2,
      orientation = "x", outliers = FALSE) +
    ggside::geom_xsidetext(data = kw1, aes(x = x, y = y, label = label),
                           size = 7.5, size.unit = "pt", vjust = -0.02, inherit.aes = FALSE) +
    ggside::geom_ysidetext(data = kw2, aes(x = x, y = y, label = label),
                           size = 7.5, size.unit = "pt", inherit.aes = FALSE) +
    scale_fill_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
    scale_color_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
    guides(color = "none") +
    ggnewscale::new_scale_color() +
    geom_segment(data = vc,
                 aes(x = 0, y = 0, xend = PC1 * scale_f, yend = PC2 * scale_f,
                     color = type, linewidth = Source),
                 arrow = arrow(length = unit(0.10, "cm")), lineend = "round") +
    ggrepel::geom_text_repel(data = arrow_label_df(vc, scale_f),
                             aes(x = x, y = y, label = abrv, color = type),
                             size=9 * 0.353, segment.color = "grey50", segment.size = 0.1,
                             seed = 58, show.legend = FALSE) +
    scale_color_brewer(name = "Variable type", palette = "Dark2") +
    scale_linewidth_manual(name = "Data source",
                           values = c("Lidar-derived" = 0.45, "Inventory-derived" = 0.2)) +
    theme_pca_paper() +
    theme(ggside.panel.scale = 0.15,
          ggside.panel.grid = element_blank(),
          ggside.panel.border = element_blank(),
          ggside.panel.background = element_blank(),
          ggside.axis.text = element_blank(),
          ggside.axis.ticks = element_blank()) +
    labs(x = pca_axis_lab(res$pca, 1), y = pca_axis_lab(res$pca, 2)) +
    scale_x_continuous(limits = ax$limits, breaks = ax$breaks) +
    scale_y_continuous(limits = ax$limits, breaks = ax$breaks)
}



# --- Light PCA biplot (plot or subplot) ----------------------------------------
plot_pca_light <- function(res, scale_f = 5, point_size = 0.8) {
  ax <- pca_axes(res$scores, res$var_coord, scale_f)
  ggplot() +
    geom_point(data = res$scores, aes(PC1, PC2, color = reg_category), size = point_size) +
    scale_color_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
    ggnewscale::new_scale_color() +
    geom_segment(data = res$var_coord,
                 aes(x = 0, y = 0, xend = PC1 * scale_f, yend = PC2 * scale_f, color = type),
                 arrow = arrow(length = unit(0.10, "cm")), lineend = "round", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = arrow_label_df(res$var_coord, scale_f),
                             aes(x = x, y = y, label = abrv, color = type),
                             size=9 * 0.353, show.legend = FALSE) +
    scale_color_brewer(name = "Variable type", palette = "Dark2") +
    # guides(color = "none") +
    theme_pca_paper() +
    labs(x = pca_axis_lab(res$pca, 1), y = pca_axis_lab(res$pca, 2)) +
    scale_x_continuous(limits = ax$limits, breaks = ax$breaks) +
    scale_y_continuous(limits = ax$limits, breaks = ax$breaks)
}

# --- Light descriptors by regeneration class -----------------------------------
plot_light_boxes <- function(box_df, kw = NULL, cld = NULL) {
  p <- ggplot(box_df, aes(label, value, fill = reg_category)) +
    geom_boxplot(position = position_dodge(width = 0.8, preserve = "single"),
                 whisker.linewidth = 0.2, median.linewidth = 0.2, box.linewidth = 0.2,
                 width = 0.5, outlier.shape = NA) +
    geom_point(aes(color = reg_category),
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1),
               size = 1, alpha = 0.8, shape = 16, show.legend = FALSE)
  if (!is.null(kw))
    p <- p + geom_text(data = kw, aes(x = label, y = 2.5, label = stars),
                       inherit.aes = FALSE, size = 3, vjust = 0)
  if (!is.null(cld))
    p <- p + geom_text(data = cld,
                       aes(x = label, y = y + 0.25, label = letter, group = reg_category),
                       position = position_dodge(width = 0.8),
                       inherit.aes = FALSE, size = 3, vjust = 0)
  p +
    scale_fill_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
    scale_color_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
    facet_grid(~type, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(y = "Scaled values", x = "") +
    theme_pca_paper() +
    theme(legend.position = "none",
          base_size = 9)
}

plot_light_panel <- function(pca_plot, box_plot) {
  cowplot::plot_grid(pca_plot + theme(legend.position = "left"), box_plot,
                     ncol = 2, align = "h", axis = "tb", rel_widths = c(1, 1.5))
}

# --- Bootstrapped effect plot (light or adults, plot or subplot) ----------------
plot_effects <- function(stats, var_meta, responses, predictors, drop = NULL,
                         type_levels = c("Composition", "Quantity", "Complexity"),
                         lw = c(0.2, 0.5)) {
  predictors <- setdiff(predictors, drop)
  lab <- dplyr::select(var_meta, variable, label, type)
  ord <- unique(c(lab$label[match(predictors, lab$variable)],
                  lab$label[match(responses,  lab$variable)]))
  
  stats %>%
    filter(explaining %in% predictors) %>%
    mutate(p_class = cut(p_boot, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         labels = names(p_class_colours),
                         right = FALSE, include.lowest = TRUE)) %>%
    left_join(rename(lab, explained = variable, label_explained = label,
                     type_explained = type), by = "explained") %>%
    left_join(dplyr::select(lab, explaining = variable, label_explaining = label),
              by = "explaining") %>%
    mutate(type_explained   = factor(type_explained, levels = type_levels),
           label_explained  = factor(label_explained, levels = rev(ord)),
           label_explaining = factor(label_explaining, levels = ord)) %>%
    ggplot(aes(x = med, y = label_explained)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = .2, colour = "grey70") +
    geom_linerange(aes(xmin = lo,  xmax = hi,  colour = p_class), linewidth = lw[1]) +
    geom_linerange(aes(xmin = lo5, xmax = hi5, colour = p_class), linewidth = lw[2]) +
    geom_point(aes(colour = p_class), size = 1.3) +
    scale_colour_manual(values = p_class_colours, drop = FALSE, name = "p (bootstrap)") +
    facet_grid(type_explained ~ label_explaining, scales = "free_y", space = "free_y",switch = "y") +
    labs(x = "Standardized coefficient", y = NULL) +
    theme_minimal() +
    theme(base_size = 9,
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          panel.grid = element_line(linewidth = 0.2),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0),
          strip.text=element_text(face="bold"),
          strip.placement   = "outside",                  # strips go outside the y-axis labels
          strip.text.y.left = element_text(angle = 90),
          legend.position = "bottom",
          text = element_text(size=9))
}

# --- Temporal vs spatial variability -------------------------------------------
plot_variability <- function(res, var_meta,
                             levels = c("Mean extinction \n at 4.5m", "Gini index",
                                        "Light profile slope \nat 4.5m")) {
  res %>%
    left_join(var_meta, by = c("metric" = "variable")) %>%
    mutate(label = factor(label, levels = levels)) %>%
    ggplot(aes(x = spatial_within, y = temporal_within)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 2.5) +
    facet_wrap(~ label) +
    expand_limits(x = 0, y = 0) +
    labs(x = "Spatial variability (between subplots)",
         y = "Temporal variability (between simulation times)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_blank())
}

# --- Variance partitioning figures ---------------------------------------------
vp_group_colours    <- c("Light" = "#E0A526", "Adults" = "#2F6B4F")
vp_fraction_colours <- c("Light (unique)" = "#E0A526", "Shared" = "#A7A77A",
                         "Adults (unique)" = "#2F6B4F")

plot_varpart_venn <- function(vp) {
  multi <- vp$multi
  fmt   <- function(v) ifelse(v < 0, sprintf("%.2f (≈ 0)", v), sprintf("%.2f", v))
  circle <- function(x0, g) {
    t <- seq(0, 2 * pi, length.out = 200)
    data.frame(x = x0 + cos(t), y = sin(t), group = g)
  }
  ci_txt <- function(t) with(vp$multi_jack[vp$multi_jack$term == t, ],
                             sprintf("[%.2f; %.2f]", ci_low, ci_high))
  labs_df <- data.frame(
    x = c(-1.05, 0, 1.05), y = 0,
    text = paste0(
      c(sprintf("[a]\n%s %s", fmt(multi$light_unique), p_stars(multi$p_light_unique, "n.s.")),
        sprintf("[b]\n%s",    fmt(multi$shared)),
        sprintf("[c]\n%s %s", fmt(multi$adult_unique), p_stars(multi$p_adult_unique, "n.s."))),
      "\n", c(ci_txt("light_unique"), ci_txt("shared"), ci_txt("adult_unique"))))
  
  ggplot() +
    geom_polygon(data = rbind(circle(-0.55, "Light"), circle(0.55, "Adults")),
                 aes(x, y, fill = group, group = group),
                 alpha = 0.35, colour = "white", linewidth = 1.2) +
    geom_text(data = labs_df, aes(x, y, label = text), size = 9 * 0.353, fontface = "bold",
              colour = "grey15", lineheight = 0.9) +
    geom_text(data = data.frame(x = c(-0.95, 0.95), y = 1.2, group = c("Light", "Adults")),
              aes(x, y, label = group, colour = group), size = 10 * 0.353, fontface = "bold") +
    annotate("text", x = 0, y = -1.3, colour = "grey40", size=9 * 0.353,
             label = sprintf("Residuals = %.2f \nAdjusted R² [jackknife 95 %% CI] \n full model: adj. R² = %.2f, p = %.3f",
                             multi$residuals,multi$full_adj_r2, multi$p_full)) +
    scale_fill_manual(values = vp_group_colours) +
    scale_colour_manual(values = vp_group_colours) +
    coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.45, 1.4)) +
    theme_void() +
    theme(legend.position = "none",
          text=element_text(size=9))
}

plot_varpart_univariate <- function(vp, labels, alpha = 0.1) {
  uf <- vp$uni_fractions
  lev <- vp$regen_vars
  pn  <- function(v) pretty_name(v, labels)
  
  d <- bind_rows(
    data.frame(response = uf$response, fraction = "Light (unique)",
               value = uf$light_unique, p = uf$p_light_unique_holm),
    data.frame(response = uf$response, fraction = "Shared", value = uf$shared, p = NA),
    data.frame(response = uf$response, fraction = "Adults (unique)",
               value = uf$adult_unique, p = uf$p_adult_unique_holm)) %>%
    mutate(fraction = factor(fraction, levels = names(vp_fraction_colours)),
           response = factor(response, levels = lev),
           y_star   = pmax(value, 0, na.rm = TRUE) + 0.04)
  top <- uf %>%
    mutate(response = factor(response, levels = lev),
           label = sprintf("adj. R² = %.2f %s", full_adj_r2, p_stars(p_full, "n.s.")))
  y_top <- max(d$y_star) + 0.12
  
  fig_a <- ggplot(d, aes(x = response, y = value, fill = fraction)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(y = y_star, label = p_stars(p, "n.s."), group = fraction),
              position = position_dodge(width = 0.8), vjust = 0, size=9 * 0.353, colour = "grey20") +
    geom_text(data = top, aes(x = response, y = y_top, label = label),
              inherit.aes = FALSE, size=9 * 0.353, colour = "grey25") +
    scale_fill_manual(values = vp_fraction_colours, name = NULL) +
    scale_x_discrete(labels = pn) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = "Variance explained (adjusted R²)")+
         #title = "Variance explained by light and adults",
         #subtitle = "Stars: permutation tests of unique fractions,Holm-corrected across responses") +
    theme_varpart() +
    theme(panel.grid.major.x = element_blank())
  
  fig_b <-vp$uni_importance %>%
    mutate(response = factor(response, levels = lev),
           variable = factor(variable, levels = c(rev(vp$light_vars), rev(vp$adult_vars))),
           group    = factor(group, levels = c("Adults", "Light")),
           label    = sprintf("#%d%s", rank,
                              if_else(p_perm > alpha, "", paste0("\n", p_stars(p_perm, "n.s."))))) %>%
    ggplot(aes(x = response, y = variable, fill = individual)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = label), size=9 * 0.353, lineheight = 0.85, colour = "grey10") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y",switch = "y") +
    scale_fill_gradient2(low = "#B2443A", mid = "white", high = "#2B4C7E", midpoint = 0,
                         name = "Individual \nimportance\n(adjusted R²)") +
    scale_x_discrete(labels = pn, position = "top") +
    scale_y_discrete(labels = pn) +
    labs(x = NULL, y = NULL) + #, title = "Importance ranking of predictors",
        # subtitle = "#: rank within each response  |  Stars: permutation test") +
    theme_varpart() +
    theme(panel.grid = element_blank(),
          strip.text.y = element_text(face = "bold", angle = 0, size = 11),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.25, "cm"),
          strip.placement   = "outside",                  # strips go outside the y-axis labels
          strip.text.y.left = element_text(angle = 90),
          text=element_text(zie=8))
  
  patchwork::wrap_plots(fig_a, fig_b, widths = c(1.3, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
}

plot_rda_triplot <- function(vp, labels, regen_colour = "#2B4C7E") {
  model <- vp$model
  pn    <- function(v) pretty_name(v, labels)
  pct   <- round(100 * model$CCA$eig[1:2] / model$tot.chi, 1)
  pax   <- vp$axis_test$`Pr(>F)`[1:2]
  sc    <- vegan::scores(model, scaling = 2, choices = 1:2,
                         display = c("sites", "species", "bp"))
  ren   <- function(d) { d <- as.data.frame(d); colnames(d)[1:2] <- c("axis1", "axis2"); d$name <- rownames(d); d }
  sites <- ren(sc$sites); resp <- ren(sc$species); bp <- ren(sc$biplot)
  bp$group <- ifelse(bp$name %in% vp$adult_vars, "Adults", "Light")
  spread <- 0.85 * max(abs(sites[, 1:2]))
  bp[, 1:2]   <- bp[, 1:2]   * spread / max(abs(bp[, 1:2]))
  resp[, 1:2] <- resp[, 1:2] * spread / max(abs(resp[, 1:2]))
  
  regen_lab <- "Regeneration"
  col_vals  <- c(vp_group_colours, setNames(regen_colour, regen_lab))
  
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(data = sites, aes(axis1, axis2), colour = "grey60", size = 2, alpha = 0.8) +
    geom_segment(data = bp, aes(x = 0, y = 0, xend = axis1, yend = axis2, colour = group),
                 arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.8) +
    geom_segment(data = resp, aes(x = 0, y = 0, xend = axis1, yend = axis2, colour = regen_lab),
                 arrow = arrow(length = unit(0.22, "cm")), linewidth = 1.2) +
    ggrepel::geom_text_repel(data = bp, aes(axis1, axis2, label = pn(name), colour = group),
                             size = 3.4, show.legend = FALSE, seed = 1) +
    ggrepel::geom_text_repel(data = resp, aes(axis1, axis2, label = pn(name), colour = regen_lab),
                             fontface = "bold", size = 3.8, show.legend = FALSE, seed = 1) +
    scale_colour_manual(values = col_vals, breaks = names(col_vals), name = NULL) +
    coord_equal() +
    labs(x = sprintf("RDA axis 1 (%.1f %%, p = %.3f)", pct[1], pax[1]),
         y = sprintf("RDA axis 2 (%.1f %%, p = %.3f)", pct[2], pax[2])) +
    theme_varpart()
}