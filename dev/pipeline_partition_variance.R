# =============================================================================
# Pipeline: effects of light and adult stand structure
#           on forest regeneration structure (RDA + variance partitioning)
# =============================================================================
# Part A - Multivariate analysis: regen [richness, abundance, ri] ~ light + adults
#   A1. Global test of the full model
#   A2. Variance partitioning light / adults with permutation tests
#   Figure 1: partitioning diagram
#
# Part B - Univariate analyses: the same for each regeneration variable
#   B1. Global test + variance partitioning per response (Holm across responses)
#   B2. Hierarchical partitioning: importance ranking of predictors per response
#   Figure 2: (A) variance explained by light / adults, (B) importance ranking
#
# Part C - Jackknife (leave-one-out): uncertainty and influential plots
#   Applied to multivariate fractions, univariate fractions and importances
#
# Additional figures
#   Figure 3: multivariate RDA triplot
#   Figure 4: standardised coefficients of the univariate models
#   Figure 5: predictor importance with jackknife 95 % CI
# =============================================================================

library(vegan)
library(rdacca.hp)
library(dplyr)
library(ggplot2)
library(patchwork)   # multi-panel figures
library(ggrepel)     # non-overlapping labels

set.seed(42)         # reproducible permutations
n_perm    <- 999    # permutations for RDA tests
n_perm_hp <- 999     # permutations for permu.hp (slower)
alpha     <- 0.1

# --- Parameters: variable names ----------------------------------------------
adult_vars <- c("species_richness_alive", "mean_dbh_alive_cm", "mean_height_alive_m")
light_vars <- c("lrb_mean", "gini_h_rb", "slope_ref_agg")
regen_vars <- c("richness", "abundance", "ri")
log_vars   <- c("abundance")   # responses to log1p-transform (character(0) if none)

# Labels displayed in figures (adapt, especially gini / slope)
display_names <- c(
  species_richness_alive = "Adult species richness",
  mean_dbh_alive_cm      = "Mean DBH",
  mean_height_alive_m    = "Mean height",
  lrb_mean                = "Mean light",
  gini_h_rb              = "Horizontal heterogeneity",
  slope_ref_agg          = "Vertical heterogeneity",
  richness               = "Richness",
  abundance              = "Abundance",
  ri                     = "Complexity"
)

# =============================================================================
# 0. Data preparation
# =============================================================================
dat <- prep_data(tls_metrics,forest_data_plot,id_cols = c(regen_vars, adult_vars, light_vars),
                 plot_level = "plot") 

# Plots excluded because of missing values (to report in the paper)
na_plots <- which(!complete.cases(dat))
if (length(na_plots) > 0)
  message("Plots excluded (NA): ", paste(na_plots, collapse = ", "))
dat <- dat[complete.cases(dat), ]
n <- nrow(dat)
message("Number of plots analysed: n = ", n)

# Transform skewed responses, then standardise (mean 0, variance 1)
standardise <- function(d) d %>% mutate(across(everything(), ~ as.numeric(scale(.))))
regen  <- dat %>% select(all_of(regen_vars)) %>%
  mutate(across(all_of(log_vars), log1p)) %>% standardise()
adults <- dat %>% select(all_of(adult_vars)) %>% standardise()
light  <- dat %>% select(all_of(light_vars)) %>% standardise()
predictors <- cbind(adults, light)

# Collinearity among predictors (rule of thumb: VIF < 5 comfortable, > 10 problematic)
print(round(diag(solve(cor(predictors))), 2))

# =============================================================================
# Shared helpers
# =============================================================================
# Global test + variance partitioning for any response (matrix or single vector)
partition_analysis <- function(Y) {
  full <- rda(Y ~ ., data = predictors)
  vp   <- varpart(Y, light, adults)
  f    <- vp$part$fract$Adj.R.square   # [1] light alone, [2] adults alone, [3] both
  test <- function(m) anova(m, permutations = n_perm)$`Pr(>F)`[1]
  data.frame(
    full_adj_r2    = f[3],
    p_full         = test(full),
    light_total    = f[1],
    p_light_total  = test(rda(Y ~ ., data = light)),
    adult_total    = f[2],
    p_adult_total  = test(rda(Y ~ ., data = adults)),
    light_unique   = f[3] - f[2],
    p_light_unique = test(rda(Y ~ . + Condition(as.matrix(adults)), data = light)),
    shared         = f[1] + f[2] - f[3],             # not testable
    adult_unique   = f[3] - f[1],
    p_adult_unique = test(rda(Y ~ . + Condition(as.matrix(light)), data = adults)),
    residuals      = 1 - f[3]
  )
}

stars <- function(p, ns = "n.s.") ifelse(is.na(p), "",
                                         ifelse(p < 0.001, "***", 
                                                ifelse(p < 0.01, "**", 
                                                       ifelse(p < 0.05, "*",
                                                              ifelse(p < 0.1, "(*)", ns)))))
pretty_name <- function(v) unname(ifelse(v %in% names(display_names), display_names[v], v))

# =============================================================================
# PART A - Multivariate analysis
# =============================================================================
stopifnot(ncol(regen) == length(regen_vars))
model <- rda(regen ~ ., data = predictors)
stopifnot(abs(model$tot.chi - ncol(regen)) < 1e-6)   # total variance = number of responses

# A1. Global test
global_test <- anova(model, permutations = n_perm)
print(global_test)
model_r2 <- RsquareAdj(model)
message(sprintf("Multivariate model: raw R² = %.3f | adjusted R² = %.3f | p = %.4f",
                model_r2$r.squared, model_r2$adj.r.squared, global_test$`Pr(>F)`[1]))
if (global_test$`Pr(>F)`[1] > alpha)
  warning("Multivariate model NOT significant: partitioning results are exploratory only.")

# A2. Variance partitioning with tests
multi <- partition_analysis(regen)
multi_fractions <- data.frame(
  fraction = c("Full model [a+b+c]", "Light (total) [a+b]", "Adults (total) [b+c]",
               "Light (unique) [a]", "Shared [b]", "Adults (unique) [c]", "Residuals [d]"),
  adj_r2   = round(c(multi$full_adj_r2, multi$light_total, multi$adult_total,
                     multi$light_unique, multi$shared, multi$adult_unique, multi$residuals), 3),
  p_value  = c(multi$p_full, multi$p_light_total, multi$p_adult_total,
               multi$p_light_unique, NA, multi$p_adult_unique, NA)
)
print(multi_fractions)
# Reminder: [b] cannot be tested; a negative adjusted fraction is interpreted as 0.

# =============================================================================
# PART B - Univariate analyses (one per regeneration variable)
# =============================================================================
uni_fractions  <- list()
uni_importance <- list()

for (v in regen_vars) {
  y <- regen[[v]]
  
  # B1. Global test + partitioning
  uni_fractions[[v]] <- cbind(response = v, partition_analysis(y))
  
  # B2. Hierarchical partitioning + permutation test of each predictor
  hp   <- rdacca.hp(y, predictors, method = "RDA", type = "adjR2")$Hier.part
  perm <- permu.hp(y, predictors, method = "RDA", type = "adjR2", permutations = n_perm_hp)
  p_col <- grep("^Pr", colnames(perm), value = TRUE)[1]
  uni_importance[[v]] <- data.frame(
    response   = v,
    variable   = rownames(hp),
    group      = ifelse(rownames(hp) %in% adult_vars, "Adults", "Light"),
    unique     = hp[, "Unique"],
    shared     = hp[, "Average.share"],
    individual = hp[, "Individual"],
    percent    = hp[, "I.perc(%)"],
    p_perm     = as.numeric(gsub("[^0-9.eE-]", "", perm[rownames(hp), p_col]))
  )
}

# Holm correction across the 3 responses (one family per type of test)
uni_fractions <- bind_rows(uni_fractions) %>%
  mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "holm"), .names = "{.col}_holm"))
print(uni_fractions, digits = 3)

# Importance: Holm correction across the 6 predictors within each response, then ranking
uni_importance <- bind_rows(uni_importance) %>%
  group_by(response) %>%
  mutate(p_holm = p.adjust(p_perm, method = "holm"),
         rank   = rank(-individual, ties.method = "min")) %>%
  ungroup() %>%
  arrange(response, rank)
print(as.data.frame(uni_importance), digits = 3)
# Only interpret the ranking for responses whose global model is significant.

# =============================================================================
# PART C - Jackknife (leave-one-out)
# =============================================================================
# Each plot is dropped in turn and the analysis refitted on the n - 1 others.
# - jackknife standard error -> approximate 95 % CI (estimate +/- 1.96 SE)
# - min / max over the n refits and the most influential plot -> sensitivity
# (The classic bootstrap is unsuitable with n = 19: resamples contain ~12 distinct
#  plots, which biases the adjusted R² and sometimes saturates the model.)

drop_k <- function(Y, k) if (is.null(dim(Y))) Y[-k] else Y[-k, , drop = FALSE]

# Fractions refitted without plot k (matrix: fractions x plots)
jk_fractions <- function(Y) {
  sapply(seq_len(n), function(k) {
    Yk <- drop_k(Y, k)
    r2 <- function(X) RsquareAdj(rda(Yk ~ ., data = X[-k, , drop = FALSE]))$adj.r.squared
    f_light <- r2(light); f_adult <- r2(adults); f_full <- r2(predictors)
    c(light_unique = f_full - f_adult, shared = f_light + f_adult - f_full,
      adult_unique = f_full - f_light, full_adj_r2 = f_full)
  })
}

# Summary of a jackknife matrix given the full-data estimates (named vector)
jk_summary <- function(jk, estimate) {
  est <- estimate[rownames(jk)]
  se  <- sqrt((n - 1) / n * rowSums((jk - rowMeans(jk))^2))
  data.frame(term = rownames(jk), estimate = est, se_jack = se,
             ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
             min_loo = apply(jk, 1, min), max_loo = apply(jk, 1, max),
             most_influential_plot = rownames(dat)[apply(abs(jk - est), 1, which.max)],
             row.names = NULL)
}

# C1. Multivariate fractions
multi_jack <- jk_summary(jk_fractions(regen),
                         c(light_unique = multi$light_unique, shared = multi$shared,
                           adult_unique = multi$adult_unique, full_adj_r2 = multi$full_adj_r2))
print(multi_jack, digits = 3)

# C2. Univariate fractions
uni_frac_jack <- bind_rows(lapply(regen_vars, function(v) {
  est <- unlist(uni_fractions[uni_fractions$response == v,
                              c("light_unique", "shared", "adult_unique", "full_adj_r2")])
  cbind(response = v, jk_summary(jk_fractions(regen[[v]]), est))
}))
print(uni_frac_jack, digits = 3)

# C3. Univariate importances
uni_imp_jack <- bind_rows(lapply(regen_vars, function(v) {
  y  <- regen[[v]]
  jk <- sapply(seq_len(n), function(k)
    rdacca.hp(y[-k], predictors[-k, ], method = "RDA", type = "adjR2")$Hier.part[, "Individual"])
  imp <- uni_importance[uni_importance$response == v, ]
  cbind(response = v, jk_summary(jk, setNames(imp$individual, imp$variable)))
}))
uni_importance <- uni_importance %>%
  left_join(uni_imp_jack %>% select(response, variable = term, se_jack, ci_low, ci_high,
                                    min_loo, max_loo, most_influential_plot),
            by = c("response", "variable"))
print(as.data.frame(uni_importance), digits = 3)

# =============================================================================
# FIGURES
# =============================================================================
dir.create("figures", showWarnings = FALSE)

group_colours    <- c("Light" = "#E0A526", "Adults" = "#2F6B4F")
fraction_colours <- c("Light (unique)" = "#E0A526", "Shared" = "#A7A77A",
                      "Adults (unique)" = "#2F6B4F")

theme_set(theme_minimal(base_size = 12) + theme(
  panel.grid.minor = element_blank(),
  plot.title       = element_text(face = "bold", size = 13),
  plot.subtitle    = element_text(colour = "grey35", size = 10),
  plot.caption     = element_text(colour = "grey45", size = 9, hjust = 0),
  axis.title       = element_text(colour = "grey25"),
  legend.position  = "bottom"))

save_fig <- function(plot, name, width, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0(name, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
}

# --- Figure 1: multivariate partitioning diagram ------------------------------
# (schematic circles: their size is not proportional to the fractions)
fmt_frac <- function(v) ifelse(v < 0, sprintf("%.2f (≈ 0)", v), sprintf("%.2f", v))
circle <- function(x0, g) {
  t <- seq(0, 2 * pi, length.out = 200)
  data.frame(x = x0 + cos(t), y = sin(t), group = g)
}
circles <- rbind(circle(-0.55, "Light"), circle(0.55, "Adults"))
venn_labels <- data.frame(
  x = c(-1.05, 0, 1.05), y = 0,
  text = c(sprintf("[a]\n%s %s", fmt_frac(multi$light_unique), stars(multi$p_light_unique)),
           sprintf("[b]\n%s",    fmt_frac(multi$shared)),
           sprintf("[c]\n%s %s", fmt_frac(multi$adult_unique), stars(multi$p_adult_unique))))
ci_txt <- function(t) with(multi_jack[multi_jack$term == t, ], sprintf("[%.2f; %.2f]", ci_low, ci_high))
venn_labels$text <- paste0(venn_labels$text, "\n", c(ci_txt("light_unique"), ci_txt("shared"), ci_txt("adult_unique")))
titles <- data.frame(x = c(-0.95, 0.95), y = 1.2, group = c("Light", "Adults"))

fig_partition <- ggplot() +
  geom_polygon(data = circles, aes(x, y, fill = group, group = group),
               alpha = 0.35, colour = "white", linewidth = 1.2) +
  geom_text(data = venn_labels, aes(x, y, label = text),
            size = 4.2, fontface = "bold", colour = "grey15", lineheight = 0.9) +
  geom_text(data = titles, aes(x, y, label = group, colour = group),
            size = 5, fontface = "bold") +
  annotate("text", x = 0, y = -1.3, colour = "grey40", size = 3.8,
           label = sprintf("Residuals = %.2f", multi$residuals)) +
  scale_fill_manual(values = group_colours) +
  scale_colour_manual(values = group_colours) +
  coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.45, 1.4)) +
  labs(title = "Variance partitioning of regeneration structure (multivariate RDA)",
       subtitle = sprintf("Adjusted R² [jackknife 95 %% CI]  |  full model: adj. R² = %.2f, p = %.3f",
                          multi$full_adj_r2, multi$p_full)) +
  theme_void() +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(colour = "grey35", size = 10, hjust = 0.5))
print(fig_partition)
save_fig(fig_partition, "fig1_multivariate_partitioning", 7, 5)

# --- Figure 2A: variance explained per response -------------------------------
resp_levels <- regen_vars
fig2a_data <- bind_rows(
  data.frame(response = uni_fractions$response, fraction = "Light (unique)",
             value = uni_fractions$light_unique, p = uni_fractions$p_light_unique_holm),
  data.frame(response = uni_fractions$response, fraction = "Shared",
             value = uni_fractions$shared, p = NA),
  data.frame(response = uni_fractions$response, fraction = "Adults (unique)",
             value = uni_fractions$adult_unique, p = uni_fractions$p_adult_unique_holm)) %>%
  mutate(term = c("Light (unique)" = "light_unique", "Shared" = "shared",
                  "Adults (unique)" = "adult_unique")[fraction]) %>%
  left_join(uni_frac_jack %>% select(response, term, ci_low, ci_high), by = c("response", "term")) %>%
  mutate(fraction = factor(fraction, levels = names(fraction_colours)),
         response = factor(response, levels = resp_levels),
         negative = value < 0)

# offset <- 0.04 * max(abs(c(fig2a_data$value, fig2a_data$ci_high)), na.rm = TRUE)
# fig2a_data$y_star <- pmax(fig2a_data$ci_high, fig2a_data$value, 0, na.rm = TRUE) + offset
offset <- 0.04 
fig2a_data$y_star <- pmax( fig2a_data$value, 0, na.rm = TRUE) + offset



# Full-model label above each response
fig2a_top <- uni_fractions %>%
  mutate(response = factor(response, levels = resp_levels),
         label = sprintf("Model: adj. R² = %.2f %s", full_adj_r2, stars(p_full)))
y_top <- max(fig2a_data$y_star) + 3 * offset

fig2a <- ggplot(fig2a_data, aes(x = response, y = value, fill = fraction)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  # geom_errorbar(aes(ymin = ci_low, ymax = ci_high, group = fraction),
  #               position = position_dodge(width = 0.8), width = 0.25, colour = "grey25") +
  geom_text(aes(y = y_star, label = stars(p, ns = "n.s."), group = fraction),
            position = position_dodge(width = 0.8), vjust = 0, size = 3.6, colour = "grey20") +
  geom_text(data = fig2a_top, aes(x = response, y = y_top, label = label),
            inherit.aes = FALSE, size = 3.3, colour = "grey25") +
  scale_fill_manual(values = fraction_colours, name = NULL) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4), guide = "none") +
  scale_x_discrete(labels = pretty_name) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Variance explained (adjusted R²)",
       title = "Variance explained by light and adults",
       subtitle = "Stars: permutation tests of unique fractions,\nHolm-corrected across responses") +
  theme(panel.grid.major.x = element_blank())

# --- Figure 2B: importance ranking per response -------------------------------
var_order <- c(rev(light_vars), rev(adult_vars))   # adults on top, then light
fig2b_data <- uni_importance %>%
  mutate(response = factor(response, levels = resp_levels),
         variable = factor(variable, levels = var_order),
         group    = factor(group, levels = c("Adults", "Light")),
         label    = sprintf("#%d%s", rank, if_else(p_perm>alpha,"",paste0("\n", stars(p_perm)))))

fig2b <- ggplot(fig2b_data, aes(x = response, y = variable, fill = individual)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 3.6, lineheight = 0.85, colour = "grey10") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#B2443A", mid = "white", high = "#2B4C7E", midpoint = 0,
                       name = "Individual importance\n(adjusted R²)") +
  scale_x_discrete(labels = pretty_name, position = "top") +
  scale_y_discrete(labels = pretty_name) +
  labs(x = NULL, y = NULL,
       title = "Importance ranking of predictors",
       subtitle = "#: rank within each response  |  Stars: permutation test, Holm-corrected") +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(face = "bold", angle = 0, size = 11),
        legend.key.width = unit(1.2, "cm"))

# --- Figure 2: combined --------------------------------------------------------
fig_univariate <- fig2a + fig2b +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(tag_levels = "A")
print(fig_univariate)
save_fig(fig_univariate, "fig2_univariate_partitioning_importance", 14, 6)

# --- Figure 3: multivariate RDA triplot (scaling 2) ---------------------------
# Reading: angle between arrows ~ correlation (acute = positive, 90° = none,
# obtuse = negative); only interpret significant axes (see axis_test).
axis_test <- anova(model, by = "axis", permutations = n_perm)
print(axis_test)

regen_colour <- "#2B4C7E"
axis_pct   <- round(100 * model$CCA$eig[1:2] / model$tot.chi, 1)
axis_p     <- axis_test$`Pr(>F)`[1:2]
rda_scores <- scores(model, scaling = 2, choices = 1:2, display = c("sites", "species", "bp"))
rename_axes <- function(d) { colnames(d)[1:2] <- c("axis1", "axis2"); d }
plots     <- rename_axes(as.data.frame(rda_scores$sites))
responses <- rename_axes(as.data.frame(rda_scores$species)); responses$name <- rownames(responses)
arrows_bp <- rename_axes(as.data.frame(rda_scores$biplot));  arrows_bp$name <- rownames(arrows_bp)
arrows_bp$group <- ifelse(arrows_bp$name %in% adult_vars, "Adults", "Light")
spread <- 0.85 * max(abs(plots[, 1:2]))   # scale arrows to the spread of plots
arrows_bp[, 1:2] <- arrows_bp[, 1:2] * spread / max(abs(arrows_bp[, 1:2]))
responses[, 1:2] <- responses[, 1:2] * spread / max(abs(responses[, 1:2]))

fig_triplot <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(data = plots, aes(axis1, axis2), colour = "grey60", size = 2, alpha = 0.8) +
  geom_segment(data = arrows_bp, aes(x = 0, y = 0, xend = axis1, yend = axis2, colour = group),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.8) +
  geom_segment(data = responses, aes(x = 0, y = 0, xend = axis1, yend = axis2),
               colour = regen_colour, arrow = arrow(length = unit(0.22, "cm")), linewidth = 1.2) +
  geom_text_repel(data = arrows_bp, aes(axis1, axis2, label = pretty_name(name), colour = group),
                  size = 3.4, show.legend = FALSE, seed = 1) +
  geom_text_repel(data = responses, aes(axis1, axis2, label = pretty_name(name)),
                  colour = regen_colour, fontface = "bold", size = 3.8, seed = 1) +
  scale_colour_manual(values = group_colours, name = NULL) +
  coord_equal() +
  labs(x = sprintf("RDA axis 1 (%.1f %%, p = %.3f)", axis_pct[1], axis_p[1]),
       y = sprintf("RDA axis 2 (%.1f %%, p = %.3f)", axis_pct[2], axis_p[2]),
       title = "Multivariate RDA triplot",
       subtitle = "Blue: regeneration variables  |  Points: plots  |  Scaling 2")
print(fig_triplot)
save_fig(fig_triplot, "fig3_multivariate_triplot", 7, 6.5)

# --- Figure 4: standardised coefficients of the univariate models --------------
# A univariate RDA has a single constrained axis (no 2D triplot); it is a multiple
# regression, so its standardised coefficients give direction and strength of effects.
coef_data <- bind_rows(lapply(regen_vars, function(v) {
  m  <- lm(regen[[v]] ~ ., data = predictors)
  ci <- confint(m)[-1, , drop = FALSE]
  data.frame(response = v, variable = names(coef(m))[-1], estimate = coef(m)[-1],
             ci_low = ci[, 1], ci_high = ci[, 2], row.names = NULL)
})) %>%
  mutate(group    = ifelse(variable %in% adult_vars, "Adults", "Light"),
         variable = factor(variable, levels = var_order),
         response = factor(response, levels = resp_levels, labels = pretty_name(resp_levels)))

fig_coef <- ggplot(coef_data, aes(x = estimate, y = variable, colour = group)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.25, orientation = "y", linewidth = 0.7) +
  geom_point(size = 2.8) +
  facet_wrap(~ response, nrow = 1) +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_y_discrete(labels = pretty_name) +
  labs(x = "Standardised coefficient", y = NULL,
       title = "Direction and strength of effects (univariate models)",
       subtitle = "Points: standardised partial coefficients  |  Bars: parametric 95 % CI",
       caption = "Each coefficient is adjusted for the five other predictors; wide intervals reflect collinearity and the small sample size.") +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing = unit(1.2, "lines"))
print(fig_coef)
save_fig(fig_coef, "fig4_univariate_coefficients", 11, 4.5)

# --- Figure 5: importance with jackknife CI -----------------------------------
fig5_data <- uni_importance %>%
  mutate(variable = factor(variable, levels = var_order),
         response = factor(response, levels = resp_levels, labels = pretty_name(resp_levels)))

fig_jack <- ggplot(fig5_data, aes(x = individual, y = variable, fill = group)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_col(width = 0.65) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.25,
                orientation = "y", colour = "grey25") +
  geom_point(aes(x = min_loo), shape = 124, size = 3, colour = "grey45") +
  geom_point(aes(x = max_loo), shape = 124, size = 3, colour = "grey45") +
  facet_wrap(~ response, nrow = 1) +
  scale_fill_manual(values = group_colours, name = NULL) +
  scale_y_discrete(labels = pretty_name) +
  labs(x = "Individual importance (adjusted R²)", y = NULL,
       title = "Predictor importance and sensitivity to individual plots",
       subtitle = "Error bars: jackknife 95 % CI  |  Grey ticks: min and max when each plot is dropped in turn") +
  theme(panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.spacing = unit(1.2, "lines"))
print(fig_jack)
save_fig(fig_jack, "fig5_importance_jackknife", 11, 4.5)

# =============================================================================
# Summary outputs
# =============================================================================
# write.csv(multi_fractions, "multivariate_fractions.csv", row.names = FALSE)
# write.csv(uni_fractions,   "univariate_fractions.csv",   row.names = FALSE)
# write.csv(uni_importance,  "univariate_importance.csv",  row.names = FALSE)
# write.csv(multi_jack,      "multivariate_fractions_jackknife.csv", row.names = FALSE)
# write.csv(uni_frac_jack,   "univariate_fractions_jackknife.csv",   row.names = FALSE)
# write.csv(coef_data,       "univariate_coefficients.csv",          row.names = FALSE)