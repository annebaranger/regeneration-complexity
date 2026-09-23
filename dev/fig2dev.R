
pca_data_plot %>% 
  left_join(metric3d_all_unsliced %>% 
              filter(subplot == "plot", time == 12) %>% 
              dplyr::select(plot, rb_mean, rb_cv)) %>% 
  left_join(H_ext_all %>% 
              filter(target == 0.5, subplot == "plot", time == 12)) %>% 
  left_join(gini_all %>% filter(time == 12) %>% 
              dplyr::select(plot, subplot, gini_h_rb, gini_v_rb,
                            bl_diff_abs, bl_slope_abs)) %>% 
  mutate(rb_mean=log(rb_mean),
         slope_ref_agg=-slope_ref_agg)-> light_test  
light_test$PC1=predict(pca_res,light_test[,var_select])[,1]
# mutate(across(explaining,~log(.)))



explained <- c("richness","richness_sapling","H","cwm_shade","abundance","n_sapling","Height_increment_mn","empty_gf","H_height","fhd","ri","pad_cv")
# explaining <- c("ext_ref_agg", "slope_ref_agg", "rb_cv")#"gini_v_rb",
explaining <- c("rb_mean","gini_h_rb","slope_ref_agg","bl_diff_abs")#"gini_v_rb",

order <- c(var_meta$label[match(explaining, var_meta$variable)],
           var_meta$label[match(explained,  var_meta$variable)])

std_coef <- function(d, y, x) { #, covar = "browsing"
  d <- d %>%
    dplyr::select(all_of(c(y, x))) %>% #, covar
    drop_na() %>%
    mutate(across(where(is.numeric), ~ as.numeric(scale(.))))  # factors (e.g. browsing) left alone
  f <- reformulate(c( x), response = y) #covar,
  unname(coef(lm(f, data = d))[x])
}

# --- bootstrap distribution of that coef ---
boot_coef <- function(d, y, x, B = 1000) { #, covar = "browsing"
  d <- d %>% dplyr::select(all_of(c(y, x))) %>% drop_na() %>% # , covar
    mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) 
  n <- nrow(d)
  f <- reformulate(c( x), response = y) # covar,
  map_dbl(seq_len(B), function(b) {
    db <- d[sample(n, n, replace = TRUE), ]  # <- move this line out of the loop to scale once
    unname(coef(lm(f, data = db))[x])
  })
}

results <- expand_grid(explained = explained, explaining = explaining) %>%
  mutate(
    estimate = map2_dbl(explained, explaining, ~ std_coef(light_test, .x, .y)),
    draws    = map2(explained, explaining, ~ boot_coef(light_test, .x, .y, B = 1000))
  )


boot_long <- results %>%
  dplyr::select(explained, explaining, draws) %>%
  unnest_longer(draws) %>%
  rename(coef = draws)


boot_stats <- boot_long %>%
  group_by(explained, explaining) %>%
  summarise(
    lo   = unname(quantile(coef, 0.10)),
    lo5  = unname(quantile(coef, 0.25)),
    med  = unname(quantile(coef, 0.50)),
    hi5  = unname(quantile(coef, 0.75)),
    hi   = unname(quantile(coef, 0.90)),
    nrep = dplyr::n(),
    # test bilatéral : 2 x la proportion de tirages du "mauvais" côté de 0,
    # bornée par la résolution du bootstrap 1/(B+1)
    p_boot = max(2 * min(mean(coef <= 0), mean(coef >= 0)), 1 / (nrep + 1)),
    .groups = "drop"
  )


boot_stats %>%   filter(explaining !="bl_diff_abs") %>% 
  mutate(
    p_class = cut(p_boot,
                  breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  labels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns"),
                  right = FALSE)
  ) %>%
  pivot_longer(cols = c("explained", "explaining")) %>%
  left_join(var_meta, by = c("value" = "variable")) %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = "name", values_from = c("label", "type")) %>%
  mutate(
    type_explained   = factor(type_explained,
                              levels = c("Composition", "Quantity", "Complexity")),
    label_explained  = factor(label_explained,  levels = rev(order)),
    label_explaining = factor(label_explaining, levels = order)
  ) %>%
  ggplot(aes(x = med, y = label_explained)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, colour = "grey70") +
  geom_linerange(aes(xmin = lo,  xmax = hi,  colour = p_class), linewidth = .9) +
  geom_linerange(aes(xmin = lo5, xmax = hi5, colour = p_class), linewidth = 1.5) +
  geom_point(aes(colour = p_class), size = 1.3) +
  scale_colour_manual(
    values = c("< 0.001" = "#08306b", "< 0.01" = "#2171b5",
               "< 0.05"  = "#6baed6", "< 0.1"  = "#c6dbef",
               "ns"      = "grey78"),
    drop = FALSE, name = "p (bootstrap)"
  ) +
  facet_grid(type_explained ~ label_explaining, scales = "free_y", space = "free_y") +
  labs(x = "Standardized coefficient", y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    panel.spacing.x  = unit(.6, "lines"),
    strip.text.y     = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )



pca_data_subplot %>% 
  left_join(metric3d_all_unsliced %>% 
              filter(subplot!="plot",time==12) %>% 
              dplyr::select(plot,subplot,
                            rb_mean,rb_cv)) %>% 
  left_join(H_ext_all %>% 
              filter(target==0.5,subplot!="plot",time==12)) %>% 
  mutate(rb_mean=log(rb_mean),
         slope_ref_agg=-slope_ref_agg) %>% 
  left_join(gini_all %>% filter(time==12) %>% 
              dplyr::select(plot,subplot,gini_h_rb,gini_v_rb,
                            bl_diff_abs,bl_slope_abs))->light_test_sub

explained <- c("richness","richness_sapling","H","cwm_shade","abundance","n_sapling","Height_increment_mn","empty_gf","H_height","fhd","ri","pad_cv")
# explaining <- c("ext_ref_agg", "slope_ref_agg", "rb_cv")#"gini_v_rb",
explaining <- c("rb_mean","gini_h_rb","slope_ref_agg")#"gini_v_rb",

order <- c(var_meta$label[match(explaining, var_meta$variable)],
           var_meta$label[match(explained,  var_meta$variable)])

order=c(var_meta$label[match(explaining,var_meta$variable)],
        var_meta$label[match(explained,var_meta$variable)])

std_coef <- function(d, y, x) {#, covar = "browsing"
  d <- d %>%
    dplyr::select(all_of(c(y, x, "plot"))) %>%#, covar
    drop_na() %>%
    mutate(plot = factor(plot)) %>%                              # grouping factor, kept out of scaling
    mutate(across(where(is.numeric), ~ as.numeric(scale(.))))    # factors (browsing, plot) left alone
  f <- as.formula(paste0(y, " ~ ",  " + ", x, " + (1 | plot)"))#covar,
  m <- lmer(f, data = d, control = lmerControl(calc.derivs = FALSE))
  unname(fixef(m)[x])
}

# --- bootstrap distribution of that coef ---
boot_coef <- function(d, y, x, B = 1000) {# covar = "browsing",
  d <- d %>% dplyr::select(all_of(c(y, x, "plot"))) %>% drop_na() %>%#covar, 
    mutate(plot = factor(plot)) %>%
    mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
  n <- nrow(d)
  f <- as.formula(paste0(y, " ~ ", " + ", x, " + (1 | plot)"))# covar,
  map_dbl(seq_len(B), function(b) {
    db <- d[sample(n, n, replace = TRUE), ]
    m  <- suppressMessages(suppressWarnings(
      lmer(f, data = db, control = lmerControl(calc.derivs = FALSE))))
    unname(fixef(m)[x])
  })
}


results <- expand_grid(explained = explained, explaining = explaining) %>%
  mutate(
    estimate = map2_dbl(explained, explaining, ~ std_coef(light_test_sub, .x, .y)),
    draws    = map2(explained, explaining, ~ boot_coef(light_test_sub, .x, .y, B = 1000))
  )


boot_long <- results %>%
  dplyr::select(explained, explaining, draws) %>%
  unnest_longer(draws) %>%
  rename(coef = draws)


boot_stats <- boot_long %>%
  group_by(explained, explaining) %>%
  summarise(
    lo   = unname(quantile(coef, 0.10)),
    lo5  = unname(quantile(coef, 0.25)),
    med  = unname(quantile(coef, 0.50)),
    hi5  = unname(quantile(coef, 0.75)),
    hi   = unname(quantile(coef, 0.90)),
    nrep = dplyr::n(),
    # test bilatéral : 2 x la proportion de tirages du "mauvais" côté de 0,
    # bornée par la résolution du bootstrap 1/(B+1)
    p_boot = max(2 * min(mean(coef <= 0), mean(coef >= 0)), 1 / (nrep + 1)),
    .groups = "drop"
  )

boot_stats_sub=boot_stats
boot_stats %>%
  filter(explaining !="bl_diff_abs") %>% 
  mutate(
    p_class = cut(p_boot,
                  breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  labels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns"),
                  right = FALSE)
  ) %>%
  pivot_longer(cols = c("explained", "explaining")) %>%
  left_join(var_meta, by = c("value" = "variable")) %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = "name", values_from = c("label", "type")) %>%
  mutate(
    type_explained   = factor(type_explained,
                              levels = c("Composition", "Quantity", "Complexity")),
    label_explained  = factor(label_explained,  levels = rev(order)),
    label_explaining = factor(label_explaining, levels = order)
  ) %>%
  ggplot(aes(x = med, y = label_explained)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, colour = "grey70") +
  geom_linerange(aes(xmin = lo,  xmax = hi,  colour = p_class), linewidth = .9) +
  geom_linerange(aes(xmin = lo5, xmax = hi5, colour = p_class), linewidth = 1.5) +
  geom_point(aes(colour = p_class), size = 1.3) +
  scale_colour_manual(
    values = c("< 0.001" = "#08306b", "< 0.01" = "#2171b5",
               "< 0.05"  = "#6baed6", "< 0.1"  = "#c6dbef",
               "ns"      = "grey78"),
    drop = FALSE, name = "p (bootstrap)"
  ) +
  facet_grid(type_explained ~ label_explaining, scales = "free_y", space = "free_y") +
  labs(x = "Standardized coefficient", y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    panel.spacing.x  = unit(.6, "lines"),
    strip.text.y     = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )


pd <- position_dodge(width = 0.6)   # écart vertical entre "plot" et "subplot"

rbind(boot_stats %>% mutate(level = "plot"),
      boot_stats_sub %>% mutate(level = "subplot")) %>%
  filter(explaining != "bl_diff_abs") %>%
  filter(!explained%in%c("H_height","empty_gf") ) %>% 
  mutate(
    p_class = cut(p_boot,
                  breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  labels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns"),
                  right = FALSE)
  ) %>%
  pivot_longer(cols = c("explained", "explaining")) %>%
  left_join(var_meta, by = c("value" = "variable")) %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = "name", values_from = c("label", "type")) %>%
  mutate(
    type_explained   = factor(type_explained,
                              levels = c("Composition", "Quantity", "Complexity")),
    label_explained  = factor(label_explained,  levels = rev(order)),
    label_explaining = factor(label_explaining, levels = order),
    level            = factor(level, levels = c("subplot", "plot"))  # "plot" au-dessus
  ) %>%
  ggplot(aes(x = med, y = label_explained, group = level, alpha = level)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .3, colour = "grey70") +
  # IC 95 %
  geom_linerange(aes(xmin = lo, xmax = hi, colour = p_class),
                 linewidth = .9, position = pd) +
  # IC 50 % + médiane (point)
  geom_pointrange(aes(xmin = lo5, xmax = hi5, colour = p_class),
                  linewidth = 1.5, size = .3, position = pd) +
  scale_colour_manual(
    values = c("< 0.001" = "#08306b", "< 0.01" = "#2171b5",
               "< 0.05"  = "#6baed6", "< 0.1"  = "#c6dbef",
               "ns"      = "grey78"),
    drop = FALSE, name = "p (bootstrap)"
  ) +
  scale_alpha_manual(values = c(plot = 0.2, subplot = 1), name = "Level") +
  facet_grid(type_explained ~ label_explaining, scales = "free_y", space = "free_y") +
  labs(x = "Standardized coefficient", y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    panel.spacing.x  = unit(.6, "lines"),
    strip.text.y     = element_text(angle = 0),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )  
