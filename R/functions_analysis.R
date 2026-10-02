#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analysis.R  
#' @description R script containing helpers functions
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- PCA helpers ---------------------------------------------------------------
pca_var_coord <- function(pca, var_meta, lidar_vars = NULL) {
  v <- factoextra::get_pca_var(pca)
  tibble(
    variable = rownames(v$coord),
    PC1      = v$coord[, 1],
    PC2      = v$coord[, 2],
    contrib  = rowSums(v$contrib[, 1:2]),
    cos2     = rowSums(v$cos2[, 1:2])
  ) %>%
    left_join(var_meta, by = "variable") %>%
    mutate(label  = fix_cross_labels(label, variable),
           type   = tidyr::replace_na(type, "Other"),
           Source = if_else(variable %in% lidar_vars, "Lidar-derived", "Inventory-derived"))
}

# --- Regeneration structure PCA ------------------------------------------------
run_pca_regen <- function(tls_metrics, forest_data_plot, forest_data_subplot,
                          regen_cat, vars, lidar_vars, var_meta) {
  d_plot <- prep_data(tls_metrics, forest_data_plot,    vars, plot_level = "plot")
  d_sub  <- prep_data(tls_metrics, forest_data_subplot, vars, plot_level = "subplot")
  
  X <- d_plot %>%
    dplyr::select(-plot, -subplot) %>%
    mutate(across(everything(), as.numeric)) %>%
    tidyr::drop_na()
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  
  all_d  <- bind_rows(d_plot, d_sub) %>% left_join(regen_cat)
  scores <- bind_cols(all_d, as_tibble(predict(pca, all_d[, colnames(X)])[, 1:2]))
  
  kw <- scores %>%
    filter(subplot == "plot", !is.na(reg_category)) %>%
    mutate(reg_category = droplevels(factor(reg_category))) %>%
    tidyr::pivot_longer(c(PC1, PC2), names_to = "axis", values_to = "score") %>%
    group_by(axis) %>%
    rstatix::kruskal_test(score ~ reg_category) %>%
    rstatix::add_significance("p") %>%
    mutate(label = paste0(
      "KW :\n ",
      ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3))),
      " ", ifelse(p.signif == "ns", "", p.signif)))
  
  list(pca = pca, data = X, scores = scores,
       var_coord = pca_var_coord(pca, var_meta, lidar_vars), kw = kw)
}

# --- Light PCA (plot or subplot) -----------------------------------------------
run_pca_light <- function(tls_metrics, forest_data, vars, plot_level, var_meta) {
  d <- prep_data(tls_metrics, forest_data, c(vars, "reg_category"),
                 plot_level = plot_level, t = 12, s = 0)
  X <- dplyr::select(d, all_of(vars))
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  list(pca = pca, data = X,
       scores = bind_cols(d, as_tibble(predict(pca, X)[, 1:2])),
       var_coord = pca_var_coord(pca, var_meta))
}

# --- Tests ---------------------------------------------------------------------
cor_table <- function(df) {
  purrr::map_dfr(combn(names(df), 2, simplify = FALSE), function(v) {
    ct <- cor.test(df[[v[1]]], df[[v[2]]], method = "pearson", exact = FALSE)
    tibble(var1 = v[1], var2 = v[2], rho = unname(ct$estimate), p = ct$p.value)
  }) %>%
    mutate(p_adj  = p.adjust(p, method = "bonferroni"),
           signif = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                              p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
    arrange(desc(abs(rho)))
}

kw_by_variable <- function(box_df) {
  box_df %>%
    group_by(type, label) %>%
    summarise(
      p = if (n_distinct(reg_category) < 2) NA_real_
      else kruskal.test(value ~ droplevels(reg_category))$p.value,
      y = max(value, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(stars = p_stars(p))
}

cld_by_variable <- function(box_df) {
  cld_one <- function(df) {
    lev <- levels(droplevels(df$reg_category))
    if (length(lev) < 2) return(tibble(reg_category = lev, letter = ""))
    kw <- kruskal.test(value ~ droplevels(reg_category), data = df)
    if (is.na(kw$p.value) || kw$p.value > 0.05)
      return(tibble(reg_category = lev, letter = "a"))
    dt <- rstatix::dunn_test(df, value ~ reg_category, p.adjust.method = "BH")
    l  <- multcompView::multcompLetters(
      setNames(dt$p.adj, paste(dt$group1, dt$group2, sep = "-")))$Letters
    tibble(reg_category = names(l), letter = unname(l))
  }
  box_df %>%
    group_by(type, label) %>%
    group_modify(~ cld_one(.x)) %>%
    ungroup() %>%
    mutate(reg_category = factor(reg_category, levels = levels(box_df$reg_category))) %>%
    left_join(box_df %>% group_by(type, label, reg_category) %>%
                summarise(y = max(value, na.rm = TRUE), .groups = "drop"),
              by = c("type", "label", "reg_category"))
}

# --- Bootstrapped standardized effects -----------------------------------------
# One model per (response, predictor) pair: lm, or lmer with (1 | random).
# cluster_boot = TRUE resamples whole plots instead of rows (recommended with lmer).
boot_effects <- function(d, responses, predictors, B = 1000,
                         random = NULL, cluster_boot = FALSE) {
  fit <- function(dd, y, x) {
    if (is.null(random)) {
      unname(coef(lm(reformulate(x, y), data = dd))[2])
    } else {
      f <- as.formula(paste0(y, " ~ ", x, " + (1 | ", random, ")"))
      m <- suppressMessages(suppressWarnings(
        lme4::lmer(f, data = dd, control = lme4::lmerControl(calc.derivs = FALSE))))
      unname(lme4::fixef(m)[2])
    }
  }
  resample <- function(dd) {
    if (!cluster_boot || is.null(random)) return(dd[sample(nrow(dd), replace = TRUE), ])
    ids  <- unique(as.character(dd[[random]]))
    pick <- sample(ids, length(ids), replace = TRUE)
    out  <- bind_rows(lapply(seq_along(pick), function(i) {
      o <- dd[as.character(dd[[random]]) == pick[i], ]
      o[[random]] <- paste(pick[i], i, sep = "_")   # duplicated plots = distinct groups
      o
    }))
    out[[random]] <- factor(out[[random]])
    out
  }
  
  tidyr::expand_grid(explained = responses, explaining = predictors) %>%
    purrr::pmap_dfr(function(explained, explaining) {
      dd <- d %>%
        dplyr::select(all_of(c(explained, explaining, random))) %>%
        tidyr::drop_na() %>%
        mutate(across(any_of(random), ~ factor(.x)),
               across(where(is.numeric), ~ as.numeric(scale(.x))))
      est   <- fit(dd, explained, explaining)
      draws <- vapply(seq_len(B), function(b)
        tryCatch(fit(resample(dd), explained, explaining), error = function(e) NA_real_),
        numeric(1))
      q    <- quantile(draws, c(.10, .25, .50, .75, .90), na.rm = TRUE, names = FALSE)
      nrep <- sum(!is.na(draws))
      tibble(explained, explaining, estimate = est,
             lo = q[1], lo5 = q[2], med = q[3], hi5 = q[4], hi = q[5], nrep = nrep,
             # two-sided: 2 x share of draws on the "wrong" side of 0, floored at 1/(B+1)
             p_boot = max(2 * min(mean(draws <= 0, na.rm = TRUE),
                                  mean(draws >= 0, na.rm = TRUE)), 1 / (nrep + 1)))
    })
}

# --- Temporal vs spatial variability -------------------------------------------
compute_light_variability <- function(tls_metrics, vars) {
  lv <- tls_metrics %>%
    filter(slice_num == 0, subplot != "plot") %>%
    dplyr::select(any_of(c("plot", "subplot", "time", vars))) %>%
    unique() %>%
    mutate(across(all_of(vars), as.numeric)) %>%
    tidyr::pivot_longer(all_of(vars), names_to = "metric", values_to = "value")
  
  sd_temporal <- lv %>%
    group_by(plot, metric, subplot) %>%
    summarise(sd_t = sd(value, na.rm = TRUE), m_sub = mean(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, metric) %>%
    summarise(temporal_within = mean(sd_t, na.rm = TRUE),
              spatial_marg    = sd(m_sub, na.rm = TRUE), .groups = "drop")
  
  sd_spatial <- lv %>%
    group_by(plot, metric, time) %>%
    summarise(sd_s = sd(value, na.rm = TRUE), m_time = mean(value, na.rm = TRUE), .groups = "drop") %>%
    group_by(plot, metric) %>%
    summarise(spatial_within = mean(sd_s, na.rm = TRUE),
              temporal_marg  = sd(m_time, na.rm = TRUE), .groups = "drop")
  
  lv %>%
    group_by(plot, metric) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              n_sub = n_distinct(subplot), n_date = n_distinct(time), .groups = "drop") %>%
    left_join(sd_temporal, by = c("plot", "metric")) %>%
    left_join(sd_spatial,  by = c("plot", "metric")) %>%
    mutate(temporal_part = temporal_within / (temporal_within + spatial_within),
           cv_temporal = ifelse(abs(mean) > 1e-8, temporal_within / abs(mean), NA),
           cv_spatial  = ifelse(abs(mean) > 1e-8, spatial_within  / abs(mean), NA))
}

# --- Variance partitioning (light vs adults) ------------------------------------
run_varpart <- function(tls_metrics, forest_data_plot, adult_vars, light_vars, regen_vars,
                        log_vars = character(0), n_perm = 999, n_perm_hp = 999) {
  dat <- prep_data(tls_metrics, forest_data_plot,
                   id_cols = c(regen_vars, adult_vars, light_vars), plot_level = "plot")
  excluded <- dat$plot[!complete.cases(dat)]
  dat <- dat[complete.cases(dat), ]
  n   <- nrow(dat)
  
  standardise <- function(d) mutate(d, across(everything(), ~ as.numeric(scale(.x))))
  regen  <- dat %>% dplyr::select(all_of(regen_vars)) %>%
    mutate(across(all_of(log_vars), log1p)) %>% standardise()
  adults <- dat %>% dplyr::select(all_of(adult_vars)) %>% standardise()
  light  <- dat %>% dplyr::select(all_of(light_vars)) %>% standardise()
  predictors <- cbind(adults, light)
  vif <- diag(solve(cor(predictors)))
  
  perm_p <- function(m) anova(m, permutations = n_perm)$`Pr(>F)`[1]
  
  partition <- function(Y) {
    f <- vegan::varpart(Y, light, adults)$part$fract$Adj.R.square  # light, adults, both
    data.frame(
      full_adj_r2    = f[3],
      p_full         = perm_p(vegan::rda(Y ~ ., data = predictors)),
      light_total    = f[1],
      p_light_total  = perm_p(vegan::rda(Y ~ ., data = light)),
      adult_total    = f[2],
      p_adult_total  = perm_p(vegan::rda(Y ~ ., data = adults)),
      light_unique   = f[3] - f[2],
      p_light_unique = perm_p(vegan::rda(Y ~ . + Condition(as.matrix(adults)), data = light)),
      shared         = f[1] + f[2] - f[3],
      adult_unique   = f[3] - f[1],
      p_adult_unique = perm_p(vegan::rda(Y ~ . + Condition(as.matrix(light)), data = adults)),
      residuals      = 1 - f[3]
    )
  }
  
  # A. multivariate
  model       <- vegan::rda(regen ~ ., data = predictors)
  global_test <- anova(model, permutations = n_perm)
  axis_test   <- anova(model, by = "axis", permutations = n_perm)
  multi       <- partition(regen)
  multi_fractions <- data.frame(
    fraction = c("Full model [a+b+c]", "Light (total) [a+b]", "Adults (total) [b+c]",
                 "Light (unique) [a]", "Shared [b]", "Adults (unique) [c]", "Residuals [d]"),
    adj_r2   = round(c(multi$full_adj_r2, multi$light_total, multi$adult_total,
                       multi$light_unique, multi$shared, multi$adult_unique, multi$residuals), 3),
    p_value  = c(multi$p_full, multi$p_light_total, multi$p_adult_total,
                 multi$p_light_unique, NA, multi$p_adult_unique, NA))
  
  # B. univariate
  uni <- lapply(regen_vars, function(v) {
    y     <- regen[[v]]
    hp    <- rdacca.hp::rdacca.hp(y, predictors, method = "RDA", type = "adjR2")$Hier.part
    perm  <- rdacca.hp::permu.hp(y, predictors, method = "RDA", type = "adjR2",
                                 permutations = n_perm_hp)
    p_col <- grep("^Pr", colnames(perm), value = TRUE)[1]
    list(
      frac = cbind(response = v, partition(y)),
      imp  = data.frame(
        response = v, variable = rownames(hp),
        group = ifelse(rownames(hp) %in% adult_vars, "Adults", "Light"),
        unique = hp[, "Unique"], shared = hp[, "Average.share"],
        individual = hp[, "Individual"], percent = hp[, "I.perc(%)"],
        p_perm = as.numeric(gsub("[^0-9.eE-]", "", perm[rownames(hp), p_col]))))
  })
  uni_fractions <- bind_rows(lapply(uni, `[[`, "frac")) %>%
    mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "holm"), .names = "{.col}_holm"))
  uni_importance <- bind_rows(lapply(uni, `[[`, "imp")) %>%
    group_by(response) %>%
    mutate(p_holm = p.adjust(p_perm, method = "holm"),
           rank   = rank(-individual, ties.method = "min")) %>%
    ungroup()
  
  # C. jackknife (leave-one-plot-out)
  drop_k <- function(Y, k) if (is.null(dim(Y))) Y[-k] else Y[-k, , drop = FALSE]
  jk_fractions <- function(Y) {
    sapply(seq_len(n), function(k) {
      Yk <- drop_k(Y, k)
      r2 <- function(X) vegan::RsquareAdj(vegan::rda(Yk ~ ., data = X[-k, , drop = FALSE]))$adj.r.squared
      f_l <- r2(light); f_a <- r2(adults); f_f <- r2(predictors)
      c(light_unique = f_f - f_a, shared = f_l + f_a - f_f,
        adult_unique = f_f - f_l, full_adj_r2 = f_f)
    })
  }
  jk_summary <- function(jk, estimate) {
    est <- estimate[rownames(jk)]
    se  <- sqrt((n - 1) / n * rowSums((jk - rowMeans(jk))^2))
    data.frame(term = rownames(jk), estimate = est, se_jack = se,
               ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
               min_loo = apply(jk, 1, min), max_loo = apply(jk, 1, max),
               most_influential_plot = dat$plot[apply(abs(jk - est), 1, which.max)],
               row.names = NULL)
  }
  keep <- c("light_unique", "shared", "adult_unique", "full_adj_r2")
  multi_jack <- jk_summary(jk_fractions(regen), unlist(multi[keep]))
  uni_frac_jack <- bind_rows(lapply(regen_vars, function(v) {
    est <- unlist(uni_fractions[uni_fractions$response == v, keep])
    cbind(response = v, jk_summary(jk_fractions(regen[[v]]), est))
  }))
  uni_imp_jack <- bind_rows(lapply(regen_vars, function(v) {
    y  <- regen[[v]]
    jk <- sapply(seq_len(n), function(k)
      rdacca.hp::rdacca.hp(y[-k], predictors[-k, ], method = "RDA",
                           type = "adjR2")$Hier.part[, "Individual"])
    imp <- uni_importance[uni_importance$response == v, ]
    cbind(response = v, jk_summary(jk, setNames(imp$individual, imp$variable)))
  }))
  uni_importance <- uni_importance %>%
    left_join(uni_imp_jack %>%
                dplyr::select(response, variable = term, se_jack, ci_low, ci_high,
                              min_loo, max_loo, most_influential_plot),
              by = c("response", "variable")) %>%
    arrange(response, rank)
  
  list(n = n, excluded = excluded, vif = vif,
       model = model, model_r2 = vegan::RsquareAdj(model),
       global_test = global_test, axis_test = axis_test,
       multi = multi, multi_fractions = multi_fractions, multi_jack = multi_jack,
       uni_fractions = uni_fractions, uni_frac_jack = uni_frac_jack,
       uni_importance = uni_importance,
       adult_vars = adult_vars, light_vars = light_vars, regen_vars = regen_vars)
}
