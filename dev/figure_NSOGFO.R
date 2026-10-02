# method
tar_load(voxnorm_GER14)
load(voxnorm_GER14)
library(ggplot2)
library(patchwork)

# lum : array de dimension 125 x 136 x 185  (x, y, z)
# Exemple pour tester :
# lum <- array(runif(125 * 136 * 185), dim = c(125, 136, 185))

nx <- dim(arr_norm)[1]; ny <- dim(arr_norm)[2]; nz <- dim(arr_norm)[3]

# ---- Paramètres ----
z_couches <- 0:18                        # couches de hauteur à sommer (numérotation à partir de 0)
largeur   <- 18                          # largeur de la tranche en x
x_debut   <- floor((nx - largeur) / 2)   # tranche centrée (à modifier si besoin)
x_tranche <- x_debut:(x_debut + largeur - 1)

# ---- 1) Somme sur les couches z 0 à 18 -> carte (x, y) ----
carte_xy <- rowSums(arr_norm[, , z_couches + 1, drop = FALSE], dims = 2)   # +1 car R indexe à partir de 1
df_xy <- expand.grid(x = 0:(nx - 1), y = 0:(ny - 1))
df_xy$lum <- as.vector(carte_xy)

# ---- 2) Somme sur la tranche en x -> carte (y, z) ----
carte_yz <- colSums(arr_norm[x_tranche + 1, , , drop = FALSE], dims = 1)
df_yz <- expand.grid(y = 0:(ny - 1), z = 0:(nz - 1))
df_yz$lum <- as.vector(carte_yz)

# ---- Échelle commune ----
lim <- range(c(log(df_xy$lum+1), log(df_yz$lum+1)), na.rm = TRUE)

echelle <- scale_fill_viridis_c(name = "Light (logged)", limits = lim)
# En log (valeurs > 0 uniquement) :
# echelle <- scale_fill_viridis_c(name = "Lumière", limits = lim, trans = "log10")

# ---- Graphiques ----
p1 <- ggplot(df_xy, aes(x = x*0.25, y = y*0.25, fill = log(lum+1))) +
  geom_raster() +
  # annotate("rect", xmin = min(x_tranche) - 0.5, xmax = max(x_tranche) + 0.5,
  #          ymin = -0.5, ymax = ny - 0.5,
  #          fill = NA, colour = "white", linetype = "dashed") +
  echelle +
  coord_fixed(expand = FALSE) +
  labs(#title = sprintf("Somme sur z = %d à %d", min(z_couches), max(z_couches)),
       x = "x", y = "y") +
  theme_minimal()

p2 <- ggplot(df_yz, aes(x = y*0.25, y = z*0.25, fill = log(lum+1))) +
  geom_raster() +
  # geom_hline(yintercept = c(min(z_couches) - 0.5, max(z_couches) + 0.5),
  #            colour = "white", linetype = "dashed") +
  echelle +
  coord_fixed(expand = FALSE) +
  labs(#title = sprintf("Tranche x = %d à %d (somme)", min(x_tranche), max(x_tranche)),
       x = "y", y = "z") +
  ylim(c(0,25))+
  theme_minimal()

(p1 + p2) + plot_layout(guides = "collect")

# ggsave("distribution_lumiere.png", (p1 + p2) + plot_layout(guides = "collect"),
#        width = 12, height = 6, dpi = 300)

# Pour sauvegarder :
# ggsave("distribution_lumiere.png", p1 + p2, width = 12, height = 6, dpi = 300)


# first figure
 


ggplot() +
  geom_point(
    data = pts,
    aes(PC1, PC2, shape = plot_level, size = plot_level, color = reg_category)
  ) +
  geom_xsideboxplot(
    data = pts %>% filter(!is.na(reg_category)),
    aes(x = PC1, y = reg_category, fill = reg_category),
    orientation = "y", outliers = FALSE
  ) +
  geom_ysideboxplot(
    data = pts %>% filter(!is.na(reg_category)),
    aes(y = PC2, x = reg_category, fill = reg_category),
    orientation = "x", outliers = FALSE
  )+
  scale_fill_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
  scale_shape_manual(name = "Spatial scale", values = c("Plot" = 19, "Subplot" = 1)) +
  scale_size_manual(name = "Spatial scale", values = c("Plot" = 2.5, "Subplot" = 1.5)) +
  scale_color_brewer(name = "Regeneration cover", palette = "YlGn", na.value = "grey") +
  guides(color="none",size="none",shape="none")+
  new_scale_color() +
  geom_segment(
    data = var_coord,
    aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5, color = type),
    arrow = arrow(length = unit(0.20, "cm")), lineend = "round"
  ) +
  geom_text_repel(
    data = var_coord,
    aes(x = PC1*5, y = PC2*5, label = label, color = type),
    size = 3, show.legend = FALSE
  ) +
  scale_color_brewer(name = "Variable type", palette = "Dark2") +
  theme_minimal() +
  # shrink the side panels + hide their axes
  theme(
    # --- main panel: black frame ---
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    
    # --- side panels: small ---
    ggside.panel.scale = 0.15,          # smaller = thinner strips (try 0.10–0.20)
    
    # --- side panels: no grid ---
    ggside.panel.grid = element_blank(),
    
    # --- side panels: no frame/border of their own, hide their axes ---
    ggside.panel.border = element_blank(),
    ggside.panel.background = element_blank(),
    ggside.axis.text  = element_blank(),
    ggside.axis.ticks = element_blank()  )+
  labs(
    x = paste0("PCA axis 1 (", round(summary(pca_res)$importance["Proportion of Variance","PC1"]*100, 1), "% of variation)"),
    y = paste0("PCA axis 2 (", round(summary(pca_res)$importance["Proportion of Variance","PC2"]*100, 1), "% of variation)")
  ) + 
  coord_cartesian()->c
c

## figure variability
df=metric3d_all_unsliced %>% 
  # filter(subplot=="plot",time==12) %>% 
  dplyr::select(plot,subplot,time,
                rb_mean) %>% 
  left_join(H_ext_all %>% 
              filter(target==0.9) %>% 
              select(plot,subplot,time,slope_ref_agg,ext_ref_agg)) %>% 
  mutate(rb_mean=log(rb_mean)) %>% 
  left_join(gini_all %>% 
              dplyr::select(plot,subplot,time,gini_h_rb)) #%>% 
  # filter(time!=8)



metrics <- c('ext_ref_agg', "slope_ref_agg", "gini_h_rb")#"rb_mean",

# ---- 2. Préparation ----------------------------------------------------------
# On ne garde que les sous-placettes : la ligne "plot" est l'agrégat, elle ne
# peut pas servir à mesurer une variabilité spatiale INTRA-plot.
d <- df %>%
  filter(subplot != "plot") %>%
  mutate(across(all_of(metrics), as.numeric)) %>%
  # sécurité : si plusieurs lignes partagent la même clé plot/subplot/time
  # (c'est le cas dans votre extrait), on les moyenne.
  group_by(plot, subplot, time) %>%
  summarise(across(all_of(metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value")

# ---- 3. Les deux variabilités ------------------------------------------------
# Temporelle : dispersion au cours du temps À L'INTÉRIEUR d'une sous-placette,
#              moyennée sur les sous-placettes.
# Spatiale   : dispersion entre sous-placettes À UNE DATE DONNÉE,
#              moyennée sur les dates.
# (Les versions "marginales" sont aussi calculées : dispersion des moyennes
#  temporelles / des moyennes par sous-placette.)

sd_temporel <- d %>%
  group_by(plot, metric, subplot) %>%
  summarise(sd_t = sd(value, na.rm = TRUE),
            m_sub = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(plot, metric) %>%
  summarise(temporel_within = mean(sd_t, na.rm = TRUE),
            spatial_marg    = sd(m_sub, na.rm = TRUE), .groups = "drop")

sd_spatial <- d %>%
  group_by(plot, metric, time) %>%
  summarise(sd_s = sd(value, na.rm = TRUE),
            m_time = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(plot, metric) %>%
  summarise(spatial_within  = mean(sd_s, na.rm = TRUE),
            temporel_marg   = sd(m_time, na.rm = TRUE), .groups = "drop")

res <- d %>%
  group_by(plot, metric) %>%
  summarise(moyenne = mean(value, na.rm = TRUE),
            n_sub = n_distinct(subplot), n_date = n_distinct(time),
            .groups = "drop") %>%
  left_join(sd_temporel, by = c("plot", "metric")) %>%
  left_join(sd_spatial,  by = c("plot", "metric")) %>%
  mutate(
    # part de la variabilité totale expliquée par le temps (0 = tout spatial,
    # 1 = tout temporel) -> utile pour colorer / trier les plots
    part_temporelle = temporel_within / (temporel_within + spatial_within),
    # versions normalisées (CV) pour comparer des métriques d'échelles
    # différentes ; NA si la moyenne est proche de 0 (cas de slope_ref_agg)
    cv_temporel = ifelse(abs(moyenne) > 1e-8, temporel_within / abs(moyenne), NA),
    cv_spatial  = ifelse(abs(moyenne) > 1e-8, spatial_within  / abs(moyenne), NA)
  )

print(res)

# ---- 4. Figure principale : un point par plot --------------------------------
# x = variabilité spatiale, y = variabilité temporelle.
# La diagonale sépare les plots dominés par le temps (au-dessus) de ceux
# dominés par l'espace (en dessous). Échelles libres car les 3 métriques
# n'ont pas les mêmes unités.

p1 <- res %>% 
  left_join(var_meta,by=c("metric"="variable")) %>% 
  mutate(label=factor(label,levels=c("Mean extinction \n at 4.5m","Gini index","Light profile slope \nat 4.5m"))) %>% 
  ggplot(aes(x = spatial_within, y = temporel_within)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point( size = 2.5) +
  # geom_text(aes(label = plot), vjust = -0.9, size = 3, show.legend = FALSE) +
  scale_colour_gradient2(low = "#2c7bb6", mid = "grey80", high = "#d7191c",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Part\ntemporelle") +
  facet_wrap(~ label,scale="free") +
  expand_limits(x = 0, y = 0) +
  # coord_equal() +          # retirer si les 2 axes ont des ordres de grandeur très différents
  labs(
    # title = "Variabilité temporelle vs spatiale, par plot",
    # subtitle = "Au-dessus de la diagonale : la variabilité temporelle domine",
    x = "Spatial variability (between subplots)",
    y = "Temporal variability (between simulation times)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())

print(p1)

# ---- 5. Figure complémentaire : décomposition en barres ----------------------
# Lecture plus directe quand il y a beaucoup de plots.

res_long <- res %>%
  select(plot, metric, Temporelle = temporel_within, Spatiale = spatial_within) %>%
  pivot_longer(c(Temporelle, Spatiale),
               names_to = "source", values_to = "sd") %>%
  group_by(plot, metric) %>%
  mutate(part = sd / sum(sd, na.rm = TRUE)) %>%
  ungroup()

p2 <- ggplot(res_long, aes(x = reorder(plot, part), y = part, fill = source)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = c(Spatiale = "#2c7bb6", Temporelle = "#d7191c"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ metric) +
  coord_flip() +
  labs(title = "Part relative de chaque source de variabilité",
       x = NULL, y = "Part de la variabilité totale") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())

print(p2)

# ---- 6. Export ---------------------------------------------------------------
# ggsave("variabilite_temporelle_spatiale.png", p1, width = 10, height = 4, dpi = 300)
# ggsave("part_variabilite.png", p2, width = 10, height = 4, dpi = 300)

# ---- Variante ----------------------------------------------------------------
# Pour utiliser les CV (comparables entre métriques, mais instables quand la
# moyenne est proche de 0, ce qui est le cas de slope_ref_agg) :
# ggplot(res, aes(cv_spatial, cv_temporel)) + ... (même structure que p1,
# sans facet_wrap, avec shape = metric)