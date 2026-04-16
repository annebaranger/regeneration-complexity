
library(lidR)
library(lidRviewer)
library(RColorBrewer)
library(data.table)
library(writexl)
library(ggplot2)
library(rgl)
library(terra)
library(sf)
library(raster)
library(dplyr)
library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpmisc)



setwd("C:\\Users\\eck38\\Documents_C_Drive\\Cairngorms_24\\")

files <- c("INSHOG01", "ROTHIEOG01","GLENOG01", "FESHIEOG01","ABYOG01", "FESHIEMT01",  "ABYMT01",  "GLENMT01", "INSHMT01", "ROTHIEMT01","ABYNC01", "ROTHIENC01", "GLENNC01", "INSHNC01", "FESHIENC01", "A2", "A4", "A6", "A5", "C4", "C5", "C7","C8", "plot20", "plot21", "plot31", "plot32", "plot33","plot38","plotnewold", "abyseminat", "plotsnat2"
)

base_path <- "C:/Users/eck38/Documents_C_Drive/Cairngorms_24/AllNorm_Las_return/Las_Normalised_withground"

las_list <- setNames(lapply(files, function(name) {
  readLAS(file.path(base_path, paste0(name, ".las")))
}), files)


resolutions <- c(0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0)

resolutions <- c( 1.9, 1.95, 2.0)


make_chm_list <- function(las_list, res, k = 0.5) {
  cat("Processing resolution:", res, "\n")
  
  l <- lapply(names(las_list), function(name) {
    cat("  Processing plot:", name, "\n")
    
    tryCatch({
      las <- las_list[[name]]
      
      if (is.null(las) || lidR::is.empty(las)) {
        cat("    WARNING: No data\n")
        return(NULL)
      }
      
      subcircle <- k * res
      
      chm <- lidR::rasterize_canopy(
        las, res = res,
        algorithm = lidR::p2r(subcircle = subcircle)
      )
      
      cat("    SUCCESS (subcircle =", subcircle, ")\n")
      chm
      
    }, error = function(e) {
      cat("    ERROR:", conditionMessage(e), "\n")
      NULL
    })
  })
  
  names(l) <- names(las_list)
  l
}


# Create directory
chm_dir <- "C:/Users/eck38/Documents_C_Drive/Cairngorms_24/CHM_files_pitfree"

dir.create(chm_dir, showWarnings = FALSE, recursive = TRUE)
k <- 0.5
# Process and save as GeoTIFF files
for (i in seq_along(resolutions)) {
  res <- resolutions[i]
  
  cat("\n========================================\n")
  cat("Processing", i, "of", length(resolutions), ":", res, "m\n")
  cat("========================================\n")
  
  chm_list <- make_chm_list(las_list, res, k = k)
  
  # Save each CHM as a separate GeoTIFF
  for (plot_name in names(chm_list)) {
    chm <- chm_list[[plot_name]]
    if (is.null(chm)) next
    
    res_str <- gsub("\\.", "_", as.character(res))
    filename <- file.path(chm_dir, paste0("chm_", plot_name, "_res_", res_str, ".tif"))
    
    terra::writeRaster(chm, filename, overwrite = TRUE)
    cat("  Saved:", plot_name, "\n")
  }
  
  rm(chm_list)
  gc()
  
  cat("Completed resolution:", res, "m\n")
}


# Get all CHM files
chm_files <- list.files(chm_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract metadata from filenames
parse_filename <- function(filename) {
  base <- basename(filename)
  # Pattern: chm_PLOTNAME_res_X_XX.tif
  parts <- strsplit(gsub("\\.tif$", "", base), "_")[[1]]
  
  # Find where "res" appears
  res_idx <- which(parts == "res")
  plot_name <- paste(parts[2:(res_idx-1)], collapse = "_")
  res_parts <- parts[(res_idx+1):length(parts)]
  resolution <- as.numeric(paste(res_parts, collapse = "."))
  
  return(list(plot = plot_name, resolution = resolution))}


# Process all files
results_list <- lapply(chm_files, function(file) {
  cat("Processing:", basename(file), "\n")
  
  # Parse filename
  info <- parse_filename(file)
  
  # Load CHM as terra object
  chm_terra <- rast(file)
  
  # Convert to raster package format for lidR
  chm_raster <- raster(chm_terra)
  
  # Extract heights
  heights <- values(chm_terra, na.rm = TRUE)
  
  # Calculate metrics
  rumple <- tryCatch(rumple_index(chm_raster), error = function(e) {
    warning(paste("Rumple error for", info$plot, ":", e$message))
    NA
  })
  
  cv <- if(length(heights) > 0 && mean(heights) > 0) {
    (sd(heights) / mean(heights)) * 100
  } else NA
  
  data.frame(
    resolution = info$resolution,
    plot = info$plot,
    rumple_index = rumple,
    cv_height = cv,
    mean_height = mean(heights, na.rm = TRUE),
    sd_height = sd(heights, na.rm = TRUE),
    max_height = max(heights, na.rm = TRUE),
    n_pixels = length(heights)
  )
})

# Combine all results
results_df <- do.call(rbind, results_list)

# Sort
results_df <- results_df %>%
  arrange(resolution, plot)

category_map <- c(
  "A2"= "OPM" ,
  "A4"= "OPM",
  "A5"= "OPM" ,
  "A6"= "OP", 
  "ABYMT01"= "OP",
  "ROTHIEMT01"= "OP",
  "C4"= "OPM",
  "C5"= "OPM",
  "C7"= "OPM",
  "C8" = "OPM",
  "GLENOG01"= "CP",
  "ROTHIEOG01"= "CP",
  "ABYOG01"= "CP",
  "FESHIEOG01"= "CP",
  "INSHOG01" = "CP",
  "GLENNC01"= "NR",
  "ROTHIENC01"= "NR",
  "ABYNC01"= "NR",
  "FESHIENC01"= "NR",
  "INSHNC01" = "NR",
  "INSHMT01"= "OPM",
  "FESHIEMT01"= "OP",
  "GLENMT01" = "OPM",
  "plot20" = "NR",
  "plot21" = "NR",
  "plot31" = "OP",
  "plot32" = "OPM",
  "plot33" = "OP",
  "plot38" = "OP",
  "plotnewold" = "CP",
  "abyseminat" = "CP", 
  "plotsnat2" = "CP"
)

p2r$Category <- category_map[p2r$plot]


# Save results
write.csv(results_df, "all_chmPITFREE_metrics22.csv", row.names = FALSE)
p2r <- read.csv("all_chm_metrics.csv")

boxplot(rumple_index ~ Category, data = pitfree)

#################### plot for scale differences   ..... without modelling

p2r %>%
  ggplot(aes(x = resolution, y = rumple_log, color = Category)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "~~~~")),
               label.x = 1,
               label.y = c(0.95, 0.90, 0.85, 0.80),
               parse = TRUE, size = 6) +
  labs(x = "CHM resolution", y = "Rumple Index (ln)",
       title = "Rumple Index Sensitivity to CHM(pitfree) resolution") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = margin(t = 20),
        panel.background = element_rect(fill = "lavenderblush", color = NA),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

#####################


####modelling rumple scale effects #####

##### check residuals of data 

shapiro.test(p2r$rumple_index)
qqplot(x = pitfree$resolution, y = qnorm(ppoints(length(pitfree$resolution))))

pitfree$rumple_log <- log(pitfree$rumple_index)
qqnorm(pitfree$rumple_log)
qqline(pitfree$rumple_log)



p2r$rumple_log <- log(p2r$rumple_index)
qqnorm(p2r$rumple_log)
qqline(p2r$rumple_log)






install.packages("statmod")
library(statmod)
pitfree$log_max<- log(pitfree$max_height)

library(MASS)

bc <- boxcox(lm(rumple_index ~ 1, data = pitfree), plotit = TRUE)
lambda <- bc$x[which.max(bc$y)]
lambda
rumple_bc <- (pitfree$rumple_index^lambda - 1) / lambda

qqnorm(canopy_cover$cover_gt_2m)
qqline(canopy_cover$cover_gt_2m)
hist(canopy_cover$cover_gt_2m)
install.packages("betareg")
library(betareg)
beta_model <- betareg(cover_gt_2m ~ canopy_cover$resolution, data = canopy_cover)   




mod <- lm(rumple_log ~ resolution * Category, data = pitfree)
qqnorm(residuals(mod))
qqline(residuals(mod))
shapiro.test(residuals(mod))
print(model)
summary(model)


##### model

modpr2 <- lm(rumple_log ~ resolution * Category, data = p2r)


# mod <- lm(rumple_log ~ resolution * Category, data = pitfree)

# prediction grid
  newdat <- expand.grid(
  resolution = seq(min(p2r$resolution, na.rm = TRUE),
                   max(p2r$resolution, na.rm = TRUE),
                   length.out = 200),
  Category = levels(factor(p2r$Category))
)

# fitted values + 95% CI on the log scale
pr <- predict(modpr2, newdata = newdat, interval = "confidence")

newdat <- cbind(newdat, as.data.frame(pr)) %>%
  rename(fit = fit, lwr = lwr, upr = upr)

ggplot(p2r, aes(x = resolution, y = rumple_log, color = Category)) +
  geom_point(alpha = 0.55, size = 1.6) +
  
  geom_ribbon(
    data = newdat,
    aes(
      x = resolution,
      ymin = lwr,
      ymax = upr,
      fill = Category,
      group = Category
    ),
    alpha = 0.18,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = newdat,
    aes(
      x = resolution,
      y = fit,
      color = Category,
      group = Category
    ),
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "CHM resolution (m)",
    y = "Rumple index (ln)",
    title = "Rumple index sensitivity to p2r CHM resolution\n(linear model with category interaction)"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

summary(mod)

library(scales)
gg_colors <- hue_pal()(4)

# Now use them in your boxplot
m2boxplot <-boxplot(
  sd_height ~ Category,
  data = pitfree[pitfree$resolution == 0.1, ],
  xlab = "Category",
  ylab = "Rumple index",
  main = "Rumple index by category at 2m resolution",
  col = gg_colors
)

#### trying to get statistical significance when resolution is at  


#### predicting using the model doesn't give accurate with the data
predict(mod, newdata = data.frame(resolution = 2, Category = c("CP", "NR", "OP", "OPM")))
exp(c(0.23978504,  0.07996776, -0.02703180,  0.04266631 ))
#### [1] 1.2709759 1.0832521 0.9733303 1.0435896


anVCI <- aov(VCI25cm ~ Forest.type, data = resultz)
summary(anVCI)
plot(anVCI)

shapiro.test(residuals(anVCI))


library(ggplot2)
library(ggpubr)  # For stat_compare_means

# (Optional but recommended) ensure it's a factor
resultz$Forest.type <- as.factor(resultz$Forest.type)

# 1) Welch ANOVA (for subtitle)
welch_VCI<- oneway.test(VCI25cm ~ Forest.type, data = resultz)

Fval <- unname(welch_VCI$statistic)
df1  <- unname(welch_VCI$parameter[1])
df2  <- unname(welch_VCI$parameter[2])
pval <- welch_VCI$p.value
p_text <- ifelse(pval < 0.001, "p < 0.001",
                 paste0("p = ", signif(pval, 3)))

# 2) Games–Howell posthoc
gh <- resultz %>%
  games_howell_test(VCI25cm~ Forest.type)

# Convert GH p-values to compact letters
pvec <- gh %>%
  transmute(pair = paste(group1, group2, sep = "-"),
            p = p.adj) %>%
  deframe()

letters <- multcompLetters(pvec, threshold = 0.05)$Letters

cld_df <- data.frame(
  Forest.type = names(letters),
  letter = letters,
  row.names = NULL
)

# 3) Letter positions (use Rumple_index25cm)
letter_pos <- resultz %>%
  group_by(Forest.type) %>%
  summarise(max_y = max(VCI25cm, na.rm = TRUE), .groups = "drop") %>%
  left_join(cld_df, by = "Forest.type")

# A nice dynamic offset so letters don't collide with points
offset <- 0.05 * diff(range(resultz$VCI25cm, na.rm = TRUE))

# 4) Plot
ggplot(resultz, aes(x = Forest.type, y = VCI25cm, fill = Forest.type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = "black") +
  geom_text(
    data = letter_pos,
    aes(x = Forest.type, y = max_y + offset, label = letter),
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  labs(
    title = "Vertical Complexity Index (25cm)",
    subtitle = bquote(
      "Welch ANOVA: "~
        F[.(sprintf("%.0f", df1))*","*.(sprintf("%.2f", df2))]~"="~
        .(sprintf("%.2f", Fval))*","~
        .(p_text)
    ),
    x = "Forest Category",
    y = "VCI (25cm)"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 16, color = "gray30"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right"
  )
