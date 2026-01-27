library(VoxR)

# Setup
input_folder <- "C:/Users/eck38/Documents_C_Drive/Cairngorms_24/AllNorm_Las_return/No_ground/ASC/"
output_csv <- "C:/Users/eck38/Documents_C_Drive/Cairngorms_24/enl_results3.csv"

# Get all files
all_files <- list.files(input_folder, pattern = "\\.asc$", full.names = TRUE)
cat("Found", length(all_files), "files to process\n\n")

##function found at https://rdrr.io/github/spatial-mk/tre3d/src/R/ENL.R 
##author Matthias Kunz, last updated: 27.04.2019

ENL <- function(cloud_xyz,
                voxel_size = 0.20,     # 20 cm voxels (paper standard)
                layer_thickness = 1.0  # 1 m height layers
)
{
  cat("\nComputing ENL (Ehbrecht et al. 2016 method)\n")
  
  #----------------------------------------------------------
  # 1. Height normalisation check
  #----------------------------------------------------------
  if (min(cloud_xyz$Z) < -0.001) {
    stop("Z values must be height-normalised (height above ground).")
  }
  
  #----------------------------------------------------------
  # 2. Voxelisation (20 cm grid)
  #    VoxR::vox returns 1 point per filled voxel → perfect for ENL
  #----------------------------------------------------------
  vox <- VoxR::vox(cloud_xyz, res = voxel_size)
  names(vox)[1:3] <- c("x", "y", "z")
  
  #----------------------------------------------------------
  # 3. Convert voxel heights to layer indices
  #----------------------------------------------------------
  vox$layer <- floor(vox$z / layer_thickness)
  
  # number of layers from 0 to max
  n_layers <- max(vox$layer, na.rm = TRUE) + 1
  
  #----------------------------------------------------------
  # 4. Count filled voxels per 1-m vertical layer
  #----------------------------------------------------------
  voxel_counts <- table(factor(vox$layer, levels = 0:(n_layers - 1)))
  
  # convert to numeric vector
  Fi <- as.numeric(voxel_counts)
  
  # total number of filled voxels
  F_total <- sum(Fi)
  
  # proportions per layer
  pi <- Fi / F_total
  
  #----------------------------------------------------------
  # 5. ENL computation via Hill numbers
  #----------------------------------------------------------
  
  # ENL_0 (0D) = number of layers with >0 filled voxels
  ENL0 <- sum(Fi > 0)
  
  # ENL_1 (1D) = exponential Shannon
  ENL1 <- exp(- sum(pi[pi > 0] * log(pi[pi > 0])))
  
  # ENL_2 (2D) = inverse Simpson
  ENL2 <- 1 / sum(pi^2)
  
  #----------------------------------------------------------
  # 6. Output
  #----------------------------------------------------------
  data.frame(
    ENL_0 = ENL0,
    ENL_1 = ENL1,
    ENL_2 = ENL2,
    voxel_size = voxel_size,
    layer_thickness = layer_thickness,
    n_voxels = F_total,
    n_layers = n_layers
  )
}

# results data frame
ENL_results3 <- data.frame(
  filename = character(),
  n_points = integer(),
  ENL_0 = numeric(),
  ENL_1 = numeric(),
  ENL_2 = numeric(),
  voxel_size = numeric(),
  layer_thickness = numeric(),
  n_voxels = integer(),
  n_layers = integer(),
  error_message = character(),
  stringsAsFactors = FALSE
)

# Loop through files
for (i in seq_along(all_files)) {
  file <- all_files[i]
  filename <- basename(file)
  
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(all_files), filename))
  
  tryCatch({
    # Read first few lines to check format
    test_lines <- readLines(file, n = 5)
    cat("  First line preview:", substr(test_lines[1], 1, 100), "\n")
    
    # Try to read the cloud - skip any header lines
    cloud <- read.table(file, header = FALSE, skip = 0, 
                        colClasses = c("numeric", "numeric", "numeric"),
                        comment.char = "")
    
    # If that fails, try skipping first line
    if (nrow(cloud) == 0 || any(is.na(cloud[1, 1:3]))) {
      cat("  Retrying with skip=1...\n")
      cloud <- read.table(file, header = FALSE, skip = 1,
                          colClasses = c("numeric", "numeric", "numeric"),
                          comment.char = "")
    }
    
    # Check number of columns
    cat("  Columns in file:", ncol(cloud), "\n")
    
    if (ncol(cloud) < 3) {
      stop("File has fewer than 3 columns")
    }
    
    # Set column names
    colnames(cloud)[1:3] <- c("X", "Y", "Z")
    
    # Keep only XYZ
    cloud_xyz <- cloud[, 1:3]
    colnames(cloud_xyz) <- c("X", "Y", "Z")
    
    # Force numeric conversion
    cloud_xyz$X <- as.numeric(cloud_xyz$X)
    cloud_xyz$Y <- as.numeric(cloud_xyz$Y)
    cloud_xyz$Z <- as.numeric(cloud_xyz$Z)
    
    # Check for issues
    cat("  Points loaded:", nrow(cloud_xyz), "\n")
    
    # Check for NA or infinite values AFTER conversion
    na_count <- sum(is.na(cloud_xyz))
    if (na_count > 0) {
      cat("  WARNING: Found", na_count, "NA values - removing them\n")
      cloud_xyz <- na.omit(cloud_xyz)
    }
    
    if (any(!is.finite(as.matrix(cloud_xyz)))) {
      stop("Data contains infinite values")
    }
    
    if (nrow(cloud_xyz) == 0) {
      stop("No valid points remaining after cleaning")
    }
    
    cat("  Z range:", min(cloud_xyz$Z), "to", max(cloud_xyz$Z), "\n")
    
    # Check if there's height variation
    z_range <- max(cloud_xyz$Z) - min(cloud_xyz$Z)
    if (z_range < 0.01) {
      stop(paste("No height variation - Z range is only", round(z_range, 4), "m"))
    }
    
    cat("  Height range:", round(z_range, 2), "m\n")
    
    # Run ENL function
    output <- ENL(cloud_xyz)
    
    # Add results - FIXED to match output structure
    ENL_results3 <- rbind(ENL_results3, data.frame(
      filename = filename,
      n_points = nrow(cloud_xyz),
      ENL_0 = output$ENL_0,
      ENL_1 = output$ENL_1,
      ENL_2 = output$ENL_2,
      voxel_size = output$voxel_size,
      layer_thickness = output$layer_thickness,
      n_voxels = output$n_voxels,
      n_layers = output$n_layers,
      error_message = "",
      stringsAsFactors = FALSE
    ))
    
    # FIXED: Print correct ENL values
    cat(sprintf("  ✅ ENL_0: %.0f, ENL_1: %.3f, ENL_2: %.3f\n\n", 
                output$ENL_0, output$ENL_1, output$ENL_2))
    
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    cat("  ❌ ERROR:", error_msg, "\n\n")
    
    # FIXED: Add failed row with correct structure and variable name
    ENL_results2 <<- rbind(ENL_results2, data.frame(
      filename = filename,
      n_points = NA,
      ENL_0 = NA,
      ENL_1 = NA,
      ENL_2 = NA,
      voxel_size = NA,
      layer_thickness = NA,
      n_voxels = NA,
      n_layers = NA,
      error_message = error_msg,
      stringsAsFactors = FALSE
    ))
  })
}

# Save results
write.csv(ENL_results3, output_csv, row.names = FALSE)



ENL<-read.csv("C:\\Users\\eck38\\OneDrive - University of Cambridge\\Documents\\Analysis\\Results\\Official\\enl_results3.csv")

ENL$filename <- sub("_.*","", ENL$filename)

# Run ANOVA of box dimension results 
modelENL <- aov(ENL_2 ~ Forest.type, data = resultz)
anova_sum <- summary(modelENL)
summary(modelENL)
###check model for normality assumptions
par(mfrow = c(2, 2))
plot(modelENL)

shapiro.test(residuals(modelBD))  ###data:  residuals(modelENL) W = 0.92071, p-value = 0.02173. Unfortunately the pvalueis below 0.05 which 
### means that it is like not normaly distributed. The QQ plots show there is likely heteroskedacity so will need to transform

ENLlog <- aov(log(ENL_2) ~Forest.type, data=resultz)
anova_sum<-summary(ENLlog)
summary(ENLlog)
shapiro.test(residuals(ENLlog))
plot(ENLlog)
any(resultz$ENL_2 <= 0)

welch_ENL <- oneway.test(ENL_2 ~ Forest.type, data = resultz)
welch_ENL
library(rstatix)

resultz %>%
  games_howell_test(ENL_2 ~ Forest.type)



# Pairwise Games–Howell p-values
gh <- resultz %>%
  games_howell_test(ENL_2 ~ Forest.type) %>%
  mutate(p.adj = p.adj)  # already there; keeps naming explicit

# Make a named vector of adjusted p-values like "A-B" = p
pvec <- gh %>%
  transmute(pair = paste(group1, group2, sep = "-"),
            p = p.adj) %>%
  deframe()

# Convert p-values to letters (threshold can be changed, e.g. 0.01)
letters <- multcompLetters(pvec, threshold = 0.05)$Letters

cld_df <- tibble(
  Forest.type = names(letters),
  letter = unname(letters)
)
letter_pos <- resultz %>%
  group_by(Forest.type) %>%
  summarise(max_y = max(ENL_2, na.rm = TRUE), .groups = "drop") %>%
  left_join(cld_df, by = "Forest.type")
library(ggplot2)

ggplot(resultz, aes(x = Forest.type, y = ENL_2, fill = Forest.type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.4, color = "black") +
  geom_text(
    data = letter_pos,
    aes(x = Forest.type, y = max_y + 0.05, label = letter),
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  labs(
    title = "ENL",
    x = "Forest Category",
    y = "ENL_2"
  ) +
  theme_minimal()
