# ============================================================
# TLS regeneration complexity metrics from a LAS/LAZ point cloud
# ============================================================

# ---- Packages ----
req <- c("lidR", "terra", "dplyr", "data.table", "spdep", "gstat", "sf")
to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(req, library, character.only = TRUE)

# ---- User parameters ----
las_path <- "path/to/your_pointcloud.laz"   # <-- change

# regeneration height band (meters above ground)
zmin_regen <- 0.2
zmax_regen <- 3.0

# raster grid resolution for horizontal metrics (meters)
grid_res <- 0.5

# voxel size for occupancy-based metrics (meters)
voxel_xy <- 0.25
voxel_z  <- 0.20

# height bin size for vertical profiles (meters)
zbin <- 0.20

# terrain normalization settings (if Z already normalized, set normalize = FALSE)
normalize <- TRUE
dtm_res <- 0.5     # DTM resolution for normalization
csf_cloth_resolution <- 0.5

# lacunarity settings (box sizes in cells on the regen presence raster)
lac_box_sizes <- c(2, 3, 5, 8, 13, 21)  # in raster cells (not meters)

# variogram settings
variogram_cutoff <- 15     # meters
variogram_width  <- 1.0    # meters

set.seed(1)

# ============================================================
# Helper functions
# ============================================================

# Safe Shannon entropy (natural log)
shannon_entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

# Foliage Height Diversity (FHD): Shannon entropy of normalized vertical profile
compute_fhd <- function(z, bin = 0.2, zmin = NULL, zmax = NULL) {
  if (!is.null(zmin)) z <- z[z >= zmin]
  if (!is.null(zmax)) z <- z[z <= zmax]
  if (length(z) < 10) return(NA_real_)
  brks <- seq(floor(min(z)/bin)*bin, ceiling(max(z)/bin)*bin, by = bin)
  h <- hist(z, breaks = brks, plot = FALSE)$counts
  if (sum(h) == 0) return(NA_real_)
  p <- h / sum(h)
  shannon_entropy(p)
}

# Gini coefficient (simple)
gini <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  if (sum(x) == 0) return(0)
  (2 * sum(x * seq_len(n)) / (n * sum(x))) - ((n + 1) / n)
}

# Coefficient of variation
cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  m <- mean(x)
  if (m == 0) return(NA_real_)
  sd(x) / m
}

# Moran's I from a raster with NA handling
moran_I_raster <- function(r, queen = TRUE) {
  # r: terra SpatRaster
  m <- as.matrix(r, wide = TRUE)
  # build neighbor list for all cells
  nr <- nrow(m); nc <- ncol(m)
  idx <- which(!is.na(m), arr.ind = TRUE)
  if (nrow(idx) < 10) return(list(I = NA_real_, p = NA_real_))
  
  # map matrix coords -> 1..k
  key <- matrix(NA_integer_, nr, nc)
  key[!is.na(m)] <- seq_len(sum(!is.na(m)))
  
  # neighbors for each valid cell
  nb <- vector("list", length = sum(!is.na(m)))
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]; j <- idx[k, 2]
    # neighbor offsets
    offs <- if (queen) expand.grid(di = -1:1, dj = -1:1) else rbind(
      c(-1,0), c(1,0), c(0,-1), c(0,1)
    )
    if (queen) offs <- offs[!(offs$di == 0 & offs$dj == 0), , drop = FALSE]
    neigh <- lapply(seq_len(nrow(offs)), function(t) c(i + offs[t,1], j + offs[t,2]))
    neigh <- do.call(rbind, neigh)
    ok <- neigh[,1] >= 1 & neigh[,1] <= nr & neigh[,2] >= 1 & neigh[,2] <= nc
    neigh <- neigh[ok, , drop = FALSE]
    nids <- key[cbind(neigh[,1], neigh[,2])]
    nids <- nids[!is.na(nids)]
    nb[[ key[i,j] ]] <- nids
  }
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_along(nb))
  
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  x <- m[!is.na(m)]
  mt <- spdep::moran.test(x, lw, zero.policy = TRUE)
  list(I = unname(mt$estimate[["Moran I statistic"]]),
       p = mt$p.value)
}

# Empirical variogram + simple range/sill extraction using a fitted model
variogram_metrics <- function(r, cutoff = 15, width = 1) {
  # Convert raster to points
  pts <- as.data.frame(terra::as.points(r), na.rm = TRUE)
  if (nrow(pts) < 30) return(list(range = NA_real_, sill = NA_real_))
  names(pts)[1:2] <- c("x", "y")
  names(pts)[3]   <- "v"
  sfpts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = NA)
  sppts <- as(sfpts, "Spatial")
  
  vg <- gstat::variogram(v ~ 1, sppts, cutoff = cutoff, width = width)
  # Fit a common model (exponential). Fall back if fit fails.
  fit <- try(gstat::fit.variogram(vg, model = gstat::vgm(psill = var(pts$v), model = "Exp", range = cutoff/3, nugget = 0)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(list(range = NA_real_, sill = NA_real_))
  # "range" is model range parameter; sill = nugget + psill
  rng <- fit$range[fit$model != "Nug"]
  nug <- fit$psill[fit$model == "Nug"]
  ps  <- fit$psill[fit$model != "Nug"]
  list(range = ifelse(length(rng) > 0, rng[1], NA_real_),
       sill  = ifelse(length(ps) > 0, ps[1] + ifelse(length(nug) > 0, nug[1], 0), NA_real_))
}

# Lacunarity (gliding box) for a binary matrix (0/1) across box sizes
lacunarity_gliding <- function(binmat, box_sizes) {
  # binmat: matrix with 0/1, NA allowed (treated as NA)
  out <- data.frame(box = box_sizes, lacunarity = NA_real_)
  nr <- nrow(binmat); nc <- ncol(binmat)
  
  for (b in box_sizes) {
    if (b > nr || b > nc) next
    sums <- c()
    for (i in 1:(nr - b + 1)) {
      for (j in 1:(nc - b + 1)) {
        w <- binmat[i:(i+b-1), j:(j+b-1)]
        if (all(is.na(w))) next
        # treat NA as 0 (or skip) — here: skip windows with too many NA
        if (mean(is.na(w)) > 0.5) next
        w[is.na(w)] <- 0
        sums <- c(sums, sum(w))
      }
    }
    if (length(sums) < 10) next
    mu <- mean(sums)
    if (mu == 0) {
      out$lacunarity[out$box == b] <- NA_real_
    } else {
      out$lacunarity[out$box == b] <- (var(sums) / (mu^2)) + 1
    }
  }
  out
}

# ============================================================
# Load + normalize point cloud
# ============================================================

las <- lidR::readLAS(las_path)
stopifnot(!lidR::is.empty(las))

if (normalize) {
  # Classify ground and build DTM, then normalize heights
  las <- lidR::lasground(las, algorithm = lidR::csf(cloth_resolution = csf_cloth_resolution))
  dtm <- lidR::rasterize_terrain(las, res = dtm_res, algorithm = lidR::tin())
  las <- lidR::normalize_height(las, dtm)
  # Remove negative heights (often artefacts)
  las <- lidR::lasfilter(las, Z >= -0.05)
}

# ============================================================
# Filter regeneration band
# ============================================================

las_regen <- lidR::lasfilter(las, Z >= zmin_regen & Z <= zmax_regen)
if (lidR::is.empty(las_regen)) stop("No points found in the specified regeneration height band.")

# ============================================================
# Build a horizontal grid and compute per-cell structural proxies
# ============================================================

# Use terra raster template covering point cloud extent
ext <- terra::ext(las)
r_template <- terra::rast(extent = ext, resolution = grid_res, crs = NA)

# Per-cell: point count in regen band (raw density proxy)
# lidR grid_metrics returns a raster-like object; convert to terra
r_n <- lidR::grid_metrics(las_regen, ~length(Z), res = grid_res)
r_n <- terra::rast(r_n)

# Per-cell: height percentiles / IQR in regen band
r_h50 <- terra::rast(lidR::grid_metrics(las_regen, ~quantile(Z, 0.50, na.rm=TRUE), res = grid_res))
r_h95 <- terra::rast(lidR::grid_metrics(las_regen, ~quantile(Z, 0.95, na.rm=TRUE), res = grid_res))
r_iqr <- terra::rast(lidR::grid_metrics(las_regen, ~IQR(Z, na.rm=TRUE), res = grid_res))

# Per-cell: FHD of vertical profile in regen band (entropy)
r_fhd <- terra::rast(lidR::grid_metrics(
  las_regen,
  ~compute_fhd(Z, bin = zbin, zmin = zmin_regen, zmax = zmax_regen),
  res = grid_res
))

# ============================================================
# Voxel occupancy fraction (regen band) as a robust structure proxy
# ============================================================

# Create voxel indices
dt <- data.table::as.data.table(las_regen@data[, c("X","Y","Z")])
dt[, vx := floor(X / voxel_xy)]
dt[, vy := floor(Y / voxel_xy)]
dt[, vz := floor(Z / voxel_z)]

# Unique occupied voxels
vox <- unique(dt[, .(vx, vy, vz)])

# Map voxels to coarse grid cells (for horizontal mapping)
# convert voxel XY back to meters at voxel origin
vox[, X0 := vx * voxel_xy]
vox[, Y0 := vy * voxel_xy]

# grid cell indices
vox[, gx := floor(X0 / grid_res)]
vox[, gy := floor(Y0 / grid_res)]

# occupied voxel count per grid cell
occ_vox <- vox[, .N, by = .(gx, gy)]
setnames(occ_vox, "N", "n_occ_vox")

# total possible voxels in regen band per grid cell (approx)
n_z_bins <- ceiling((zmax_regen - zmin_regen) / voxel_z)
# number of voxels in one grid cell horizontally:
n_vx_per_cell <- ceiling(grid_res / voxel_xy)
n_vy_per_cell <- ceiling(grid_res / voxel_xy)
n_vox_possible <- n_vx_per_cell * n_vy_per_cell * n_z_bins

occ_vox[, occ_vox_frac := n_occ_vox / n_vox_possible]

# Put occ_vox_frac into a raster aligned with r_n
# Create raster cell centers and map gx/gy back to coordinates
# We'll build an xy table at grid cell origins
occ_vox[, x := (gx + 0.5) * grid_res]
occ_vox[, y := (gy + 0.5) * grid_res]
r_occ <- terra::rast(r_template)
# terra::cellFromXY needs coordinates within extent
cells <- terra::cellFromXY(r_occ, as.matrix(occ_vox[, .(x, y)]))
r_occ[cells] <- occ_vox$occ_vox_frac

# ============================================================
# Plot-level (whole cloud) height distribution + vertical heterogeneity
# ============================================================

z_all <- las_regen@data$Z
plot_metrics <- data.frame(
  zmin = zmin_regen,
  zmax = zmax_regen,
  n_points_regen = length(z_all),
  H10 = unname(quantile(z_all, 0.10, na.rm = TRUE)),
  H50 = unname(quantile(z_all, 0.50, na.rm = TRUE)),
  H90 = unname(quantile(z_all, 0.90, na.rm = TRUE)),
  H95 = unname(quantile(z_all, 0.95, na.rm = TRUE)),
  IQR = IQR(z_all, na.rm = TRUE),
  MAD = mad(z_all, constant = 1, na.rm = TRUE),
  FHD = compute_fhd(z_all, bin = zbin, zmin = zmin_regen, zmax = zmax_regen),
  FHD_effective_layers = exp(compute_fhd(z_all, bin = zbin, zmin = zmin_regen, zmax = zmax_regen))
)

# ============================================================
# Horizontal heterogeneity metrics
# (Moran's I, variogram range/sill, lacunarity)
# ============================================================

# Choose a primary "regen density" map:
# - occupied voxel fraction is robust to point density effects
primary_map <- r_occ

moran_primary <- moran_I_raster(primary_map, queen = TRUE)
vg_primary    <- variogram_metrics(primary_map, cutoff = variogram_cutoff, width = variogram_width)

# Lacunarity on a binary "regen present" raster:
# threshold: present if occ_vox_frac > 0 (or a small threshold like 0.01)
r_bin <- primary_map > 0
binmat <- as.matrix(r_bin, wide = TRUE)
# convert TRUE/FALSE/NA to 1/0/NA
binmat <- ifelse(is.na(binmat), NA, ifelse(binmat, 1, 0))
lac <- lacunarity_gliding(binmat, lac_box_sizes)

# Additional simple inequality summaries across cells
vals_primary <- terra::values(primary_map, na.rm = TRUE)
horiz_summ <- data.frame(
  mean_primary = mean(vals_primary, na.rm = TRUE),
  sd_primary   = sd(vals_primary, na.rm = TRUE),
  cv_primary   = cv(vals_primary),
  gini_primary = gini(vals_primary),
  moran_I      = moran_primary$I,
  moran_p      = moran_primary$p,
  variogram_range = vg_primary$range,
  variogram_sill  = vg_primary$sill
)

# ============================================================
# "PAI-like" proxy from voxel transmittance (CAUTION)
# ============================================================
# This is NOT true PAI unless you have ray/pulse info.
# Proxy idea: for each XY column, compute fraction of empty voxels in regen band (gap proxy),
# then apply Beer–Lambert: ePAI_proxy = -ln(gap)
# where gap = empty / total (clipped to avoid log(0))

# Build per-(vx,vy) vertical occupancy presence
col_occ <- unique(vox[, .(vx, vy, vz)])
# total z bins in regen band for this voxelization
# we need to align z bins: shift by zmin_regen
dt_all <- data.table::as.data.table(las@data[, c("X","Y","Z")])
dt_all <- dt_all[Z >= zmin_regen & Z <= zmax_regen]
dt_all[, vx := floor(X / voxel_xy)]
dt_all[, vy := floor(Y / voxel_xy)]
dt_all[, vz := floor((Z - zmin_regen) / voxel_z)]
col_occ <- unique(dt_all[, .(vx, vy, vz)])

# occupancy count per column
col_counts <- col_occ[, .N, by = .(vx, vy)]
setnames(col_counts, "N", "n_occ_z")
col_counts[, n_total_z := ceiling((zmax_regen - zmin_regen) / voxel_z)]
col_counts[, gap := pmax(1e-6, pmin(1 - 1e-6, (n_total_z - n_occ_z) / n_total_z))]
col_counts[, ePAI_proxy := -log(gap)]

# Aggregate to grid cells (horizontal map)
col_counts[, X0 := vx * voxel_xy]
col_counts[, Y0 := vy * voxel_xy]
col_counts[, gx := floor(X0 / grid_res)]
col_counts[, gy := floor(Y0 / grid_res)]
pai_cell <- col_counts[, .(ePAI_proxy = mean(ePAI_proxy, na.rm = TRUE)), by = .(gx, gy)]
pai_cell[, x := (gx + 0.5) * grid_res]
pai_cell[, y := (gy + 0.5) * grid_res]
r_pai <- terra::rast(r_template)
cells2 <- terra::cellFromXY(r_pai, as.matrix(pai_cell[, .(x, y)]))
r_pai[cells2] <- pai_cell$ePAI_proxy

# Spatial heterogeneity of leafiness proxy
moran_pai <- moran_I_raster(r_pai, queen = TRUE)
vg_pai    <- variogram_metrics(r_pai, cutoff = variogram_cutoff, width = variogram_width)
vals_pai  <- terra::values(r_pai, na.rm = TRUE)

pai_summ <- data.frame(
  mean_ePAI_proxy = mean(vals_pai, na.rm = TRUE),
  cv_ePAI_proxy   = cv(vals_pai),
  gini_ePAI_proxy = gini(vals_pai),
  moran_I_ePAI_proxy = moran_pai$I,
  moran_p_ePAI_proxy = moran_pai$p,
  variogram_range_ePAI_proxy = vg_pai$range,
  variogram_sill_ePAI_proxy  = vg_pai$sill
)

# ============================================================
# Save outputs
# ============================================================

dir.create("tls_metrics_out", showWarnings = FALSE)

write.csv(plot_metrics, "tls_metrics_out/plot_level_height_metrics.csv", row.names = FALSE)
write.csv(horiz_summ,   "tls_metrics_out/horizontal_heterogeneity_primary.csv", row.names = FALSE)
write.csv(lac,          "tls_metrics_out/lacunarity_primary_binary.csv", row.names = FALSE)
write.csv(pai_summ,     "tls_metrics_out/leafiness_ePAI_proxy_summary.csv", row.names = FALSE)

# Save rasters for mapping/QA
terra::writeRaster(r_n,    "tls_metrics_out/regen_pointcount.tif", overwrite = TRUE)
terra::writeRaster(r_occ,  "tls_metrics_out/regen_occ_vox_frac.tif", overwrite = TRUE)
terra::writeRaster(r_h50,  "tls_metrics_out/regen_H50.tif", overwrite = TRUE)
terra::writeRaster(r_h95,  "tls_metrics_out/regen_H95.tif", overwrite = TRUE)
terra::writeRaster(r_iqr,  "tls_metrics_out/regen_IQR.tif", overwrite = TRUE)
terra::writeRaster(r_fhd,  "tls_metrics_out/regen_FHD.tif", overwrite = TRUE)
terra::writeRaster(r_pai,  "tls_metrics_out/ePAI_proxy.tif", overwrite = TRUE)

cat("Done. Outputs written to tls_metrics_out/\n")
