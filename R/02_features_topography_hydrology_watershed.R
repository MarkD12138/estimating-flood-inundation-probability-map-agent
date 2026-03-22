## Calgary feature engineering worker B
## Exports one function:
##   compute_topography_hydrology_watershed(config)
## The function reads the target-valid fishnet plus canonical Calgary inputs,
## computes topography, hydrology, and watershed features, then writes the
## requested CSV outputs and hydrology cache artifacts under outputs/intermediate.

`%||%` <- function(x, y) if (!is.null(x)) x else y

resolve_path <- function(path, project_root) {
  if (is.null(path) || identical(path, "")) return(NULL)
  is_abs <- grepl("^[A-Za-z]:[\\/]", path) || startsWith(path, "\\\\")
  if (is_abs) normalizePath(path, winslash = "\\", mustWork = FALSE) else normalizePath(file.path(project_root, path), winslash = "\\", mustWork = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_csv_checked <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_csv(x, path, na = "")
  message("Wrote CSV: ", path)
  invisible(path)
}

safe_st_read <- function(path, quiet = TRUE) {
  old <- Sys.getenv("SHAPE_RESTORE_SHX", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("SHAPE_RESTORE_SHX")
    } else {
      Sys.setenv(SHAPE_RESTORE_SHX = old)
    }
  }, add = TRUE)
  Sys.setenv(SHAPE_RESTORE_SHX = "YES")
  sf::st_read(path, quiet = quiet)
}

drop_empty_geoms <- function(x) {
  x <- sf::st_make_valid(x)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  x
}

weighted_mean <- function(values, coverage_fraction) {
  ok <- !is.na(values) & !is.na(coverage_fraction) & coverage_fraction > 0
  if (!any(ok)) return(NA_real_)
  values <- values[ok]
  coverage_fraction <- coverage_fraction[ok]
  wsum <- sum(coverage_fraction)
  if (!is.finite(wsum) || wsum <= 0) return(NA_real_)
  sum(values * coverage_fraction) / wsum
}

weighted_sd <- function(values, coverage_fraction) {
  ok <- !is.na(values) & !is.na(coverage_fraction) & coverage_fraction > 0
  if (!any(ok)) return(NA_real_)
  values <- values[ok]
  coverage_fraction <- coverage_fraction[ok]
  wsum <- sum(coverage_fraction)
  if (!is.finite(wsum) || wsum <= 0) return(NA_real_)
  mu <- sum(values * coverage_fraction) / wsum
  sqrt(sum(coverage_fraction * (values - mu)^2) / wsum)
}

weighted_stat <- function(rast, polys, stat, fallback_name = stat) {
  exact_fun <- switch(
    stat,
    mean = function(values, coverage_fraction) weighted_mean(values, coverage_fraction),
    min = function(values, coverage_fraction) {
      ok <- !is.na(values)
      if (!any(ok)) return(NA_real_)
      min(values[ok])
    },
    max = function(values, coverage_fraction) {
      ok <- !is.na(values)
      if (!any(ok)) return(NA_real_)
      max(values[ok])
    },
    sd = function(values, coverage_fraction) weighted_sd(values, coverage_fraction),
    stop("Unsupported stat: ", stat)
  )

  if (requireNamespace("exactextractr", quietly = TRUE)) {
    out <- tryCatch(
      exactextractr::exact_extract(rast, polys, exact_fun),
      error = function(e) {
        message("exactextractr failed for ", fallback_name, "; using terra::extract fallback. Reason: ", e$message)
        NULL
      }
    )
    if (!is.null(out)) return(out)
  }

  message("Using terra::extract fallback for ", fallback_name)
  tmp <- terra::extract(rast, terra::vect(polys))
  if (ncol(tmp) < 2) return(rep(NA_real_, nrow(polys)))
  values <- tmp[[2]]
  ids <- tmp[[1]]
  out <- tapply(values, ids, function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_real_)
    switch(
      stat,
      mean = mean(v),
      min = min(v),
      max = max(v),
      sd = stats::sd(v)
    )
  })
  as.numeric(out[match(seq_len(nrow(polys)), as.integer(names(out)))])
}

water_area_stat <- function(rast, polys, cell_area_m2) {
  if (requireNamespace("exactextractr", quietly = TRUE)) {
    out <- tryCatch(
      exactextractr::exact_extract(
        rast,
        polys,
        function(values, coverage_fraction) {
          ok <- !is.na(values) & values == 1
          sum(coverage_fraction[ok], na.rm = TRUE)
        }
      ),
      error = function(e) {
        message("exactextractr failed for water_cover_area; using terra::extract fallback. Reason: ", e$message)
        NULL
      }
    )
    if (!is.null(out)) return(out * cell_area_m2)
  }

  message("Using terra::extract fallback for water_cover_area")
  tmp <- terra::extract(rast, terra::vect(polys))
  if (ncol(tmp) < 2) return(rep(NA_real_, nrow(polys)))
  values <- tmp[[2]]
  ids <- tmp[[1]]
  out <- tapply(values, ids, function(v) sum(v == 1, na.rm = TRUE) * cell_area_m2)
  as.numeric(out[match(seq_len(nrow(polys)), as.integer(names(out)))])
}

sum_intersection_lengths <- function(fishnet, waterways) {
  inter <- suppressWarnings(sf::st_intersection(fishnet["id"], waterways))
  if (nrow(inter) == 0) {
    return(dplyr::tibble(id = fishnet$id, stream_len_m = 0))
  }
  inter$stream_len_m <- as.numeric(sf::st_length(inter))
  inter |>
    sf::st_drop_geometry() |>
    dplyr::group_by(id) |>
    dplyr::summarise(stream_len_m = sum(stream_len_m, na.rm = TRUE), .groups = "drop")
}

run_whitebox_watershed <- function(hyd_dem_path, centroids, whitebox_dir, snap_dist_m = 100, stream_threshold = 250) {
  if (!requireNamespace("whitebox", quietly = TRUE)) {
    return(list(available = FALSE, reason = "whitebox package is not installed."))
  }
  if (!isTRUE(whitebox::check_whitebox_binary(silent = TRUE))) {
    return(list(available = FALSE, reason = "WhiteboxTools binary is not available."))
  }

  ensure_dir(whitebox_dir)
  filled_dem <- file.path(whitebox_dir, "hydro_dem_filled.tif")
  d8_pointer <- file.path(whitebox_dir, "hydro_d8_pointer.tif")
  flow_accum <- file.path(whitebox_dir, "hydro_flow_accum_cells.tif")
  stream_raster <- file.path(whitebox_dir, "hydro_streams.tif")
  centroid_pts <- file.path(whitebox_dir, "hydro_pour_points_centroids.shp")
  snapped_pts <- file.path(whitebox_dir, "hydro_pour_points_snapped.shp")
  pour_pts_raster <- file.path(whitebox_dir, "hydro_pour_points_raster.tif")
  watershed_raster <- file.path(whitebox_dir, "hydro_watershed_ids.tif")

  message("Running Whitebox watershed workflow in: ", whitebox_dir)
  whitebox::wbt_fill_depressions(dem = hyd_dem_path, output = filled_dem, fix_flats = TRUE)
  whitebox::wbt_d8_pointer(dem = filled_dem, output = d8_pointer)
  whitebox::wbt_d8_flow_accumulation(input = filled_dem, output = flow_accum, out_type = "cells")
  whitebox::wbt_extract_streams(flow_accum = flow_accum, output = stream_raster, threshold = stream_threshold)
  sf::st_write(centroids, centroid_pts, delete_dsn = TRUE, quiet = TRUE)
  whitebox::wbt_snap_pour_points(pour_pts = centroid_pts, flow_accum = flow_accum, output = snapped_pts, snap_dist = snap_dist_m)
  whitebox::wbt_vector_points_to_raster(input = snapped_pts, output = pour_pts_raster, field = "id", assign = "last", nodata = TRUE, base = flow_accum)
  whitebox::wbt_watershed(d8_pntr = d8_pointer, pour_pts = pour_pts_raster, output = watershed_raster)

  list(
    available = TRUE,
    filled_dem = filled_dem,
    d8_pointer = d8_pointer,
    flow_accum = flow_accum,
    stream_raster = stream_raster,
    snapped_pts = snapped_pts,
    pour_pts_raster = pour_pts_raster,
    watershed_raster = watershed_raster
  )
}

run_terra_watershed_fallback <- function(hyd_dem, centroids, whitebox_dir) {
  ensure_dir(whitebox_dir)
  dem_dir <- terra::direction(hyd_dem, from = FALSE, degrees = FALSE)
  flow_accum <- terra::flowAccumulation(dem_dir)
  watershed_raster <- terra::watershed(dem_dir, terra::vect(centroids))

  flow_accum_path <- file.path(whitebox_dir, "terra_flow_accum.tif")
  d8_path <- file.path(whitebox_dir, "terra_direction.tif")
  watershed_path <- file.path(whitebox_dir, "terra_watershed.tif")
  terra::writeRaster(dem_dir, d8_path, overwrite = TRUE)
  terra::writeRaster(flow_accum, flow_accum_path, overwrite = TRUE)
  terra::writeRaster(watershed_raster, watershed_path, overwrite = TRUE)

  list(
    available = TRUE,
    filled_dem = NULL,
    d8_pointer = d8_path,
    flow_accum = flow_accum_path,
    stream_raster = NULL,
    snapped_pts = NULL,
    pour_pts_raster = NULL,
    watershed_raster = watershed_path
  )
}

compute_topography_hydrology_watershed <- function(config) {
  if (!requireNamespace("sf", quietly = TRUE) ||
      !requireNamespace("terra", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("readr", quietly = TRUE)) {
    stop("Required packages are missing from the local library.", call. = FALSE)
  }

  sf::sf_use_s2(FALSE)
  config <- config %||% list()
  paths <- config$paths %||% list()

  project_root <- resolve_path(config$project_root %||% getwd(), getwd())
  data_dir <- resolve_path(paths$data_dir %||% file.path(project_root, "data"), project_root)
  intermediate_dir <- resolve_path(paths$intermediate_dir %||% file.path(project_root, "outputs", "intermediate"), project_root)
  whitebox_dir <- resolve_path(paths$whitebox_dir %||% file.path(intermediate_dir, "whitebox_cache"), project_root)

  target_valid_fishnet_path <- resolve_path(paths$fishnet_target %||% paths$fishnet_path %||% file.path(intermediate_dir, "fishnet_target.gpkg"), project_root)
  features_topography_path <- resolve_path(paths$features_topography %||% file.path(intermediate_dir, "features_topography.csv"), project_root)
  features_hydrology_path <- resolve_path(paths$features_hydrology %||% file.path(intermediate_dir, "features_hydrology.csv"), project_root)
  features_watershed_path <- resolve_path(paths$features_watershed %||% file.path(intermediate_dir, "features_watershed.csv"), project_root)

  boundary_path <- file.path(data_dir, "calgary_boundary.shp")
  waterways_path <- file.path(data_dir, "calgary_waterways.shp")
  dem_path <- file.path(data_dir, "calgary_dem.tif")
  surf_water_path <- file.path(data_dir, "calgary_surf_water.tif")

  ensure_dir(intermediate_dir)
  ensure_dir(whitebox_dir)

  if (!file.exists(target_valid_fishnet_path)) stop("Target-valid fishnet not found: ", target_valid_fishnet_path, call. = FALSE)
  if (!file.exists(boundary_path)) stop("Missing boundary shapefile: ", boundary_path, call. = FALSE)
  if (!file.exists(waterways_path)) stop("Missing waterways shapefile: ", waterways_path, call. = FALSE)
  if (!file.exists(dem_path)) stop("Missing DEM: ", dem_path, call. = FALSE)
  if (!file.exists(surf_water_path)) stop("Missing surface water raster: ", surf_water_path, call. = FALSE)

  message("Reading Calgary boundary, target-valid fishnet, and waterways...")
  boundary_sf <- drop_empty_geoms(safe_st_read(boundary_path, quiet = TRUE))
  fishnet_sf <- drop_empty_geoms(safe_st_read(target_valid_fishnet_path, quiet = TRUE))
  waterways_sf <- drop_empty_geoms(safe_st_read(waterways_path, quiet = TRUE))

  if (!"id" %in% names(fishnet_sf)) stop("Fishnet must contain an `id` column.", call. = FALSE)
  if (anyDuplicated(fishnet_sf$id) > 0) stop("Fishnet `id` values are not unique.", call. = FALSE)

  analysis_crs <- if (!is.null(config$analysis_crs_input)) sf::st_crs(config$analysis_crs_input) else sf::st_crs(boundary_sf)
  if (is.na(analysis_crs)) stop("Analysis CRS is missing.", call. = FALSE)
  message("Analysis CRS: ", analysis_crs$input)

  fishnet_sf <- sf::st_transform(fishnet_sf, analysis_crs)
  waterways_sf <- sf::st_transform(waterways_sf, analysis_crs)
  boundary_sf <- sf::st_transform(boundary_sf, analysis_crs)

  fishnet_area_m2 <- as.numeric(sf::st_area(fishnet_sf))
  if (any(!is.finite(fishnet_area_m2) | fishnet_area_m2 <= 0)) stop("Invalid fishnet cell areas detected.", call. = FALSE)

  message("Reading DEM and surface-water raster...")
  dem <- terra::rast(dem_path)
  surf_water <- terra::rast(surf_water_path)

  if (terra::crs(dem) != analysis_crs$wkt) {
    message("Reprojecting DEM to analysis CRS.")
    dem <- terra::project(dem, analysis_crs$wkt, method = "bilinear")
  }
  if (terra::crs(surf_water) != analysis_crs$wkt) {
    message("Reprojecting surface-water raster to analysis CRS.")
    surf_water <- terra::project(surf_water, analysis_crs$wkt, method = "near")
  }

  dem_slope <- terra::terrain(dem, v = "slope", unit = "degrees", neighbors = 8)
  surface_vals <- terra::freq(surf_water, digits = 0)
  message("Surface-water raster values:")
  print(surface_vals)

  message("Computing topography features...")
  topography_tbl <- dplyr::tibble(
    id = fishnet_sf$id,
    mean_elev = weighted_stat(dem, fishnet_sf, "mean", "mean_elev"),
    min_elev = weighted_stat(dem, fishnet_sf, "min", "min_elev"),
    max_elev = weighted_stat(dem, fishnet_sf, "max", "max_elev"),
    sd_elev = weighted_stat(dem, fishnet_sf, "sd", "sd_elev"),
    mean_slope = weighted_stat(dem_slope, fishnet_sf, "mean", "mean_slope"),
    max_slope = weighted_stat(dem_slope, fishnet_sf, "max", "max_slope")
  ) |>
    dplyr::mutate(elev_range = max_elev - min_elev) |>
    dplyr::select(id, mean_elev, min_elev, elev_range, mean_slope, max_slope, sd_elev)
  message("Topography rows: ", nrow(topography_tbl), " | missing values: ", sum(is.na(topography_tbl)))
  message("Topography summary:")
  print(summary(dplyr::select(topography_tbl, -id)))
  write_csv_checked(topography_tbl, features_topography_path)

  message("Computing hydrology features...")
  fishnet_centroids <- sf::st_centroid(fishnet_sf)
  fishnet_centroids$id <- fishnet_sf$id
  nearest_idx <- sf::st_nearest_feature(fishnet_centroids, waterways_sf)
  dist_nearest_stream <- as.numeric(sf::st_distance(fishnet_centroids, waterways_sf[nearest_idx, ], by_element = TRUE))
  stream_len_tbl <- sum_intersection_lengths(fishnet_sf, waterways_sf)
  water_cover_area <- water_area_stat(surf_water, fishnet_sf, cell_area_m2 = prod(terra::res(surf_water)))

  hydrology_tbl <- dplyr::tibble(
    id = fishnet_sf$id,
    dist_nearest_stream = dist_nearest_stream,
    water_cover_area = water_cover_area
  ) |>
    dplyr::left_join(stream_len_tbl, by = "id") |>
    dplyr::mutate(
      stream_len_m = dplyr::coalesce(stream_len_m, 0),
      river_density = stream_len_m / fishnet_area_m2
    ) |>
    dplyr::select(id, dist_nearest_stream, water_cover_area, river_density)
  message("Hydrology rows: ", nrow(hydrology_tbl), " | missing values: ", sum(is.na(hydrology_tbl)))
  message("Hydrology summary:")
  print(summary(dplyr::select(hydrology_tbl, -id)))
  write_csv_checked(hydrology_tbl, features_hydrology_path)

  message("Preparing square-grid DEM for watershed calculations...")
  hydro_res_m <- config$hydrology_resolution_m %||% 10
  hydro_template <- terra::rast(ext = terra::ext(dem), resolution = hydro_res_m, crs = analysis_crs$wkt)
  dem_square <- terra::resample(dem, hydro_template, method = "bilinear")
  dem_square_path <- file.path(whitebox_dir, "hydrology_dem_square_10m.tif")
  terra::writeRaster(dem_square, dem_square_path, overwrite = TRUE)
  message("Saved square hydrology DEM: ", dem_square_path)

  snap_dist_m <- config$snap_distance_m %||% 100
  stream_threshold <- config$stream_threshold_cells %||% 250
  watershed_backend <- NULL

  if (isTRUE(config$use_whitebox %||% TRUE)) {
    watershed_backend <- tryCatch(
      run_whitebox_watershed(
        hyd_dem_path = dem_square_path,
        centroids = fishnet_centroids,
        whitebox_dir = whitebox_dir,
        snap_dist_m = snap_dist_m,
        stream_threshold = stream_threshold
      ),
      error = function(e) {
        message("Whitebox watershed workflow failed: ", e$message)
        NULL
      }
    )
  }
  if (is.null(watershed_backend) || !isTRUE(watershed_backend$available)) {
    message("Falling back to terra-based watershed workflow.")
    watershed_backend <- run_terra_watershed_fallback(dem_square, fishnet_centroids, whitebox_dir)
  }

  watershed_rast <- terra::rast(watershed_backend$watershed_raster)
  watershed_cell_area_m2 <- prod(terra::res(watershed_rast))
  watershed_freq <- terra::freq(watershed_rast, digits = 0)
  watershed_tbl <- watershed_freq |>
    dplyr::as_tibble() |>
    dplyr::rename(id = value, cell_count = count) |>
    dplyr::filter(!is.na(id), id > 0) |>
    dplyr::mutate(
      id = as.integer(id),
      upstream_ws_area = cell_count * watershed_cell_area_m2
    ) |>
    dplyr::select(id, upstream_ws_area)

  if (nrow(watershed_tbl) == 0) {
    stop("Watershed computation returned no positive watershed ids.", call. = FALSE)
  }
  message("Watershed rows: ", nrow(watershed_tbl), " | missing values: ", sum(is.na(watershed_tbl)))
  message("Watershed summary:")
  print(summary(dplyr::select(watershed_tbl, -id)))
  write_csv_checked(watershed_tbl, features_watershed_path)

  message("Verifying one-row-per-id outputs...")
  stopifnot(anyDuplicated(topography_tbl$id) == 0)
  stopifnot(anyDuplicated(hydrology_tbl$id) == 0)
  stopifnot(anyDuplicated(watershed_tbl$id) == 0)
  stopifnot(nrow(topography_tbl) == nrow(fishnet_sf))
  stopifnot(nrow(hydrology_tbl) == nrow(fishnet_sf))

  message("Feature worker B completed successfully.")

  invisible(list(
    analysis_crs = analysis_crs,
    fishnet = fishnet_sf,
    topography = topography_tbl,
    hydrology = hydrology_tbl,
    watershed = watershed_tbl,
    caches = list(
      whitebox_dir = whitebox_dir,
      hydrology_dem_square = dem_square_path,
      watershed_backend = watershed_backend
    ),
    paths = list(
      features_topography = features_topography_path,
      features_hydrology = features_hydrology_path,
      features_watershed = features_watershed_path
    )
  ))
}
