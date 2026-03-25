# Calgary fishnet and target generation
# This script is intentionally self-contained so it can be sourced from the Quarto report.

build_fishnet_and_target <- function(config = list()) {
  if (!is.list(config)) {
    stop("config must be a list.", call. = FALSE)
  }

  required_pkgs <- c("sf", "terra", "exactextractr", "dplyr", "readr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them before running build_fishnet_and_target().",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(sf)
    library(terra)
    library(exactextractr)
    library(dplyr)
    library(readr)
  })

  `%||%` <- function(x, y) if (!is.null(x)) x else y
  paths <- config$paths %||% list()
  project_root <- config$project_root %||% getwd()

  # Keep shapefile reads stable even if the .shx index is missing.
  Sys.setenv(SHAPE_RESTORE_SHX = "YES")

  base_dir <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  data_dir <- normalizePath(
    paths$data_dir %||% file.path(base_dir, "data"),
    winslash = "/",
    mustWork = TRUE
  )
  intermediate_dir <- paths$intermediate_dir %||% file.path(base_dir, "outputs", "intermediate")
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

  base_fishnet_path <- paths$fishnet_base %||% file.path(intermediate_dir, "fishnet_30m.gpkg")
  target_fishnet_path <- paths$fishnet_target %||% file.path(intermediate_dir, "fishnet_target.gpkg")
  target_summary_path <- paths$target_summary %||% file.path(intermediate_dir, "fishnet_target_summary.csv")

  boundary_path <- file.path(data_dir, "calgary_boundary.shp")
  inundation_candidates <- c(
    file.path(data_dir, "inundation_calgary.tif"),
    file.path(data_dir, "inundation.tif")
  )
  inundation_path <- inundation_candidates[file.exists(inundation_candidates)][1]
  if (is.na(inundation_path) || !nzchar(inundation_path)) {
    stop("Could not find a Calgary inundation raster in data_dir.", call. = FALSE)
  }
  if (!file.exists(boundary_path)) {
    stop("Missing boundary shapefile: ", boundary_path, call. = FALSE)
  }

  fishnet_size_m <- as.numeric(config$fishnet_size_m %||% 30)
  classification_threshold <- as.numeric(config$classification_threshold %||% 0.30)
  target_col <- config$target_col %||% "target"
  if (is.na(fishnet_size_m) || fishnet_size_m <= 0) {
    stop("config$fishnet_size_m must be a positive number.", call. = FALSE)
  }
  if (fishnet_size_m %% 10 != 0) {
    stop("fishnet_size_m must be divisible by 10 for the inundation target workflow.", call. = FALSE)
  }
  agg_fact <- as.integer(fishnet_size_m / 10)

  # Read and repair the Calgary boundary, then use its CRS as the analysis CRS.
  boundary_sf <- sf::st_read(boundary_path, quiet = TRUE)
  boundary_sf <- sf::st_make_valid(boundary_sf)
  boundary_sf <- boundary_sf |>
    dplyr::filter(!sf::st_is_empty(geometry))
  analysis_crs <- sf::st_crs(boundary_sf)
  if (nrow(boundary_sf) < 1) {
    stop("Boundary layer has no valid geometries after repair.", call. = FALSE)
  }

  if (!is.null(config$analysis_crs_input) && nzchar(as.character(config$analysis_crs_input))) {
    message("Using boundary CRS for analysis. Ignoring config$analysis_crs_input in favor of the boundary CRS.")
  }

  boundary_union <- sf::st_union(boundary_sf)
  boundary_union_sf <- sf::st_sf(geometry = sf::st_sfc(boundary_union, crs = analysis_crs))
  boundary_vect <- terra::vect(boundary_union_sf)

  # Read the inundation raster and align it to the boundary CRS if needed.
  inundation_rast <- terra::rast(inundation_path)
  if (!terra::same.crs(inundation_rast, boundary_vect)) {
    inundation_rast <- terra::project(inundation_rast, boundary_vect, method = "near")
  }

  # Reclassify the raster so only 0/1 are valid and all other codes become NA.
  inundation_clean <- terra::ifel(
    inundation_rast == 0 | inundation_rast == 1,
    inundation_rast,
    NA
  )
  valid_mask_10 <- terra::ifel(is.na(inundation_clean), NA, 1)
  cell_id_10 <- terra::init(inundation_clean, "cell")

  # Aggregate the 10 m raster to the requested fishnet size so interior cells can be summarized quickly.
  valid_count_30 <- terra::aggregate(valid_mask_10, fact = agg_fact, fun = "sum", na.rm = TRUE)
  inund_count_30 <- terra::aggregate(inundation_clean, fact = agg_fact, fun = "sum", na.rm = TRUE)
  cell_id_30 <- terra::aggregate(cell_id_10, fact = agg_fact, fun = "min", na.rm = TRUE)
  names(cell_id_30) <- "id"
  names(valid_count_30) <- "valid_inund_cells"
  names(inund_count_30) <- "inundated_cells"

  # Convert the aggregated raster blocks to polygons and clip them to the Calgary boundary.
  fishnet_vect <- terra::as.polygons(cell_id_30, dissolve = FALSE, values = TRUE, na.rm = TRUE)
  fishnet_vect <- terra::intersect(fishnet_vect, boundary_vect)
  fishnet_sf <- sf::st_as_sf(fishnet_vect)
  fishnet_sf <- sf::st_make_valid(fishnet_sf)
  fishnet_sf <- fishnet_sf |>
    dplyr::filter(!sf::st_is_empty(geometry))
  fishnet_sf$id <- as.integer(fishnet_sf$id)

  if (anyDuplicated(fishnet_sf$id) > 0) {
    stop("Fishnet id values are not unique after clipping.", call. = FALSE)
  }

  # Build a fast target table from the raster aggregates.
  target_df <- data.frame(
    id = as.integer(values(cell_id_30, mat = FALSE)),
    valid_inund_cells = as.numeric(values(valid_count_30, mat = FALSE)),
    inundated_cells = as.numeric(values(inund_count_30, mat = FALSE))
  )
  target_df$inundated_share <- dplyr::if_else(target_df$valid_inund_cells > 0, target_df$inundated_cells / target_df$valid_inund_cells, NA_real_)
  target_df[[target_col]] <- ifelse(!is.na(target_df$inundated_share) & target_df$inundated_share >= classification_threshold, 1L, 0L)
  target_df[[target_col]] <- as.integer(target_df[[target_col]])

  # Correct edge cells by computing exact polygon summaries only where clipping changed the cell area.
  full_cell_area <- fishnet_size_m^2
  fishnet_sf$cell_area_m2 <- as.numeric(sf::st_area(fishnet_sf))
  edge_sf <- fishnet_sf |>
    dplyr::filter(abs(cell_area_m2 - full_cell_area) > 1e-06)

  if (nrow(edge_sf) > 0) {
    summarize_target <- function(df) {
      valid <- !is.na(df$value) & df$value %in% c(0, 1)
      valid_inund_cells <- sum(df$coverage_fraction[valid], na.rm = TRUE)
      inundated_cells <- sum(df$coverage_fraction[df$value == 1], na.rm = TRUE)
      inundated_share <- if (valid_inund_cells > 0) inundated_cells / valid_inund_cells else NA_real_
      c(valid_inund_cells = valid_inund_cells, inundated_cells = inundated_cells, inundated_share = inundated_share)
    }

    edge_matrix <- exactextractr::exact_extract(
      inundation_clean,
      edge_sf,
      summarize_target,
      summarize_df = TRUE
    )
    edge_df <- as.data.frame(t(edge_matrix), check.names = FALSE, stringsAsFactors = FALSE)
    edge_df$id <- edge_sf$id

    target_df <- target_df |>
      dplyr::left_join(edge_df |> dplyr::select(id, valid_inund_cells, inundated_cells, inundated_share), by = "id", suffix = c("_agg", "_exact")) |>
      dplyr::mutate(
        valid_inund_cells = dplyr::coalesce(valid_inund_cells_exact, valid_inund_cells_agg),
        inundated_cells = dplyr::coalesce(inundated_cells_exact, inundated_cells_agg),
        inundated_share = dplyr::coalesce(inundated_share_exact, inundated_share_agg)
      ) |>
      dplyr::select(id, valid_inund_cells, inundated_cells, inundated_share, all_of(target_col))
  }

  # Remove cells with no valid inundation coverage so downstream modeling is clean.
  fishnet_target_sf <- fishnet_sf |>
    dplyr::left_join(target_df, by = "id") |>
    dplyr::filter(!is.na(valid_inund_cells) & valid_inund_cells > 0)
  if (nrow(fishnet_target_sf) == 0) {
    stop("No fishnet cells retained after applying the valid inundation coverage filter.", call. = FALSE)
  }

  fishnet_target_sf$valid_inund_cells <- as.numeric(fishnet_target_sf$valid_inund_cells)
  fishnet_target_sf$inundated_cells <- as.numeric(fishnet_target_sf$inundated_cells)
  fishnet_target_sf$inundated_share <- as.numeric(fishnet_target_sf$inundated_share)
  fishnet_target_sf[[target_col]] <- as.integer(fishnet_target_sf[[target_col]])

  # Drop helper geometry area from the final target fishnet.
  fishnet_target_sf$cell_area_m2 <- NULL
  fishnet_sf$cell_area_m2 <- NULL

  # Write geospatial outputs for later feature engineering and modeling.
  if (file.exists(base_fishnet_path)) unlink(base_fishnet_path)
  if (file.exists(target_fishnet_path)) unlink(target_fishnet_path)
  sf::st_write(fishnet_sf, base_fishnet_path, delete_layer = TRUE, quiet = TRUE, layer_options = "SPATIAL_INDEX=NO")
  sf::st_write(fishnet_target_sf, target_fishnet_path, delete_layer = TRUE, quiet = TRUE, layer_options = "SPATIAL_INDEX=NO")

  # Save a compact QA summary for the report and later checks.
  target_summary <- fishnet_target_sf |>
    sf::st_drop_geometry() |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_target_1 = sum(.data[[target_col]] == 1L, na.rm = TRUE),
      n_target_0 = sum(.data[[target_col]] == 0L, na.rm = TRUE),
      share_target_1 = mean(.data[[target_col]] == 1L, na.rm = TRUE),
      min_valid_inund_cells = min(valid_inund_cells, na.rm = TRUE),
      max_valid_inund_cells = max(valid_inund_cells, na.rm = TRUE)
    )
  readr::write_csv(target_summary, target_summary_path)

  # Print reproducibility checks after the major step.
  message("Boundary CRS: ", analysis_crs$input)
  message("Fishnet size (m): ", fishnet_size_m)
  message("Base fishnet rows: ", nrow(fishnet_sf))
  message("Target-valid fishnet rows: ", nrow(fishnet_target_sf))
  message("Duplicate fishnet ids: ", anyDuplicated(fishnet_sf$id))
  message("Duplicate target fishnet ids: ", anyDuplicated(fishnet_target_sf$id))
  message("Missing target rows retained: ", sum(is.na(fishnet_target_sf[[target_col]])))
  message("Target class counts:")
  print(table(fishnet_target_sf[[target_col]], useNA = "ifany"))
  message("Output written to: ", base_fishnet_path)
  message("Output written to: ", target_fishnet_path)
  message("Target summary written to: ", target_summary_path)

  invisible(list(
    analysis_crs = analysis_crs,
    fishnet_base = fishnet_sf,
    fishnet_target = fishnet_target_sf,
    paths = list(
      base_fishnet_path = base_fishnet_path,
      target_fishnet_path = target_fishnet_path,
      target_summary_path = target_summary_path
    )
  ))
}




