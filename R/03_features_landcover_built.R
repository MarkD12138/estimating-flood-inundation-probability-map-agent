compute_landcover_built_features <- function(config = list()) {
  suppressPackageStartupMessages({
    library(sf)
    library(dplyr)
    library(readr)
    library(tibble)
  })

  # Allow shapefiles with missing .shx indexes to be restored on read.
  Sys.setenv(SHAPE_RESTORE_SHX = 'YES')

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  get_config_value <- function(name, default = NULL) {
    if (!is.null(config[[name]])) return(config[[name]])
    if (!is.null(config$paths) && !is.null(config$paths[[name]])) return(config$paths[[name]])
    default
  }

  resolve_path <- function(path_value, project_root = NULL) {
    if (is.null(path_value) || identical(path_value, "")) return(path_value)
    if (grepl("^(?:[A-Za-z]:[\\/]|/)", path_value) || is.null(project_root)) return(path_value)
    file.path(project_root, path_value)
  }

  clean_polygon_layer <- function(x, target_crs) {
    x <- st_as_sf(x)
    x <- x[!is.na(st_geometry(x)), , drop = FALSE]
    x <- x[!st_is_empty(x), , drop = FALSE]
    x <- st_make_valid(x)
    x <- x[!st_is_empty(x), , drop = FALSE]
    x <- suppressWarnings(st_collection_extract(x, "POLYGON", warn = FALSE))
    x <- x[!st_is_empty(x), , drop = FALSE]
    st_transform(x, target_crs)
  }

  summarize_overlap_ratio <- function(fishnet, layer, class_col, keep_values, out_name, area_lookup) {
    layer <- layer[layer[[class_col]] %in% keep_values, , drop = FALSE]
    if (nrow(layer) == 0L) {
      message(sprintf("No matched geometries for %s; returning zero ratios.", out_name))
      return(tibble(id = fishnet$id, value = 0))
    }

    overlap <- st_intersection(fishnet["id"], layer[c(class_col)])
    if (nrow(overlap) == 0L) {
      return(tibble(id = fishnet$id, value = 0))
    }

    overlap$area_m2 <- as.numeric(st_area(overlap))
    out <- overlap |>
      st_drop_geometry() |>
      group_by(id) |>
      summarise(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") |>
      mutate(value = area_m2 / area_lookup[match(id, fishnet$id)]) |>
      select(id, value)

    tibble(id = fishnet$id) |>
      left_join(out, by = "id") |>
      mutate(value = coalesce(value, 0)) |>
      rename(!!out_name := value)
  }

  summarize_building_ratio <- function(fishnet, buildings, out_name, area_lookup) {
    if (nrow(buildings) == 0L) {
      message("No matched building geometries; returning zero ratios.")
      return(tibble(id = fishnet$id, value = 0))
    }

    overlap <- st_intersection(fishnet["id"], buildings)
    if (nrow(overlap) == 0L) {
      return(tibble(id = fishnet$id, value = 0))
    }

    overlap$area_m2 <- as.numeric(st_area(overlap))
    out <- overlap |>
      st_drop_geometry() |>
      group_by(id) |>
      summarise(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop") |>
      mutate(value = area_m2 / area_lookup[match(id, fishnet$id)]) |>
      select(id, value)

    tibble(id = fishnet$id) |>
      left_join(out, by = "id") |>
      mutate(value = coalesce(value, 0)) |>
      rename(!!out_name := value)
  }

  project_root <- get_config_value("project_root", getwd())
  fishnet_path <- resolve_path(get_config_value("fishnet_target", get_config_value("fishnet_path", file.path("outputs", "intermediate", "fishnet_target.gpkg"))), project_root)
  landcover_path <- resolve_path(get_config_value("landcover_path", file.path("data", "calgary_landcover.geojson")), project_root)
  impervious_path <- resolve_path(get_config_value("impervious_path", file.path("data", "calgary_Imp_surf.geojson")), project_root)
  building_path <- resolve_path(get_config_value("building_path", file.path("data", "calgary_building.shp")), project_root)
  landcover_out <- resolve_path(get_config_value("features_landcover", get_config_value("landcover_out", file.path("outputs", "intermediate", "features_landcover.csv"))), project_root)
  built_out <- resolve_path(get_config_value("features_built", get_config_value("built_out", file.path("outputs", "intermediate", "features_built.csv"))), project_root)

  landcover_groups <- get_config_value("landcover_groups", list())
  vegetation_values <- landcover_groups$vegetation %||% c("Forest", "Grassland", "Shrubland", "Manicured", "GolfCourse", "Natural Wetland")
  open_soil_values <- landcover_groups$open_soil %||% c("Bare Ground")
  agriculture_values <- landcover_groups$agriculture %||% c("Agricultural_Annual Crops", "Agricultural_Pasture/Fallow")
  impervious_values <- get_config_value("impervious_values", c("Concrete", "Pavement"))

  if (!file.exists(fishnet_path)) stop("Fishnet file not found: ", fishnet_path)
  fishnet <- st_read(fishnet_path, quiet = TRUE, stringsAsFactors = FALSE)
  if (!"id" %in% names(fishnet)) stop("The fishnet layer must contain an 'id' field.")

  analysis_crs <- st_crs(fishnet)
  analysis_crs_input <- get_config_value("analysis_crs_input", NULL)
  if (!is.null(analysis_crs_input)) {
    analysis_crs_input <- st_crs(analysis_crs_input)
    if (!is.na(analysis_crs) && !is.na(analysis_crs_input) && analysis_crs$wkt != analysis_crs_input$wkt) {
      message("analysis_crs_input differs from the fishnet CRS; using the fishnet CRS as authoritative.")
    }
  }
  if (is.na(analysis_crs)) stop("Fishnet CRS is missing; cannot continue.")
  if (isTRUE(st_is_longlat(fishnet))) stop("Fishnet CRS must be projected in meters for area ratios.")

  fishnet <- fishnet[!is.na(st_geometry(fishnet)), , drop = FALSE]
  fishnet <- fishnet[!st_is_empty(fishnet), , drop = FALSE]
  fishnet <- st_make_valid(fishnet)
  fishnet <- fishnet[!st_is_empty(fishnet), , drop = FALSE]

  fishnet_area_by_id <- as.numeric(st_area(fishnet))
  if (any(!is.finite(fishnet_area_by_id))) stop("Fishnet area calculation returned non-finite values.")

  landcover <- clean_polygon_layer(st_read(landcover_path, quiet = TRUE, stringsAsFactors = FALSE), analysis_crs)
  impervious <- clean_polygon_layer(st_read(impervious_path, quiet = TRUE, stringsAsFactors = FALSE), analysis_crs)
  buildings <- clean_polygon_layer(st_read(building_path, quiet = TRUE, stringsAsFactors = FALSE), analysis_crs)

  if (!"lc_cat" %in% names(landcover)) stop("Land cover layer must contain an 'lc_cat' field.")
  if (!"gen_surface" %in% names(impervious)) stop("Impervious layer must contain a 'gen_surface' field.")

  landcover$lc_cat <- as.character(landcover$lc_cat)
  impervious$gen_surface <- as.character(impervious$gen_surface)

  vegetation_ratio <- summarize_overlap_ratio(fishnet, landcover, "lc_cat", vegetation_values, "vegetation_ratio", fishnet_area_by_id)
  open_soil_ratio <- summarize_overlap_ratio(fishnet, landcover, "lc_cat", open_soil_values, "open_soil_ratio", fishnet_area_by_id)
  agriculture_ratio <- summarize_overlap_ratio(fishnet, landcover, "lc_cat", agriculture_values, "agriculture_ratio", fishnet_area_by_id)
  impervious_ratio <- summarize_overlap_ratio(fishnet, impervious, "gen_surface", impervious_values, "impervious_ratio", fishnet_area_by_id)
  building_cover_ratio <- summarize_building_ratio(fishnet, buildings, "building_cover_ratio", fishnet_area_by_id)

  landcover_features <- vegetation_ratio |>
    left_join(open_soil_ratio, by = "id") |>
    left_join(agriculture_ratio, by = "id")

  built_features <- impervious_ratio |>
    left_join(building_cover_ratio, by = "id")

  fishnet_features <- fishnet |>
    left_join(landcover_features, by = "id") |>
    left_join(built_features, by = "id") |>
    mutate(
      impervious_ratio = coalesce(impervious_ratio, 0),
      vegetation_ratio = coalesce(vegetation_ratio, 0),
      open_soil_ratio = coalesce(open_soil_ratio, 0),
      agriculture_ratio = coalesce(agriculture_ratio, 0),
      building_cover_ratio = coalesce(building_cover_ratio, 0)
    )

  message("Land cover and built-environment feature computation complete.")
  message("Fishnet rows: ", nrow(fishnet_features))
  message("Duplicate id count: ", sum(duplicated(fishnet_features$id)))
  message("Missing values by feature:")
  print(colSums(is.na(st_drop_geometry(fishnet_features)[, c(
    "impervious_ratio",
    "vegetation_ratio",
    "open_soil_ratio",
    "agriculture_ratio",
    "building_cover_ratio"
  )])))
  message("Feature summaries:")
  print(summary(st_drop_geometry(fishnet_features)[, c(
    "impervious_ratio",
    "vegetation_ratio",
    "open_soil_ratio",
    "agriculture_ratio",
    "building_cover_ratio"
  )]))

  dir.create(dirname(landcover_out), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(built_out), recursive = TRUE, showWarnings = FALSE)

  landcover_csv <- landcover_features |>
    select(id, vegetation_ratio, open_soil_ratio, agriculture_ratio)
  built_csv <- built_features |>
    select(id, impervious_ratio, building_cover_ratio)

  write_csv(landcover_csv, landcover_out, na = "")
  write_csv(built_csv, built_out, na = "")

  invisible(list(
    fishnet = fishnet_features,
    landcover = landcover_csv,
    built = built_csv,
    paths = list(
      landcover_csv = landcover_out,
      built_csv = built_out
    )
  ))
}




