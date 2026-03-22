suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

resolve_path <- function(paths, primary, fallback = NULL) {
  if (!is.null(paths[[primary]]) && nzchar(paths[[primary]])) return(paths[[primary]])
  if (!is.null(fallback) && !is.null(paths[[fallback]]) && nzchar(paths[[fallback]])) return(paths[[fallback]])
  NULL
}

normalize_merge_eda_config <- function(config) {
  defaults <- list(
    project_root = ".",
    fishnet_size_m = 30,
    seed = 42,
    classification_threshold = 0.5,
    target_col = "target",
    predictor_cols = NULL,
    analysis_crs_input = NULL,
    landcover_groups = list(),
    impervious_values = c("Concrete", "Pavement"),
    paths = list(
      data_dir = "data",
      intermediate_dir = file.path("outputs", "intermediate"),
      diagnostics_dir = file.path("outputs", "diagnostics"),
      figures_dir = file.path("outputs", "figures"),
      eda_fig_dir = file.path("outputs", "figures", "eda"),
      diag_fig_dir = file.path("outputs", "figures", "diagnostics"),
      fishnet_base = file.path("outputs", "intermediate", "fishnet_30m.gpkg"),
      fishnet_target = file.path("outputs", "intermediate", "fishnet_target.gpkg"),
      features_topography = file.path("outputs", "intermediate", "features_topography.csv"),
      features_hydrology = file.path("outputs", "intermediate", "features_hydrology.csv"),
      features_watershed = file.path("outputs", "intermediate", "features_watershed.csv"),
      features_landcover = file.path("outputs", "intermediate", "features_landcover.csv"),
      features_built = file.path("outputs", "intermediate", "features_built.csv"),
      features_all = file.path("outputs", "intermediate", "features_all.gpkg"),
      modeling_raw = file.path("outputs", "intermediate", "modeling_table_raw.csv"),
      modeling_clean = file.path("outputs", "intermediate", "modeling_table_clean.csv"),
      model_rds = file.path("outputs", "models", "glm_logistic_calgary.rds"),
      scaling_params = file.path("outputs", "intermediate", "scaling_params.csv"),
      model_summary = file.path("outputs", "diagnostics", "model_summary.csv"),
      confusion_metrics = file.path("outputs", "diagnostics", "confusion_metrics.csv"),
      outcome_definitions = file.path("outputs", "diagnostics", "outcome_definitions.csv"),
      roc_data = file.path("outputs", "diagnostics", "roc_data.csv"),
      roc_plot = file.path("outputs", "figures", "diagnostics", "roc_curve.png"),
      predictions_full = file.path("outputs", "intermediate", "predictions_full.csv"),
      train_eval = file.path("outputs", "intermediate", "train_eval.csv"),
      test_eval = file.path("outputs", "intermediate", "test_eval.csv"),
      prediction_map = file.path("outputs", "figures", "diagnostics", "prediction_map.png"),
      train_outcome_map = file.path("outputs", "figures", "diagnostics", "train_outcome_map.png"),
      test_outcome_map = file.path("outputs", "figures", "diagnostics", "test_outcome_map.png"),
      missingness_csv = file.path("outputs", "diagnostics", "missingness_summary.csv"),
      corr_csv = file.path("outputs", "diagnostics", "correlation_matrix.csv"),
      variance_csv = file.path("outputs", "diagnostics", "variance_check.csv"),
      outlier_csv = file.path("outputs", "diagnostics", "outlier_summary.csv"),
      class_balance_csv = file.path("outputs", "diagnostics", "class_balance.csv"),
      whitebox_dir = file.path("outputs", "intermediate", "whitebox")
    ),
    eda = list(
      map_priority = c("target", "mean_elev", "upstream_ws_area", "impervious_ratio", "water_cover_area", "river_density"),
      histogram_exclude = c("id"),
      iqr_multiplier = 1.5
    )
  )

  if (is.null(config) || !is.list(config)) return(defaults)
  if (!is.null(config$paths)) defaults$paths[names(config$paths)] <- config$paths
  if (!is.null(config$eda)) defaults$eda[names(config$eda)] <- config$eda
  if (!is.null(config$project_root)) defaults$project_root <- config$project_root
  if (!is.null(config$fishnet_size_m)) defaults$fishnet_size_m <- config$fishnet_size_m
  if (!is.null(config$seed)) defaults$seed <- config$seed
  if (!is.null(config$classification_threshold)) defaults$classification_threshold <- config$classification_threshold
  if (!is.null(config$target_col)) defaults$target_col <- config$target_col
  if (!is.null(config$predictor_cols)) defaults$predictor_cols <- config$predictor_cols
  if (!is.null(config$analysis_crs_input)) defaults$analysis_crs_input <- config$analysis_crs_input
  if (!is.null(config$landcover_groups)) defaults$landcover_groups <- config$landcover_groups
  if (!is.null(config$impervious_values)) defaults$impervious_values <- config$impervious_values
  defaults
}

ensure_merge_eda_dirs <- function(paths) {
  dirs <- unique(c(paths$intermediate_dir, paths$diagnostics_dir, paths$figures_dir, paths$eda_fig_dir, paths$diag_fig_dir, dirname(paths$fishnet_base), dirname(paths$fishnet_target), dirname(paths$features_all), dirname(paths$modeling_raw), dirname(paths$modeling_clean), dirname(paths$model_rds), dirname(paths$scaling_params), dirname(paths$model_summary), dirname(paths$confusion_metrics), dirname(paths$outcome_definitions), dirname(paths$roc_data), dirname(paths$roc_plot), dirname(paths$predictions_full), dirname(paths$train_eval), dirname(paths$test_eval), dirname(paths$prediction_map), dirname(paths$train_outcome_map), dirname(paths$test_outcome_map), dirname(paths$missingness_csv), dirname(paths$corr_csv), dirname(paths$variance_csv), dirname(paths$outlier_csv), dirname(paths$class_balance_csv), paths$whitebox_dir))
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

validate_unique_id <- function(x, label) {
  if (!"id" %in% names(x)) stop("Object ", label, " is missing `id`.")
  dup_n <- sum(duplicated(x$id))
  if (dup_n > 0) stop("Object ", label, " has duplicated `id` values: ", dup_n)
  invisible(TRUE)
}

read_target_fishnet <- function(path) {
  if (!file.exists(path)) stop("Target-valid fishnet not found: ", path)
  st_read(path, quiet = TRUE)
}

read_feature_tables <- function(paths_list) {
  out <- list()
  for (nm in names(paths_list)) {
    path <- paths_list[[nm]]
    if (!file.exists(path)) stop("Missing feature CSV for ", nm, ": ", path)
    tbl <- read_csv(path, show_col_types = FALSE, progress = FALSE)
    validate_unique_id(tbl, paste0("feature table: ", nm))
    out[[nm]] <- tbl
  }
  out
}

compute_missingness_summary <- function(df) {
  dplyr::tibble(
    variable = names(df),
    n_missing = vapply(df, function(x) sum(is.na(x)), numeric(1)),
    pct_missing = vapply(df, function(x) mean(is.na(x)), numeric(1))
  ) %>% arrange(desc(n_missing), variable)
}

compute_class_balance <- function(df, target_col) {
  if (!target_col %in% names(df)) stop("Class balance summary requires target column: ", target_col)
  total_n <- nrow(df)
  value_counts <- sort(table(df[[target_col]]), decreasing = TRUE)
  dplyr::tibble(
    target_value = as.integer(names(value_counts)),
    n = as.integer(value_counts),
    pct = as.integer(value_counts) / total_n,
    class_label = ifelse(as.integer(names(value_counts)) == 1, "Inundated", "Not inundated")
  ) %>%
    select(class_label, target_value, n, pct) %>%
    arrange(desc(target_value))
}

compute_variance_check <- function(df, predictors) {
  dplyr::tibble(
    variable = predictors,
    n_unique = vapply(predictors, function(v) length(unique(df[[v]][!is.na(df[[v]])])), integer(1)),
    variance = vapply(predictors, function(v) stats::var(df[[v]], na.rm = TRUE), numeric(1))
  ) %>%
    mutate(
      zero_variance = is.na(variance) | variance == 0,
      unique_prop = n_unique / nrow(df),
      near_zero_variance = zero_variance | n_unique <= 2 | unique_prop <= 0.05,
      drop_flag = zero_variance | near_zero_variance
    ) %>% arrange(desc(drop_flag), variable)
}

compute_correlation_matrix <- function(df, predictors) {
  if (length(predictors) == 0) stop("No predictors available for correlation matrix.")
  if (length(predictors) == 1) {
    mat <- matrix(1, nrow = 1, ncol = 1, dimnames = list(predictors, predictors))
  } else {
    mat <- stats::cor(df[, predictors, drop = FALSE], use = "pairwise.complete.obs")
  }
  wide_df <- as.data.frame(mat)
  wide_df <- tibble::rownames_to_column(wide_df, var = "variable")
  long_df <- as.data.frame(as.table(mat))
  names(long_df) <- c("x", "y", "correlation")
  list(wide_df = wide_df, long_df = long_df)
}

save_correlation_heatmap <- function(long_df, path) {
  p <- ggplot(long_df, aes(x = x, y = y, fill = correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
    coord_fixed() +
    labs(title = "Correlation matrix for retained numeric predictors", x = NULL, y = NULL, fill = "r") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(path, plot = p, width = 9, height = 7, dpi = 300)
}

compute_outlier_summary <- function(df, predictors, iqr_multiplier = 1.5) {
  rows <- lapply(predictors, function(v) {
    x <- df[[v]]
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(dplyr::tibble(variable = v, n_outliers = 0, pct_outliers = NA_real_, lower_fence = NA_real_, upper_fence = NA_real_, q1 = NA_real_, q3 = NA_real_, min = NA_real_, max = NA_real_))
    }
    q1 <- as.numeric(stats::quantile(x, 0.25, names = FALSE))
    q3 <- as.numeric(stats::quantile(x, 0.75, names = FALSE))
    iqr <- q3 - q1
    lower <- q1 - iqr_multiplier * iqr
    upper <- q3 + iqr_multiplier * iqr
    n_out <- sum(x < lower | x > upper)
    dplyr::tibble(
      variable = v,
      n_outliers = n_out,
      pct_outliers = n_out / length(x),
      lower_fence = lower,
      upper_fence = upper,
      q1 = q1,
      q3 = q3,
      min = min(x),
      max = max(x)
    )
  })
  bind_rows(rows) %>% arrange(desc(n_outliers), variable)
}

save_histograms <- function(df, vars, out_dir, exclude = character(), target_col = "target") {
  vars <- setdiff(vars, exclude)
  for (v in vars) {
    if (!v %in% names(df)) next
    p <- if (v == target_col) {
      ggplot(df, aes(x = factor(.data[[v]]))) +
        geom_bar(fill = "#2b8cbe") +
        labs(title = "Target class balance", x = "Target", y = "Count") +
        theme_minimal(base_size = 11)
    } else {
      ggplot(df, aes(x = .data[[v]])) +
        geom_histogram(bins = 30, fill = "#2b8cbe", color = "white") +
        labs(title = paste("Histogram:", v), x = v, y = "Count") +
        theme_minimal(base_size = 11)
    }
    safe_name <- gsub("[^A-Za-z0-9_]+", "_", v)
    ggsave(file.path(out_dir, paste0("hist_", safe_name, ".png")), plot = p, width = 7, height = 5, dpi = 300)
  }
}

is_binaryish <- function(x) {
  x <- x[!is.na(x)]
  is.numeric(x) && length(x) > 0 && all(sort(unique(x)) %in% c(0, 1))
}

make_sf_map <- function(sfobj, var, title = NULL) {
  p <- ggplot(sfobj) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    labs(title = title %||% var, fill = var) +
    theme_void(base_size = 11) +
    theme(legend.position = "right")
  if (is_binaryish(sfobj[[var]])) {
    p + scale_fill_manual(values = c("0" = "#2b8cbe", "1" = "#de2d26"), na.value = "grey90")
  } else {
    p + scale_fill_gradient(low = "#edf8fb", high = "#006d2c", na.value = "grey90")
  }
}

save_selected_maps <- function(sfobj, out_dir, priority_vars) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  vars <- priority_vars[priority_vars %in% names(sfobj)]
  for (v in vars) {
    ggsave(file.path(out_dir, paste0("map_", v, ".png")), plot = make_sf_map(sfobj, v, paste("Selected map:", v)), width = 7.5, height = 6.5, dpi = 300)
  }
}

save_binary_target_map <- function(sfobj, path, target_col) {
  if (!target_col %in% names(sfobj)) return(invisible(NULL))
  p <- ggplot(sfobj) +
    geom_sf(aes(fill = factor(.data[[target_col]])), color = NA) +
    scale_fill_manual(values = c("0" = "#2b8cbe", "1" = "#de2d26"), name = "Target") +
    labs(title = "Observed binary inundation target") +
    theme_void(base_size = 11) +
    theme(legend.position = "right")
  ggsave(path, plot = p, width = 7.5, height = 6.5, dpi = 300)
}

save_target_valid_map <- function(sfobj, path) {
  if (!"inundated_share" %in% names(sfobj)) return(invisible(NULL))
  p <- ggplot(sfobj) +
    geom_sf(aes(fill = inundated_share), color = NA) +
    scale_fill_gradient(low = "#edf8fb", high = "#006d2c", na.value = "grey90") +
    labs(title = "Inundated share within fishnet cells", fill = "Share") +
    theme_void(base_size = 11) +
    theme(legend.position = "right")
  ggsave(path, plot = p, width = 7.5, height = 6.5, dpi = 300)
}

merge_and_clean_features <- function(config = list()) {
  config <- normalize_merge_eda_config(config)
  ensure_merge_eda_dirs(config$paths)

  target_col <- config$target_col
  if (is.null(target_col) || !nzchar(target_col)) target_col <- "target"

  fishnet_path <- resolve_path(config$paths, "fishnet_target", "fishnet_base")
  if (is.null(fishnet_path)) stop("No fishnet path found in config$paths.")

  feature_paths <- list(
    topography = config$paths$features_topography,
    hydrology = config$paths$features_hydrology,
    watershed = config$paths$features_watershed,
    landcover = config$paths$features_landcover,
    built = config$paths$features_built
  )

  fishnet_sf <- read_target_fishnet(fishnet_path)
  validate_unique_id(fishnet_sf, "target-valid fishnet")
  message("Target-valid fishnet rows: ", nrow(fishnet_sf))
  message("Target-valid fishnet CRS: ", sf::st_crs(fishnet_sf)$input)

  feature_tables <- read_feature_tables(feature_paths)
  fishnet_features_sf <- fishnet_sf
  for (nm in names(feature_tables)) {
    ft <- feature_tables[[nm]]
    message("Joining feature table: ", nm, " (rows = ", nrow(ft), ")")
    fishnet_features_sf <- left_join(fishnet_features_sf, ft, by = "id")
    validate_unique_id(fishnet_features_sf, paste0("fishnet after join: ", nm))
  }

  st_write(fishnet_features_sf, dsn = config$paths$features_all, delete_dsn = TRUE, quiet = TRUE)
  message("Saved merged geospatial feature layer: ", config$paths$features_all)

  modeling_raw <- st_drop_geometry(fishnet_features_sf)
  name_order <- c("id", target_col, "inundated_share", setdiff(names(modeling_raw), c("id", target_col, "inundated_share")))
  name_order <- name_order[name_order %in% names(modeling_raw)]
  modeling_raw <- modeling_raw[, name_order, drop = FALSE]
  modeling_raw$id <- as.integer(modeling_raw$id)
  modeling_raw[[target_col]] <- as.integer(modeling_raw[[target_col]])
  modeling_raw$inundated_share <- as.numeric(modeling_raw$inundated_share)
  modeling_raw <- modeling_raw[modeling_raw[[target_col]] %in% c(0L, 1L), , drop = FALSE]

  write_csv(modeling_raw, config$paths$modeling_raw, na = "")
  message("Saved raw modeling table: ", config$paths$modeling_raw)

  missingness_summary <- compute_missingness_summary(modeling_raw)
  write_csv(missingness_summary, config$paths$missingness_csv, na = "")

  class_balance <- compute_class_balance(modeling_raw, target_col)
  write_csv(class_balance, config$paths$class_balance_csv, na = "")

  if (!is.null(config$predictor_cols) && length(config$predictor_cols) > 0) {
    candidate_predictors <- intersect(config$predictor_cols, names(modeling_raw))
  } else {
    candidate_predictors <- setdiff(names(modeling_raw), c("id", target_col, "inundated_share"))
  }
  candidate_predictors <- candidate_predictors[sapply(modeling_raw[candidate_predictors], is.numeric)]
  if (length(candidate_predictors) == 0) stop("No numeric predictors found after merging feature tables.")

  variance_check <- compute_variance_check(modeling_raw, candidate_predictors)
  write_csv(variance_check, config$paths$variance_csv, na = "")
  drop_predictors <- variance_check$variable[variance_check$drop_flag]
  retained_predictors <- setdiff(candidate_predictors, drop_predictors)
  message("Candidate predictors: ", length(candidate_predictors))
  message("Dropped zero/near-zero variance predictors: ", length(drop_predictors))
  if (length(drop_predictors) > 0) message("Dropped predictors: ", paste(drop_predictors, collapse = ", "))

  corr_predictors <- retained_predictors[sapply(modeling_raw[retained_predictors], function(x) stats::sd(x, na.rm = TRUE) > 0)]
  if (length(corr_predictors) > 0) {
    corr_out <- compute_correlation_matrix(modeling_raw, corr_predictors)
    write_csv(corr_out$wide_df, config$paths$corr_csv, na = "")
    save_correlation_heatmap(corr_out$long_df, config$paths$corr_png)
  } else {
    corr_out <- list(wide_df = NULL, long_df = NULL)
  }

  outlier_summary <- compute_outlier_summary(modeling_raw, retained_predictors, config$eda$iqr_multiplier)
  write_csv(outlier_summary, config$paths$outlier_csv, na = "")

  save_histograms(
    modeling_raw,
    vars = c(target_col, retained_predictors),
    out_dir = config$paths$eda_fig_dir,
    exclude = config$eda$histogram_exclude,
    target_col = target_col
  )
  save_selected_maps(fishnet_features_sf, config$paths$eda_fig_dir, config$eda$map_priority)
  save_binary_target_map(fishnet_features_sf, config$paths$target_map_png, target_col)
  save_target_valid_map(fishnet_features_sf, config$paths$target_valid_map_png)

  modeling_raw_for_clean <- modeling_raw[!is.na(modeling_raw[[target_col]]) & modeling_raw[[target_col]] %in% c(0L, 1L), , drop = FALSE]
  retained_predictors <- retained_predictors[retained_predictors %in% names(modeling_raw_for_clean)]
  retained_predictors <- retained_predictors[sapply(modeling_raw_for_clean[retained_predictors], is.numeric)]
  if (length(retained_predictors) == 0) stop("No predictors remain after variance screening.")

  keep_rows <- complete.cases(modeling_raw_for_clean[, c(target_col, retained_predictors), drop = FALSE])
  modeling_clean <- modeling_raw_for_clean[keep_rows, c("id", target_col, "inundated_share", retained_predictors), drop = FALSE]

  final_variance_check <- compute_variance_check(modeling_clean, retained_predictors)
  final_drop <- final_variance_check$variable[final_variance_check$drop_flag]
  if (length(final_drop) > 0) {
    retained_predictors <- setdiff(retained_predictors, final_drop)
    if (length(retained_predictors) == 0) stop("All predictors were removed after complete-case filtering.")
    modeling_clean <- modeling_clean[, c("id", target_col, "inundated_share", retained_predictors), drop = FALSE]
    message("Dropped predictors after complete-case filtering: ", paste(final_drop, collapse = ", "))
  }

  if (nrow(modeling_clean) == 0) stop("Cleaning removed all rows from the modeling table.")
  write_csv(modeling_clean, config$paths$modeling_clean, na = "")
  message("Saved cleaned modeling table: ", config$paths$modeling_clean)

  message("Modeling rows before cleaning: ", nrow(modeling_raw))
  message("Modeling rows after cleaning: ", nrow(modeling_clean))
  message("Final retained predictors: ", paste(retained_predictors, collapse = ", "))

  invisible(list(
    fishnet_features_sf = fishnet_features_sf,
    modeling_raw = modeling_raw,
    modeling_clean = modeling_clean,
    feature_tables = feature_tables,
    missingness_summary = missingness_summary,
    class_balance = class_balance,
    variance_check = variance_check,
    final_variance_check = final_variance_check,
    correlation_matrix = corr_out$wide_df,
    outlier_summary = outlier_summary,
    retained_predictors = retained_predictors,
    dropped_predictors = drop_predictors,
    paths = config$paths,
    target_col = target_col
  ))
}

