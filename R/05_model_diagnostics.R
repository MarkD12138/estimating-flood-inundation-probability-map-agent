run_glm_diagnostics <- function(config) {
  # Diagnostics and glm modeling for the Calgary flood workflow.
  # Expected config fields:
  # - config$project_root
  # - config$seed
  # - config$classification_threshold
  # - config$target_col
  # - config$predictor_cols
  # - config$analysis_crs_input
  # - config$paths with the keys documented in the user contract

  ensure_libs <- function() {
    # Prefer a project-local library if one exists.
    local_lib <- file.path(getwd(), "r_libs")
    if (dir.exists(local_lib)) {
      .libPaths(c(normalizePath(local_lib), .libPaths()))
    }

    required <- c(
      "sf", "dplyr", "readr", "ggplot2", "pROC", "car", "broom",
      "tidyr", "stringr", "forcats", "janitor", "scales"
    )
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) {
      stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  clean_name <- function(x) {
    if (is.null(x)) return(NULL)
    janitor::make_clean_names(x)
  }

  is_abs_path <- function(path) {
    grepl("^(?:[A-Za-z]:[\\\\/]|/)", path)
  }

  resolve_path <- function(path, project_root = NULL) {
    if (is.null(path) || identical(path, "")) return(path)
    if (!is.null(project_root) && !is_abs_path(path)) {
      file.path(project_root, path)
    } else {
      path
    }
  }

  get_path <- function(x, ..., default = NULL) {
    val <- x
    for (nm in list(...)) {
      if (is.null(val) || !is.list(val) || is.null(val[[nm]])) {
        val <- NULL
        break
      }
      val <- val[[nm]]
    }
    if (is.null(val) || identical(val, "")) default else val
  }

  ensure_dirs <- function(paths) {
    for (p in paths) {
      d <- dirname(p)
      if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(TRUE)
  }

  write_csv_safe <- function(x, path) {
    readr::write_csv(x, path, na = "")
  }

  save_gg <- function(plot, path, width = 10, height = 7, dpi = 300) {
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      units = "in",
      bg = "white"
    )
  }

  standardize_df <- function(df, params) {
    out <- df
    for (i in seq_len(nrow(params))) {
      p <- params$predictor[i]
      out[[p]] <- (out[[p]] - params$mean[i]) / params$sd[i]
    }
    out
  }

  drop_zero_variance <- function(df, predictors) {
    keep <- predictors[vapply(predictors, function(p) {
      vals <- df[[p]]
      vals <- vals[!is.na(vals)]
      length(unique(vals)) > 1 && stats::sd(vals) > 0
    }, logical(1))]
    list(keep = keep, dropped = setdiff(predictors, keep))
  }

  reduce_high_correlation <- function(df, predictors, protected = character(0), threshold = 0.90) {
    if (length(predictors) < 2) {
      return(list(keep = predictors, dropped = character(0), cor = NULL))
    }

    mat <- stats::cor(df[, predictors, drop = FALSE], use = "pairwise.complete.obs")
    keep <- predictors
    dropped <- character(0)

    repeat {
      if (length(keep) < 2) break
      submat <- mat[keep, keep, drop = FALSE]
      diag(submat) <- 0
      mx <- suppressWarnings(max(abs(submat), na.rm = TRUE))
      if (!is.finite(mx) || mx <= threshold) break

      idx <- which(abs(submat) == mx, arr.ind = TRUE)[1, ]
      pair <- rownames(submat)[idx]
      cand <- pair[!pair %in% protected]
      if (length(cand) == 0) break

      if (length(cand) == 2) {
        scores <- vapply(cand, function(p) mean(abs(submat[p, setdiff(keep, p)]), na.rm = TRUE), numeric(1))
        drop <- cand[which.max(scores)]
      } else {
        drop <- cand[1]
      }

      keep <- setdiff(keep, drop)
      dropped <- unique(c(dropped, drop))
    }

    list(keep = keep, dropped = dropped, cor = mat)
  }

  select_forced_watershed <- function(df, watershed_vars) {
    watershed_vars <- intersect(watershed_vars, names(df))
    watershed_vars <- watershed_vars[vapply(watershed_vars, function(p) {
      vals <- df[[p]]
      vals <- vals[!is.na(vals)]
      length(unique(vals)) > 1 && stats::sd(vals) > 0
    }, logical(1))]
    if (length(watershed_vars) == 0) return(character(0))
    if (length(watershed_vars) == 1) return(watershed_vars)

    scores <- vapply(watershed_vars, function(p) {
      fit <- stats::glm(stats::as.formula(paste("target ~", p)), data = df, family = stats::binomial())
      stats::AIC(fit)
    }, numeric(1))
    watershed_vars[which.min(scores)][1]
  }

  calc_classification_metrics <- function(truth, prob, threshold = 0.5) {
    truth <- as.integer(truth)
    pred <- as.integer(prob >= threshold)

    tp <- sum(pred == 1L & truth == 1L, na.rm = TRUE)
    tn <- sum(pred == 0L & truth == 0L, na.rm = TRUE)
    fp <- sum(pred == 1L & truth == 0L, na.rm = TRUE)
    fn <- sum(pred == 0L & truth == 1L, na.rm = TRUE)

    accuracy <- if ((tp + tn + fp + fn) > 0) (tp + tn) / (tp + tn + fp + fn) else NA_real_
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    f1 <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA_real_
    }

    data.frame(
      TP = tp,
      TN = tn,
      FP = fp,
      FN = fn,
      Accuracy = accuracy,
      Precision = precision,
      Recall = recall,
      F1 = f1,
      Specificity = specificity,
      stringsAsFactors = FALSE
    )
  }

  plot_roc <- function(roc_df, auc_value, path) {
    p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
      ggplot2::geom_line(color = "#2c7fb8", linewidth = 1) +
      ggplot2::geom_abline(linetype = "dashed", color = "grey50") +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = paste0("ROC Curve (AUC = ", formatC(auc_value, digits = 3, format = "f"), ")"),
        x = "False positive rate",
        y = "True positive rate"
      ) +
      ggplot2::theme_minimal(base_size = 12)
    save_gg(p, path, width = 7, height = 7)
  }

  plot_prob_map <- function(sf_obj, path, title) {
    p <- ggplot2::ggplot(sf_obj) +
      ggplot2::geom_sf(ggplot2::aes(fill = prob), color = NA, linewidth = 0) +
      ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90", name = "Predicted\nprobability") +
      ggplot2::labs(title = title, caption = "Calgary flood inundation probability") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    save_gg(p, path, width = 9, height = 8)
  }

  plot_outcome_map <- function(sf_obj, path, title) {
    outcome_cols <- c(TP = "#1b9e77", TN = "#7570b3", FP = "#d95f02", FN = "#e7298a")
    p <- ggplot2::ggplot(sf_obj) +
      ggplot2::geom_sf(ggplot2::aes(fill = outcome), color = NA, linewidth = 0) +
      ggplot2::scale_fill_manual(values = outcome_cols, drop = FALSE, na.value = "grey90", name = "Outcome") +
      ggplot2::labs(title = title, caption = "Threshold = 0.5") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    save_gg(p, path, width = 9, height = 8)
  }

  ensure_libs()
  sf::sf_use_s2(FALSE)

  project_root <- config$project_root %||% getwd()
  paths <- config$paths %||% list()
  model_cfg <- config$model %||% list()
  seed <- config$seed %||% 20250322L
  threshold <- config$classification_threshold %||% 0.5
  target_col <- clean_name(config$target_col %||% "target")
  predictor_cols <- clean_name(config$predictor_cols %||% NULL)
  analysis_crs_input <- config$analysis_crs_input %||% NULL

  cleaned_modeling_csv <- resolve_path(get_path(paths, "modeling_clean", default = get_path(paths, "cleaned_modeling_table", default = get_path(paths, "modeling_table_clean", default = "outputs/intermediate/modeling_table_clean.csv"))), project_root)
  merged_features_gpkg <- resolve_path(get_path(paths, "features_all", default = get_path(paths, "merged_feature_gpkg", default = get_path(paths, "merged_features_gpkg", default = get_path(paths, "fishnet_features_gpkg", default = "outputs/intermediate/features_all.gpkg")))), project_root)
  diagnostics_dir <- resolve_path(get_path(paths, "diagnostics_dir", default = "outputs/diagnostics"), project_root)
  models_dir <- resolve_path(get_path(paths, "models_dir", default = "outputs/models"), project_root)
  figures_dir <- resolve_path(get_path(paths, "figures_dir", default = "outputs/figures"), project_root)
  diag_fig_dir <- resolve_path(get_path(paths, "diag_fig_dir", default = get_path(paths, "diagnostics_fig_dir", default = file.path(figures_dir, "diagnostics"))), project_root)

  model_rds <- resolve_path(get_path(paths, "model_rds", default = file.path(models_dir, "glm_logistic_calgary.rds")), project_root)
  scaling_csv <- resolve_path(get_path(paths, "scaling_params", default = get_path(paths, "scaling_params_csv", default = file.path(diagnostics_dir, "scaling_parameters.csv"))), project_root)
  model_summary_csv <- resolve_path(get_path(paths, "model_summary", default = get_path(paths, "model_summary_csv", default = file.path(diagnostics_dir, "model_summary.csv"))), project_root)
  confusion_csv <- resolve_path(get_path(paths, "confusion_metrics", default = get_path(paths, "confusion_metrics_csv", default = file.path(diagnostics_dir, "confusion_metrics.csv"))), project_root)
  tpdef_csv <- resolve_path(get_path(paths, "outcome_definitions", default = get_path(paths, "tp_tn_fp_fn_definitions_csv", default = file.path(diagnostics_dir, "tp_tn_fp_fn_definitions.csv"))), project_root)
  roc_csv <- resolve_path(get_path(paths, "roc_data", default = get_path(paths, "roc_csv", default = file.path(diagnostics_dir, "roc_curve.csv"))), project_root)
  roc_png <- resolve_path(get_path(paths, "roc_plot", default = get_path(paths, "roc_png", default = file.path(diag_fig_dir, "roc_curve.png"))), project_root)
  prediction_map <- resolve_path(get_path(paths, "prediction_map", default = get_path(paths, "probability_map_png", default = file.path(diag_fig_dir, "calgary_prediction_probability.png"))), project_root)
  train_outcome_png <- resolve_path(get_path(paths, "train_outcome_map", default = file.path(diag_fig_dir, "train_outcome_map.png")), project_root)
  test_outcome_png <- resolve_path(get_path(paths, "test_outcome_map", default = file.path(diag_fig_dir, "test_outcome_map.png")), project_root)
  predictions_full_path <- resolve_path(get_path(paths, "predictions_full", default = file.path(diagnostics_dir, "calgary_prediction_probabilities.csv")), project_root)
  train_eval_path <- resolve_path(get_path(paths, "train_eval", default = file.path(diagnostics_dir, "train_eval.csv")), project_root)
  test_eval_path <- resolve_path(get_path(paths, "test_eval", default = file.path(diagnostics_dir, "test_eval.csv")), project_root)

  ensure_dirs(c(
    dirname(cleaned_modeling_csv), dirname(merged_features_gpkg), dirname(model_rds),
    dirname(scaling_csv), dirname(model_summary_csv), dirname(confusion_csv),
    dirname(tpdef_csv), dirname(roc_csv), dirname(roc_png), dirname(prediction_map),
    dirname(train_outcome_png), dirname(test_outcome_png), dirname(predictions_full_path),
    dirname(train_eval_path), dirname(test_eval_path)
  ))

  message("Reading modeling data and geospatial features...")
  modeling_df <- readr::read_csv(cleaned_modeling_csv, show_col_types = FALSE) |>
    janitor::clean_names()
  feature_sf <- sf::st_read(merged_features_gpkg, quiet = TRUE) |>
    janitor::clean_names()

  if (!target_col %in% names(modeling_df)) {
    stop("Cleaned modeling CSV must contain target column `", target_col, "`.", call. = FALSE)
  }
  if (!"id" %in% names(modeling_df)) stop("Cleaned modeling CSV must contain `id`.", call. = FALSE)
  if (!"id" %in% names(feature_sf)) stop("Merged feature GPKG must contain `id`.", call. = FALSE)

  modeling_df[[target_col]] <- as.integer(modeling_df[[target_col]])

  if (!is.null(predictor_cols) && length(predictor_cols) > 0) {
    missing_predictors <- setdiff(predictor_cols, names(modeling_df))
    if (length(missing_predictors) > 0) {
      stop("Configured predictor columns are missing from the cleaned modeling table: ", paste(missing_predictors, collapse = ", "), call. = FALSE)
    }
    numeric_predictors <- predictor_cols
  } else {
    numeric_predictors <- names(modeling_df)[vapply(modeling_df, is.numeric, logical(1))]
    numeric_predictors <- setdiff(numeric_predictors, c("id", target_col))
  }
  if (length(numeric_predictors) == 0) stop("No numeric predictors found in the cleaned modeling table.", call. = FALSE)

  watershed_candidates <- intersect(c("max_flow_accum", "upstream_ws_area"), numeric_predictors)
  if (length(watershed_candidates) == 0) {
    watershed_candidates <- intersect(numeric_predictors[stringr::str_detect(numeric_predictors, "flow|ws")], numeric_predictors)
  }

  if (!is.null(analysis_crs_input)) {
    crs_target <- tryCatch(sf::st_crs(analysis_crs_input), error = function(e) NULL)
    if (!is.null(crs_target) && !is.na(crs_target$wkt)) {
      feature_sf <- sf::st_transform(feature_sf, crs_target)
    }
  }

  model_vars <- c(target_col, numeric_predictors)
  modeling_df <- modeling_df |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(model_vars), ~ !is.na(.x)))

  feature_sf <- feature_sf |>
    dplyr::filter(.data$id %in% modeling_df$id)

  modeling_sf <- feature_sf |>
    dplyr::left_join(modeling_df, by = "id")

  if (any(is.na(modeling_sf[[target_col]]))) {
    stop("Target is missing after the join between geometry and modeling data.", call. = FALSE)
  }

  message("Rows in modeling data after complete-case filter: ", nrow(modeling_df))
  message("Rows in geospatial modeling layer: ", nrow(modeling_sf))
  message("Target balance:")
  print(table(modeling_df[[target_col]]))

  geospatial_modeling_gpkg <- resolve_path(file.path(dirname(merged_features_gpkg), "modeling_table_clean.gpkg"), project_root)
  sf::st_write(modeling_sf, geospatial_modeling_gpkg, delete_dsn = TRUE, quiet = TRUE)

  set.seed(seed)
  train_id <- sample(seq_len(nrow(modeling_df)), size = floor(0.7 * nrow(modeling_df)))
  train_df <- modeling_df[train_id, , drop = FALSE]
  test_df <- modeling_df[-train_id, , drop = FALSE]

  message("Train rows: ", nrow(train_df), " | Test rows: ", nrow(test_df))

  drop_zero <- drop_zero_variance(train_df, numeric_predictors)
  if (length(drop_zero$dropped) > 0) {
    message("Dropping zero-variance predictors: ", paste(drop_zero$dropped, collapse = ", "))
  }
  candidate_predictors <- drop_zero$keep

  corr_filter <- reduce_high_correlation(
    train_df,
    candidate_predictors,
    protected = intersect(candidate_predictors, watershed_candidates),
    threshold = 0.90
  )
  if (length(corr_filter$dropped) > 0) {
    message("Dropping highly correlated predictors: ", paste(corr_filter$dropped, collapse = ", "))
  }
  candidate_predictors <- corr_filter$keep

  forced_watershed <- select_forced_watershed(train_df, watershed_candidates)
  if (length(forced_watershed) == 0) {
    message("No estimable watershed predictor found; continuing without a forced watershed term.")
  } else {
    message("Forced watershed predictor retained for model selection: ", paste(forced_watershed, collapse = ", "))
    candidate_predictors <- union(candidate_predictors, forced_watershed)
  }

  selected_candidates <- candidate_predictors
  if (length(selected_candidates) == 0) stop("No predictors remain after variance and correlation screening.", call. = FALSE)

  scaling_params <- dplyr::tibble(
    predictor = selected_candidates,
    mean = vapply(selected_candidates, function(p) mean(train_df[[p]], na.rm = TRUE), numeric(1)),
    sd = vapply(selected_candidates, function(p) stats::sd(train_df[[p]], na.rm = TRUE), numeric(1))
  ) |>
    dplyr::mutate(sd = ifelse(sd == 0 | is.na(sd), 1, sd))
  write_csv_safe(scaling_params, scaling_csv)

  train_std <- standardize_df(train_df, scaling_params)
  test_std <- standardize_df(test_df, scaling_params)
  full_std <- standardize_df(modeling_df, scaling_params)

  upper_formula <- stats::as.formula(paste(target_col, "~", paste(selected_candidates, collapse = " + ")))
  lower_terms <- if (length(forced_watershed) > 0) forced_watershed else character(0)
  lower_formula <- if (length(lower_terms) > 0) {
    stats::as.formula(paste(target_col, "~", paste(lower_terms, collapse = " + ")))
  } else {
    stats::as.formula(paste(target_col, "~ 1"))
  }

  full_model <- stats::glm(upper_formula, data = train_std, family = stats::binomial())
  aliased <- stats::alias(full_model)$Complete
  if (!is.null(aliased) && length(aliased) > 0) {
    aliased_names <- names(aliased)
    message("Removing aliased predictors: ", paste(aliased_names, collapse = ", "))
    selected_candidates <- setdiff(selected_candidates, aliased_names)
    upper_formula <- stats::as.formula(paste(target_col, "~", paste(selected_candidates, collapse = " + ")))
    full_model <- stats::glm(upper_formula, data = train_std, family = stats::binomial())
  }

  vif_threshold <- model_cfg$vif_threshold %||% 10
  active_predictors <- selected_candidates
  repeat {
    if (length(active_predictors) <= 1) break
    vif_fit <- stats::glm(stats::as.formula(paste(target_col, "~", paste(active_predictors, collapse = " + "))), data = train_std, family = stats::binomial())
    vif_vals <- tryCatch(car::vif(vif_fit), error = function(e) NULL)
    if (is.null(vif_vals) || length(vif_vals) == 0) break
    vif_vals <- sort(vif_vals, decreasing = TRUE)
    if (max(vif_vals, na.rm = TRUE) <= vif_threshold) break

    drop_cand <- names(vif_vals)[1]
    if (drop_cand %in% lower_terms) {
      vif_vals <- vif_vals[setdiff(names(vif_vals), lower_terms)]
      if (length(vif_vals) == 0 || max(vif_vals, na.rm = TRUE) <= vif_threshold) break
      drop_cand <- names(vif_vals)[1]
    }
    message("Dropping high-VIF predictor: ", drop_cand, " (VIF=", round(vif_vals[drop_cand], 2), ")")
    active_predictors <- setdiff(active_predictors, drop_cand)
  }
  if (length(active_predictors) == 0) stop("All predictors were removed during VIF screening.", call. = FALSE)

  step_upper_formula <- stats::as.formula(paste(target_col, "~", paste(active_predictors, collapse = " + ")))
  step_model <- stats::glm(step_upper_formula, data = train_std, family = stats::binomial())
  final_glm <- suppressWarnings(
    stats::step(
      step_model,
      scope = list(lower = lower_formula, upper = step_upper_formula),
      direction = "backward",
      trace = 0
    )
  )

  final_terms <- attr(stats::terms(final_glm), "term.labels")
  if (length(intersect(final_terms, lower_terms)) == 0 && length(lower_terms) > 0) {
    final_terms <- unique(c(lower_terms[1], final_terms))
    final_formula <- stats::as.formula(paste(target_col, "~", paste(final_terms, collapse = " + ")))
    final_glm <- stats::glm(final_formula, data = train_std, family = stats::binomial())
  }

  message("Final model terms: ", paste(attr(stats::terms(final_glm), "term.labels"), collapse = ", "))

  model_summary <- broom::tidy(final_glm) |>
    dplyr::mutate(model_auc = as.numeric(pROC::auc(pROC::roc(train_std[[target_col]], stats::predict(final_glm, type = "response"), quiet = TRUE, levels = c(0, 1), direction = "<"))))
  write_csv_safe(model_summary, model_summary_csv)
  saveRDS(final_glm, model_rds)

  train_prob <- stats::predict(final_glm, newdata = train_std, type = "response")
  test_prob <- stats::predict(final_glm, newdata = test_std, type = "response")
  full_prob <- stats::predict(final_glm, newdata = full_std, type = "response")

  train_pred <- train_std |>
    dplyr::mutate(prob = train_prob, predicted = as.integer(prob >= threshold))
  test_pred <- test_std |>
    dplyr::mutate(prob = test_prob, predicted = as.integer(prob >= threshold))
  full_pred <- full_std |>
    dplyr::mutate(prob = full_prob, predicted = as.integer(prob >= threshold))

  train_pred <- train_pred |>
    dplyr::mutate(outcome = dplyr::case_when(
      .data[[target_col]] == 1L & predicted == 1L ~ "TP",
      .data[[target_col]] == 0L & predicted == 0L ~ "TN",
      .data[[target_col]] == 0L & predicted == 1L ~ "FP",
      .data[[target_col]] == 1L & predicted == 0L ~ "FN",
      TRUE ~ NA_character_
    ))
  test_pred <- test_pred |>
    dplyr::mutate(outcome = dplyr::case_when(
      .data[[target_col]] == 1L & predicted == 1L ~ "TP",
      .data[[target_col]] == 0L & predicted == 0L ~ "TN",
      .data[[target_col]] == 0L & predicted == 1L ~ "FP",
      .data[[target_col]] == 1L & predicted == 0L ~ "FN",
      TRUE ~ NA_character_
    ))
  full_pred <- full_pred |>
    dplyr::mutate(outcome = dplyr::case_when(
      .data[[target_col]] == 1L & predicted == 1L ~ "TP",
      .data[[target_col]] == 0L & predicted == 0L ~ "TN",
      .data[[target_col]] == 0L & predicted == 1L ~ "FP",
      .data[[target_col]] == 1L & predicted == 0L ~ "FN",
      TRUE ~ NA_character_
    ))

  write_csv_safe(train_pred, train_eval_path)
  write_csv_safe(test_pred, test_eval_path)
  write_csv_safe(full_pred |> dplyr::select(id, dplyr::all_of(target_col), prob, predicted, outcome), predictions_full_path)

  roc_obj <- pROC::roc(response = test_pred[[target_col]], predictor = test_pred$prob, quiet = TRUE, levels = c(0, 1), direction = "<")
  roc_df <- data.frame(
    specificity = pROC::specificities(roc_obj),
    sensitivity = pROC::sensitivities(roc_obj)
  )
  write_csv_safe(roc_df, roc_csv)
  plot_roc(roc_df, as.numeric(pROC::auc(roc_obj)), roc_png)

  confusion_tbl <- calc_classification_metrics(test_pred[[target_col]], test_pred$prob, threshold = threshold)
  write_csv_safe(confusion_tbl, confusion_csv)

  tp_defs <- data.frame(
    outcome = c("TP", "TN", "FP", "FN"),
    definition = c(
      "Observed flood cell predicted as flood.",
      "Observed dry cell predicted as dry.",
      "Observed dry cell predicted as flood.",
      "Observed flood cell predicted as dry."
    ),
    stringsAsFactors = FALSE
  )
  write_csv_safe(tp_defs, tpdef_csv)

  model_map_sf <- modeling_sf |>
    dplyr::left_join(full_pred |> dplyr::select(id, prob, predicted, outcome), by = "id")

  train_map_sf <- model_map_sf |>
    dplyr::filter(.data$id %in% train_pred$id)
  test_map_sf <- model_map_sf |>
    dplyr::filter(.data$id %in% test_pred$id)

  plot_prob_map(model_map_sf, prediction_map, "Calgary predicted inundation probability")
  plot_outcome_map(train_map_sf, train_outcome_png, "Train-set classification outcomes")
  plot_outcome_map(test_map_sf, test_outcome_png, "Test-set classification outcomes")

  prediction_summary <- data.frame(
    dataset = c("train", "test", "full"),
    n = c(nrow(train_pred), nrow(test_pred), nrow(full_pred)),
    mean_prob = c(mean(train_prob), mean(test_prob), mean(full_prob)),
    stringsAsFactors = FALSE
  )
  write_csv_safe(prediction_summary, file.path(diagnostics_dir, "prediction_summary.csv"))

  invisible(list(
    model = final_glm,
    scaling_params = scaling_params,
    confusion_metrics = confusion_tbl,
    roc = roc_obj,
    roc_data = roc_df,
    retained_predictors = attr(stats::terms(final_glm), "term.labels"),
    train_predictions = train_pred,
    test_predictions = test_pred,
    full_predictions = full_pred,
    train_map = train_map_sf,
    test_map = test_map_sf,
    geospatial_modeling = model_map_sf
  ))
}