# Calgary Workflow Completion Note

Completed in this run:
- Created the Calgary-only Quarto workflow at `reports/calgary_flood_workflow.qmd`.
- Implemented helper R scripts for fishnet and target generation, feature engineering, merge and EDA, and glm modeling and diagnostics.
- Generated the Calgary fishnet and binary target outputs in `outputs/intermediate/`.

Explicit assumptions used:
- `data/` is the only authoritative input folder for this run.
- `inundation_calgary.tif` uses `0 = no inundation`, `1 = inundation`, and other values are excluded from the denominator.
- `calgary_surf_water.tif` is binary and can be used directly for `water_cover_area`.
- Conservative land-cover mapping is used for vegetation, open soil, agriculture, and impervious classes.

Parts that may still need user confirmation or runtime follow-up:
- The watershed step is computationally heavy and should be monitored closely on the first full end-to-end run.
- Large polygon-overlap steps may need performance tuning if they are run repeatedly at full Calgary scale.

Next step for City X deployment:
- Reproduce the exact same predictor names and preprocessing contract for City X, then apply the saved Calgary scaling rules and glm model to the City X feature table.
