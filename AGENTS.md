You are the lead coding agent for a flood inundation probability modeling project in R.

Your goal in this run is to build a fully reproducible Calgary-only workflow that:
1. creates the fishnet,
2. generates the binary flood target,
3. computes all listed features,
4. performs EDA and minor data cleaning,
5. builds a glm logistic regression model,
6. performs validation and diagnostics,
7. writes everything into a single well-structured Quarto file (.qmd).

Important scope limit:
- Do NOT implement City X prediction in this run.
- Stop after the Calgary model, validation, and diagnostics are complete.
- However, keep the code structure reusable so that the same feature-engineering pipeline can later be applied to City X with the exact same variable names.

Project constraints:
- All code must be in R.
- Final deliverable of this run is a .qmd file.
- The QMD must contain:
  - markdown section headings,
  - code chunks for every step,
  - English comments inside code chunks explaining what each step does,
  - short markdown discussion paragraphs interpreting outputs.
- Make the workflow reproducible and transparent.
- Use glm logistic regression only.
- Use ggplot for all plots and maps in the report.
- Use local projected CRS appropriate for Calgary and preserve consistent CRS across layers.
- Use subagents to split major sections, but keep one master agent orchestrating dependencies between sections.

Working style:
- First inspect the available files and data structure.
- If file names, land-cover class names, CRS, raster encoding, or geometry validity create ambiguity, pause and ask the user a concise clarifying question instead of guessing.
- If one implementation path fails, debug and try a second reasonable approach.
- Keep all assumptions explicit in comments.
- Save intermediate outputs when useful so later sections can reuse them.

Suggested subagent structure:
- Subagent A: environment setup, file inspection, CRS handling, fishnet creation, target variable generation.
- Subagent B: topography and hydrology feature engineering.
- Subagent C: land-cover and built-environment feature engineering.
- Subagent D: merge features, EDA, data cleaning, and exploratory plots/maps.
- Subagent E: logistic regression, validation, diagnostics, and report assembly.
- Master agent: manage handoffs, ensure consistent variable names, ensure all outputs are joined by fishnet id, and assemble the final QMD.

Expected final objects:
- A fishnet sf object or equivalent geodataframe with one row per 30 m grid cell.
- Required columns:
  - id
  - target variable (binary inundation)
  - all engineered features
- A cleaned modeling dataset derived from that fishnet.
- A final fitted glm object.
- Validation outputs:
  - model summary table,
  - ROC curve,
  - confusion metrics table,
  - simple explanatory table for TP, TN, FP, FN.

Data and target construction:
1. Create fishnet from the Calgary city boundary shapefile.
   - Grid size = 30 m.
   - Use the same projected CRS as the boundary, or reproject everything to an appropriate local projected CRS for Calgary before analysis.
   - Assign a unique integer id to each grid cell.

2. Generate the binary target variable from inundation.tif.
   - inundation.tif is 10 m resolution.
   - For each 30 m fishnet cell, calculate the share of inundation raster cells inside that grid.
   - Valid raster values:
     - 0 = no inundation
     - 1 = inundation
     - other values = nodata / exclude from denominator
   - Define target:
     - if inundated share >= 0.30, target = 1
     - if inundated share < 0.30, target = 0
   - Remove or flag fishnet cells with no valid inundation raster coverage.
   - Be careful with cells north of raster coverage; do not keep invalid target rows.

Feature list to compute for each fishnet cell:
Topography
- mean_elev
- min_elev
- elev_range
- mean_slope
- max_slope
- sd_elev

Hydrology
- dist_nearest_stream
- water_cover_area
- river_density

Watershed
- max_flow_accum
- upstream_ws_area

Land cover
- impervious_ratio
- vegetation_ratio
- open_soil_ratio
- agriculture_ratio

Built environment
- building_cover_ratio

3.1 Detailed feature calculation rules
Implement the following exactly and document the method in the QMD.

A. Mean elevation
- Source: calgary_dem.tif
- Method: zonal statistics over each fishnet cell using DEM pixels intersecting the grid.
- Output: mean of DEM values in meters for each grid cell.

B. Min. elevation
- Source: calgary_dem.tif
- Method: zonal minimum over each fishnet cell.
- Output: minimum DEM value in meters.

C. Elevation range
- Source: calgary_dem.tif
- Method: zonal maximum minus zonal minimum within each fishnet cell.
- Output: max_elev - min_elev.

D. Mean slope
- Source: derived from calgary_dem.tif
- Method:
  1. derive slope raster from DEM in degrees or percent,
  2. compute zonal mean slope for each fishnet cell.
- Output: average slope within the grid.
- Use the same unit consistently for all slope outputs and clearly state it in the QMD.

E. Max slope
- Source: slope raster derived from calgary_dem.tif
- Method: zonal maximum slope within each fishnet cell.

F. Std. of elevation
- Source: calgary_dem.tif
- Method: zonal standard deviation of DEM values within each fishnet cell.

G. Dist. to nearest stream
- Source: calgary_waterways.shp
- Method:
  1. ensure stream layer is in same projected CRS as fishnet,
  2. compute centroid of each fishnet cell,
  3. calculate Euclidean distance from centroid to nearest stream geometry.
- Output: distance in meters.
- If needed, use sf::st_distance with nearest feature logic.

H. Water cover area
- Source: calgary_surf_water.tif
- Method:
  1. interpret raster so water cells are distinguished from non-water / nodata,
  2. within each fishnet cell, sum the area of raster cells classified as surface water.
- Output: water-covered area in square meters per grid.
- If the raster is binary, multiply count of water cells by cell area.
- If classification codes are ambiguous, inspect raster values first and ask user if needed.

I. River Density
- Source: calgary_waterways.shp
- Method:
  1. intersect stream lines with each fishnet polygon,
  2. sum total stream length inside each grid cell,
  3. divide by fishnet cell area.
- Output: stream length per square meter (or optionally per hectare, but keep unit explicit and consistent).
- Prefer a simple length/area definition.

J. Max flow accumulation
- Source: calgary_dem.tif
- Method:
  1. hydrologically preprocess DEM: fill sinks/depressions,
  2. derive flow direction,
  3. derive flow accumulation raster,
  4. compute zonal maximum flow accumulation value within each fishnet cell.
- Output: maximum flow accumulation value in the grid.
- Use an R-based hydrology workflow. Suitable packages may include terra plus whitebox, RSAGA, or another robust R-accessible hydrology backend.
- Document the exact algorithm and units/count interpretation.

K. Upstream watershed area of fishnet centroid point
- Source: calgary_dem.tif
- This feature must follow the pour-point interpretation:
  1. compute centroid for each fishnet cell,
  2. hydrologically preprocess DEM,
  3. derive flow direction and flow accumulation,
  4. optionally derive stream network or use high-flow cells as candidate drainage paths,
  5. snap each centroid to the nearest appropriate drainage cell / stream cell,
  6. delineate the upstream watershed for that snapped point,
  7. calculate watershed area.
- Output: upstream contributing area in square meters for the fishnet centroid point.
- This is NOT simply local cell area and NOT just raw flow accumulation unless flow accumulation is explicitly converted to area.
- Entire workflow must be implemented in R.
- If delineating a watershed for every centroid is computationally heavy, implement a robust batch workflow and cache intermediate rasters/results.
- Document clearly that this is “upstream watershed area of fishnet centroid point” based on a pour-point approach.

L. Impervious surface ratio
- Source: calgary_Imp_surf.geojson
- Method:
  1. intersect impervious polygons with fishnet,
  2. sum impervious area within each grid cell,
  3. divide by fishnet cell area.
- Output: proportion from 0 to 1.

M. Vegetation cover proportion
- Source: calgary_landcover.geojson
- Method:
  1. identify vegetation-related class(es) from the land-cover dataset,
  2. intersect with fishnet,
  3. sum vegetation area in each grid,
  4. divide by fishnet cell area.
- Output: proportion from 0 to 1.
- Inspect actual class names first. Do not guess silently.

N. Open soil land proportion
- Source: calgary_landcover.geojson
- Method:
  1. identify open soil / bare soil / barren-related class(es),
  2. intersect with fishnet,
  3. sum area,
  4. divide by fishnet cell area.
- Output: proportion from 0 to 1.

O. Agricultural land proportion
- Source: calgary_landcover.geojson
- Method:
  1. identify agriculture-related class(es),
  2. intersect with fishnet,
  3. sum area,
  4. divide by fishnet cell area.
- Output: proportion from 0 to 1.

P. Building cover ratio
- Source: calgary_building.shp
- Method:
  1. intersect building footprints with fishnet,
  2. sum building footprint area inside each grid,
  3. divide by fishnet cell area.
- Output: proportion from 0 to 1.

Implementation guidance for feature engineering:
- Prefer terra for raster handling and sf for vector handling.
- Use exactextractr when helpful for efficient raster summaries over polygons.
- Use tidyverse-style data management when it improves readability.
- Before computing ratios, ensure all geometries are valid and CRS are aligned.
- For area-based ratios, use projected CRS in meters.
- For any feature with missing values introduced by absent intersections, replace with 0 only when that is substantively correct (for example, no building intersection = 0 building cover ratio). Otherwise keep NA and discuss.
- Keep variable names compact, consistent, and ready for later City X transfer.

EDA and data cleaning:
4. Perform exploratory data analysis and minor cleaning.
Required outputs:
4.1 Histograms and spatial distribution maps for target and features.
- Use ggplot for all plots and maps.
- Produce a reasonable subset in the QMD if the full set is too long, but save all figures if practical.

4.2 Multicollinearity checks
- Correlation matrix for numeric predictors.
- VIF test for candidate logistic model predictors.
- Discuss whether any variables should be removed or retained despite correlation, especially watershed variables.

4.3 Minor, realistic cleaning only
- Do not over-clean.
- Apply practical steps that improve robustness and later transferability:
  - check missingness,
  - check zero-variance or near-zero-variance predictors,
  - inspect extreme outliers,
  - confirm binary target balance,
  - remove invalid rows created by missing target coverage,
  - consider simple transformations only if justified and documented.
- Keep the dataset usable for later prediction in another city.

Modeling:
5. Build the glm logistic regression model.
- Standardize all predictor variables with z-score normalization before modeling.
- Keep the target as binary 0/1.
- Create a training/testing split.
- Use a reproducible random seed.
- Try a rational candidate model selection process, but remain within glm logistic regression.
- Preserve final chosen feature names clearly.
- Save both:
  - the preprocessed modeling table,
  - the fitted model object.

Validation and diagnostics:
6. Produce validation and diagnostic outputs.
Must include:
- final logistic regression model summary table,
- ROC curve plot,
- table of confusion metrics including at least:
  - True Positives,
  - True Negatives,
  - False Positives,
  - False Negatives,
  - Accuracy,
  - Precision,
  - Recall,
  - F1 if straightforward,
  - Specificity if straightforward.
- a simple explanatory table with plain-language definitions of TP, TN, FP, FN.
- short markdown discussion of:
  - modeling workflow,
  - train/test process,
  - threshold used for classification,
  - model strengths,
  - limitations,
  - how reliable the model appears for future transfer.

Recommended outputs for report readiness:
- A map of observed target in Calgary.
- A map of selected feature examples.
- A map of model prediction probabilities for Calgary.
- A map of classification outcome categories on the test or evaluation set if practical.
These are useful even if the full final course report will later expand them.

QMD structure to generate:
- Title
- Setup / libraries / reproducibility notes
- Data inputs
- Fishnet creation
- Target variable generation
- Feature engineering
  - topography
  - hydrology
  - watershed
  - land cover
  - built environment
- Merge final fishnet feature table
- EDA
- Data cleaning
- Logistic regression modeling
- Validation and diagnostics
- Key findings / discussion
- Limitations and next steps

Code quality requirements:
- Every code chunk should include English comments.
- Keep functions modular where possible.
- Use relative paths where possible.
- Avoid hard-coding values unless documented.
- Print informative checks after major steps:
  - CRS checks,
  - row counts,
  - number of missing values,
  - summary statistics,
  - class balance,
  - successful joins.

Error-handling expectations:
- If a dataset name in the folder differs slightly from the expected one, inspect and adapt.
- If land-cover class names do not clearly match vegetation/open soil/agriculture, stop and ask the user to confirm mapping.
- If hydrology tools for watershed delineation fail in one package, try an alternative R-accessible method.
- If the watershed-centroid feature is too slow, optimize rather than dropping it.
- Do not silently omit required features.

Final deliverables of this run:
1. One completed .qmd file.
2. Any helper R scripts only if they improve modularity.
3. A final fishnet feature dataset saved to disk.
4. A brief completion note stating:
   - what was completed,
   - what assumptions were made,
   - any parts needing user confirmation,
   - what should be done next for City X deployment.