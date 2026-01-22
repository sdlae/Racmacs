# ==============================================================================
# RACMACS BASELINE VALIDATION SCRIPT
# ==============================================================================
# Purpose: Compare your modified Racmacs baseline against the original package
#          to verify that stress calculations are consistent
#
# Instructions:
#   1. First, run this script with the ORIGINAL Racmacs package installed
#   2. Save the results
#   3. Then run with your MODIFIED Racmacs package
#   4. Compare the outputs
#
# ==============================================================================

# ==============================================================================
# CONFIGURATION - MODIFY THESE TO MATCH YOUR EXPERIMENT
# ==============================================================================

INPUT_TITER_FILE <- "shannen/h3map2004_hitable.csv"  # Path to your titer table
RANDOM_SEED <- 42                               # MUST match your experiment seed
N_DIM <- 2
N_OPTIMIZATIONS <- 500
MIN_COLUMN_BASIS <- "none"

# Output file for results (change name between original/modified runs)
OUTPUT_FILE <- "baseline_validation_results.rds"

# ==============================================================================
# SETUP
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("RACMACS BASELINE VALIDATION\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Load Racmacs
library(Racmacs)

# Print package info
cat("Package Information:\n")
cat("  Racmacs version:", as.character(packageVersion("Racmacs")), "\n")
cat("  Package location:", find.package("Racmacs"), "\n")

# Check if this is a modified version (look for custom options)
tryCatch({
  default_opts <- RacOptimizer.options()
  if ("optimize_colbases" %in% names(default_opts)) {
    cat("  Version type: MODIFIED (has optimize_colbases option)\n")
    is_modified <- TRUE
  } else {
    cat("  Version type: ORIGINAL (standard Racmacs)\n")
    is_modified <- FALSE
  }
}, error = function(e) {
  cat("  Version type: UNKNOWN (could not check options)\n")
  is_modified <- NA
})

cat("\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data from:", INPUT_TITER_FILE, "\n")

titer_table <- read.titerTable(INPUT_TITER_FILE)
map_base <- acmap(titer_table = titer_table)

cat("  Antigens:", numAntigens(map_base), "\n")
cat("  Sera:", numSera(map_base), "\n")
cat("  Measured titers:", sum(!is.na(titerTable(map_base)) & titerTable(map_base) != "*"), "\n")

# ==============================================================================
# RUN BASELINE OPTIMIZATION
# ==============================================================================

cat("\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat("Running baseline optimization...\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

set.seed(RANDOM_SEED)
cat("Random seed set to:", RANDOM_SEED, "\n\n")

start_time <- Sys.time()

# Run with explicit options to ensure consistency
if (is_modified) {
  # For modified version, explicitly disable all custom features
  cat("Using MODIFIED Racmacs - disabling all custom column base options\n")
  
  opts <- RacOptimizer.options(
    dim_annealing = FALSE,
    optimize_colbases = FALSE,
    bound_colbases = FALSE,
    regularize_colbases = FALSE
  )
  
  map <- optimizeMap(
    map_base,
    number_of_dimensions = N_DIM,
    number_of_optimizations = N_OPTIMIZATIONS,
    minimum_column_basis = MIN_COLUMN_BASIS,
    options = opts
  )
} else {
  # For original version, use defaults
  cat("Using ORIGINAL Racmacs - standard optimization\n")
  
  map <- optimizeMap(
    map_base,
    number_of_dimensions = N_DIM,
    number_of_optimizations = N_OPTIMIZATIONS,
    minimum_column_basis = MIN_COLUMN_BASIS
  )
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")

# ==============================================================================
# EXTRACT RESULTS
# ==============================================================================

cat("\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat("RESULTS\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

# Core metrics
stress <- mapStress(map)
colbases <- colBases(map)

cat("Optimization completed in", round(as.numeric(elapsed), 2), "minutes\n\n")

cat("STRESS:\n")
cat("  Total stress:", round(stress, 4), "\n")
cat("\n")

cat("COLUMN BASES:\n")
cat("  Min:", round(min(colbases), 4), "\n")
cat("  Max:", round(max(colbases), 4), "\n")
cat("  Mean:", round(mean(colbases), 4), "\n")
cat("  SD:", round(sd(colbases), 4), "\n")
cat("\n")

# Get coordinates
ag_coords <- agCoords(map)
sr_coords <- srCoords(map)

cat("COORDINATES:\n")
cat("  AG coords range X:", round(range(ag_coords[,1]), 4), "\n")
cat("  AG coords range Y:", round(range(ag_coords[,2]), 4), "\n")
cat("  SR coords range X:", round(range(sr_coords[,1]), 4), "\n")
cat("  SR coords range Y:", round(range(sr_coords[,2]), 4), "\n")
cat("\n")

# Get distances
table_dists <- tableDistances(map)
map_dists <- mapDistances(map)

# Calculate point-wise stress
valid_idx <- !is.na(table_dists) & !is.na(map_dists)
residuals <- (as.numeric(table_dists[valid_idx]) - as.numeric(map_dists[valid_idx]))

cat("RESIDUALS (Table Distance - Map Distance):\n")
cat("  Mean:", round(mean(residuals), 4), "\n")
cat("  SD:", round(sd(residuals), 4), "\n")
cat("  Min:", round(min(residuals), 4), "\n")
cat("  Max:", round(max(residuals), 4), "\n")
cat("  RMSD:", round(sqrt(mean(residuals^2)), 4), "\n")
cat("\n")

# ==============================================================================
# DETAILED COLUMN BASE OUTPUT
# ==============================================================================

cat(paste(rep("-", 80), collapse = ""), "\n")
cat("COLUMN BASES BY SERUM (first 20)\n")
cat(paste(rep("-", 80), collapse = ""), "\n\n")

cb_df <- data.frame(
  Serum = srNames(map),
  ColBase = round(colbases, 4)
)
print(head(cb_df, 20), row.names = FALSE)

# ==============================================================================
# SAVE RESULTS FOR COMPARISON
# ==============================================================================

results <- list(
  package_version = as.character(packageVersion("Racmacs")),
  package_location = find.package("Racmacs"),
  is_modified = is_modified,
  random_seed = RANDOM_SEED,
  n_optimizations = N_OPTIMIZATIONS,
  n_dimensions = N_DIM,
  min_column_basis = MIN_COLUMN_BASIS,
  stress = stress,
  colbases = colbases,
  colbase_stats = list(
    min = min(colbases),
    max = max(colbases),
    mean = mean(colbases),
    sd = sd(colbases)
  ),
  ag_coords = ag_coords,
  sr_coords = sr_coords,
  residuals = residuals,
  residual_stats = list(
    mean = mean(residuals),
    sd = sd(residuals),
    rmsd = sqrt(mean(residuals^2))
  ),
  elapsed_mins = as.numeric(elapsed),
  timestamp = Sys.time()
)

saveRDS(results, OUTPUT_FILE)
cat("\nResults saved to:", OUTPUT_FILE, "\n")

# ==============================================================================
# COMPARISON FUNCTION (run after you have both results)
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("TO COMPARE RESULTS:\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat('
# After running with BOTH original and modified packages, use this code:

original <- readRDS("baseline_validation_original.rds")
modified <- readRDS("baseline_validation_modified.rds")

cat("\\n=== COMPARISON ===\\n\\n")
cat("STRESS DIFFERENCE:\\n")
cat("  Original:", original$stress, "\\n")
cat("  Modified:", modified$stress, "\\n")
cat("  Difference:", modified$stress - original$stress, "\\n")
cat("  Percent diff:", 100 * (modified$stress - original$stress) / original$stress, "%\\n")

cat("\\nCOLUMN BASE COMPARISON:\\n")
cb_diff <- modified$colbases - original$colbases
cat("  Max absolute difference:", max(abs(cb_diff)), "\\n")
cat("  Mean absolute difference:", mean(abs(cb_diff)), "\\n")
cat("  All equal (tolerance 0.01):", all(abs(cb_diff) < 0.01), "\\n")

cat("\\nCOORDINATE COMPARISON:\\n")
# Note: Coordinates may differ due to rotation/reflection but distances should match
')

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("VALIDATION COMPLETE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")