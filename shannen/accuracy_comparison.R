source("shannen/cpp_compile.R")
reload_cpp(".", compile_attrs = TRUE)

library(Racmacs)
options(RacOptimizer.num_cores = parallel::detectCores())

# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_DIM <- 5                    # Number of dimensions
N_OPTIMIZATIONS <- 500        # Number of optimization runs
BLIND_FRACTION <- 0.15        # Fraction of data held out for testing
MIN_COLUMN_BASIS <- "none"    # Minimum column basis
RANDOM_SEED <- 42             # For reproducibility

source("shannen/accuracy_helpers.R")

titer_table <- read.titerTable("shannen/h3map2004_hitable.csv")
cat(sprintf("Dataset: %d antigens x %d sera\n", nrow(titer_table), ncol(titer_table)))
cat(sprintf("Dimensions: %d\n", N_DIM))
cat(sprintf("Blind fraction: %.0f%%\n", BLIND_FRACTION * 100))

set.seed(RANDOM_SEED)
blind_test <- prepare_blind_test(
  as.matrix(titer_table), 
  blind_fraction = BLIND_FRACTION
)

# Create base map with training data only
map_base <- acmap(titer_table = blind_test$train)

# ==============================================================================
# FIXED COLUMN BASES
# ==============================================================================

set.seed(RANDOM_SEED)

map_fixed <- optimizeMap(
  map_base,
  number_of_dimensions = N_DIM,
  number_of_optimizations = N_OPTIMIZATIONS,
  minimum_column_basis = MIN_COLUMN_BASIS,
  options = make_options(optimize_colbases = FALSE)
)

stress_fixed <- mapStress(map_fixed)
colbases_fixed <- colBases(map_fixed)

cat(sprintf("\nStress: %.4f\n", stress_fixed))
cat("Column bases summary:\n")
print(summary(colbases_fixed))

cat("\nPredicting blind titers...\n")
pred_fixed <- predict_titers(map_fixed, blind_test$blind)
stats_fixed <- calculate_error_stats(pred_fixed, blind_test$blind)

cat(sprintf("Prediction accuracy (fixed colbases):\n"))
cat(sprintf("  Valid predictions: %d/%d\n", stats_fixed$n_valid, nrow(blind_test$blind)))
cat(sprintf("  Mean error:   %.4f log2 units\n", stats_fixed$mean_error))
cat(sprintf("  Median error: %.4f log2 units\n", stats_fixed$median_error))
cat(sprintf("  RMSE:         %.4f log2 units\n", stats_fixed$rmse))
cat(sprintf("  Correlation:  %.4f\n", stats_fixed$correlation))


# ==============================================================================
# OPTIMIZED COLUMN BASES
# ==============================================================================


set.seed(RANDOM_SEED)

map_optimized <- optimizeMap(
  map_base,
  number_of_dimensions = N_DIM,
  number_of_optimizations = N_OPTIMIZATIONS,
  minimum_column_basis = MIN_COLUMN_BASIS,
  options = make_options(optimize_colbases = TRUE)
)

stress_optimized <- mapStress(map_optimized)
colbases_optimized <- colBases(map_optimized)

cat(sprintf("\nStress: %.4f\n", stress_optimized))
cat("Column bases summary:\n")
print(summary(colbases_optimized))

pred_optimized <- predict_titers(map_optimized, blind_test$blind)
stats_optimized <- calculate_error_stats(pred_optimized, blind_test$blind)

cat(sprintf("Prediction accuracy (optimized colbases):\n"))
cat(sprintf("  Valid predictions: %d/%d\n", stats_optimized$n_valid, nrow(blind_test$blind)))
cat(sprintf("  Mean error:   %.4f log2 units\n", stats_optimized$mean_error))
cat(sprintf("  Median error: %.4f log2 units\n", stats_optimized$median_error))
cat(sprintf("  RMSE:         %.4f log2 units\n", stats_optimized$rmse))
cat(sprintf("  Correlation:  %.4f\n", stats_optimized$correlation))



# ==============================================================================
# RESULTS
# ==============================================================================


stress_reduction <- stress_fixed - stress_optimized
stress_reduction_pct <- 100 * stress_reduction / stress_fixed

cat(sprintf("\n--- STRESS ---\n"))
cat(sprintf("Fixed colbases:     %.4f\n", stress_fixed))
cat(sprintf("Optimized colbases: %.4f\n", stress_optimized))
cat(sprintf("Reduction:          %.4f (%.2f%%)\n", stress_reduction, stress_reduction_pct))

error_reduction <- stats_fixed$mean_error - stats_optimized$mean_error
error_reduction_pct <- 100 * error_reduction / stats_fixed$mean_error

cat(sprintf("\n--- PREDICTION ACCURACY ---\n"))
cat(sprintf("                        Fixed      Optimized    Change\n"))
cat(sprintf("Mean error (log2):      %.4f     %.4f       %+.4f (%.1f%%)\n", 
            stats_fixed$mean_error, stats_optimized$mean_error, 
            -error_reduction, -error_reduction_pct))
cat(sprintf("Median error (log2):    %.4f     %.4f       %+.4f\n", 
            stats_fixed$median_error, stats_optimized$median_error,
            stats_optimized$median_error - stats_fixed$median_error))
cat(sprintf("RMSE (log2):            %.4f     %.4f       %+.4f\n", 
            stats_fixed$rmse, stats_optimized$rmse,
            stats_optimized$rmse - stats_fixed$rmse))
cat(sprintf("Correlation:            %.4f     %.4f       %+.4f\n", 
            stats_fixed$correlation, stats_optimized$correlation,
            stats_optimized$correlation - stats_fixed$correlation))

colbase_diff <- colbases_optimized - colbases_fixed
cat(sprintf("\n--- COLUMN BASE CHANGES ---\n"))
cat(sprintf("Mean change:   %+.4f\n", mean(colbase_diff)))
cat(sprintf("SD of change:  %.4f\n", sd(colbase_diff)))
cat(sprintf("Min change:    %+.4f\n", min(colbase_diff)))
cat(sprintf("Max change:    %+.4f\n", max(colbase_diff)))

# Interpretation
cat(sprintf("\n--- INTERPRETATION ---\n"))
if (stress_optimized < stress_fixed) {
  cat("✓ Optimized colbases achieved LOWER stress (better fit to training data)\n")
} else {
  cat("✗ Optimized colbases achieved HIGHER stress (worse fit to training data)\n")
}

if (!is.na(stats_optimized$mean_error) && !is.na(stats_fixed$mean_error)) {
  if (stats_optimized$mean_error < stats_fixed$mean_error) {
    cat("✓ Optimized colbases achieved LOWER prediction error (better generalization)\n")
  } else if (stats_optimized$mean_error > stats_fixed$mean_error) {
    cat("✗ Optimized colbases achieved HIGHER prediction error (possible overfitting)\n")
  } else {
    cat("= Optimized colbases achieved SAME prediction error\n")
  }
}

# ==============================================================================
# DETAILED COLUMN BASE COMPARISON
# ==============================================================================

cat("\n========================================\n")
cat("Top 10 Sera with Largest Column Base Changes\n")
cat("========================================\n")

comparison_df <- data.frame(
  serum = srNames(map_base),
  fixed = colbases_fixed,
  optimized = colbases_optimized,
  difference = colbase_diff
)
comparison_df <- comparison_df[order(abs(comparison_df$difference), decreasing = TRUE), ]
print(head(comparison_df, 10))

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

results <- list(
  config = list(
    n_dim = N_DIM,
    n_optimizations = N_OPTIMIZATIONS,
    blind_fraction = BLIND_FRACTION,
    random_seed = RANDOM_SEED
  ),
  fixed = list(
    map = map_fixed,
    stress = stress_fixed,
    colbases = colbases_fixed,
    predictions = pred_fixed,
    stats = stats_fixed
  ),
  optimized = list(
    map = map_optimized,
    stress = stress_optimized,
    colbases = colbases_optimized,
    predictions = pred_optimized,
    stats = stats_optimized
  ),
  comparison = list(
    stress_reduction = stress_reduction,
    stress_reduction_pct = stress_reduction_pct,
    error_reduction = error_reduction,
    error_reduction_pct = error_reduction_pct,
    colbase_changes = comparison_df
  ),
  blind_test = blind_test
)

cat("\n========================================\n")
cat("Results saved in 'results' variable\n")
cat("========================================\n")

# Final summary table
cat("\n")
cat("=== FINAL SUMMARY TABLE ===\n")
cat(sprintf("Dimensions: %d\n", N_DIM))
cat(sprintf("%-25s %12s %12s %12s\n", "Metric", "Fixed", "Optimized", "Change"))
cat(sprintf("%-25s %12.2f %12.2f %12.2f%%\n", "Stress", stress_fixed, stress_optimized, -stress_reduction_pct))
cat(sprintf("%-25s %12.4f %12.4f %12.2f%%\n", "Mean Error (log2)", stats_fixed$mean_error, stats_optimized$mean_error, -error_reduction_pct))
cat(sprintf("%-25s %12.4f %12.4f\n", "RMSE (log2)", stats_fixed$rmse, stats_optimized$rmse))
cat(sprintf("%-25s %12.4f %12.4f\n", "Correlation", stats_fixed$correlation, stats_optimized$correlation))
cat(sprintf("%-25s %12.2f %12.2f\n", "Mean ColBase", mean(colbases_fixed), mean(colbases_optimized)))