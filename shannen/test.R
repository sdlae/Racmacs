source("shannen/helpers.R")
reload_cpp(".", compile_attrs = TRUE)

library(Racmacs)
options(RacOptimizer.num_cores = parallel::detectCores())
titer_table <- read.titerTable("shannen/h3map2004_hitable.csv")

cat("Number of antigens:", nrow(titer_table), "\n")
cat("Number of sera:", ncol(titer_table), "\n")

map_base <- acmap(titer_table = titer_table)

make_options <- function(optimize_colbases = FALSE) {
  opts <- RacOptimizer.options()
  opts$optimize_colbases <- optimize_colbases
  opts
}

n_dim <- 2

cat("\n========================================\n")
cat("Test 1: Fixed Column Bases (Original)\n")
cat("========================================\n")

set.seed(42)

map_fixed <- optimizeMap(
  map_base,
  number_of_dimensions = n_dim,
  number_of_optimizations = 500,
  minimum_column_basis = "none",
  options = make_options(optimize_colbases = FALSE)
)

stress_fixed <- mapStress(map_fixed)
colbases_fixed <- colBases(map_fixed)

cat("\nStress (fixed colbases):", stress_fixed, "\n")
cat("Column bases (first 10):\n")
print(head(colbases_fixed, 10))
cat("Column bases summary:\n")
print(summary(colbases_fixed))



cat("\n========================================\n")
cat("Test 2: Optimized Column Bases (New)\n")
cat("========================================\n")

set.seed(42)

map_optimized <- optimizeMap(
  map_base,
  number_of_dimensions = n_dim,
  number_of_optimizations = 500,
  minimum_column_basis = "none",
  options = make_options(optimize_colbases = TRUE)
)

stress_optimized <- mapStress(map_optimized)
colbases_optimized <- colBases(map_optimized)

cat("\nStress (optimized colbases):", stress_optimized, "\n")
cat("Column bases (first 10):\n")
print(head(colbases_optimized, 10))
cat("Column bases summary:\n")
print(summary(colbases_optimized))



cat("\n========================================\n")
cat("Comparison Summary\n")
cat("========================================\n")

stress_reduction <- stress_fixed - stress_optimized
stress_reduction_pct <- 100 * stress_reduction / stress_fixed

cat("\nStress (fixed colbases):    ", stress_fixed, "\n")
cat("Stress (optimized colbases):", stress_optimized, "\n")
cat("Stress reduction:           ", stress_reduction, "\n")
cat("Stress % reduction:         ", sprintf("%.2f%%", stress_reduction_pct), "\n")

cat("\nColumn base changes:\n")
colbase_diff <- colbases_optimized - colbases_fixed
cat("  Mean change:   ", mean(colbase_diff), "\n")
cat("  SD of change:  ", sd(colbase_diff), "\n")
cat("  Min change:    ", min(colbase_diff), "\n")
cat("  Max change:    ", max(colbase_diff), "\n")



cat("\n========================================\n")
cat("Detailed Column Base Comparison\n")
cat("========================================\n")

comparison_df <- data.frame(
  serum = srNames(map_base),
  fixed = colbases_fixed,
  optimized = colbases_optimized,
  difference = colbase_diff
)

comparison_df <- comparison_df[order(abs(comparison_df$difference), decreasing = TRUE), ]
cat("\nTop 10 sera with largest column base changes:\n")
print(head(comparison_df, 10))


results <- list(
  stress_fixed = stress_fixed,
  stress_optimized = stress_optimized,
  colbases_fixed = colbases_fixed,
  colbases_optimized = colbases_optimized,
  stress_reduction = stress_reduction,
  stress_reduction_pct = stress_reduction_pct,
  colbase_comparison = comparison_df
)

cat("\n")
cat("=== FINAL SUMMARY ===\n")
cat(sprintf("Dataset: %d antigens x %d sera\n", nrow(titer_table), ncol(titer_table)))
cat(sprintf("Stress (fixed):     %.4f\n", stress_fixed))
cat(sprintf("Stress (optimized): %.4f\n", stress_optimized))
cat(sprintf("Improvement:        %.4f (%.2f%%)\n", stress_reduction, stress_reduction_pct))

view(map_fixed)
view(map_optimized)