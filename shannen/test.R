# shannen/test.R
# Test script for comparing maps with fixed vs optimized column bases

source("shannen/helpers.R")

# Reload the C++ code with new changes
reload_cpp(".", compile_attrs = TRUE)

# ============================================
# Create a simple test titer table
# ============================================

# Example: Create a small titer table for testing
# Replace this with your actual data
titer_matrix <- matrix(
  c(
    "1280", "640", "320", "160",
    "640", "1280", "640", "320",
    "320", "640", "1280", "640",
    "160", "320", "640", "1280"
  ),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(
    c("AG1", "AG2", "AG3", "AG4"),
    c("SR1", "SR2", "SR3", "SR4")
  )
)

# Create an acmap from the titer table
map_base <- acmap(
  titer_table = titer_matrix
)

# ============================================
# Define custom optimizer options function 
# (if not already in the package)
# ============================================

# This function creates the options list that gets passed to C++
# You may need to modify your existing RacOptimizer.options function
# or use this one
my_optimizer_options <- function(
    dim_annealing = FALSE,
    method = "L-BFGS",
    maxit = 1000,
    num_basis = 10,
    armijo_constant = 1e-4,
    wolfe = 0.9,
    min_gradient_norm = 1e-6,
    factr = 1e7,
    max_line_search_trials = 50,
    min_step = 1e-20,
    max_step = 1e20,
    num_cores = 1,
    report_progress = TRUE,
    progress_bar_length = 50,
    optimize_colbases = FALSE  # NEW PARAMETER
) {
  list(
    dim_annealing = dim_annealing,
    method = method,
    maxit = maxit,
    num_basis = num_basis,
    armijo_constant = armijo_constant,
    wolfe = wolfe,
    min_gradient_norm = min_gradient_norm,
    factr = factr,
    max_line_search_trials = max_line_search_trials,
    min_step = min_step,
    max_step = max_step,
    num_cores = num_cores,
    report_progress = report_progress,
    progress_bar_length = progress_bar_length,
    optimize_colbases = optimize_colbases
  )
}

# ============================================
# Test 1: Fixed column bases (original behavior)
# ============================================

cat("\n========================================\n")
cat("Test 1: Fixed Column Bases (Original)\n")
cat("========================================\n")

set.seed(42)

map_fixed <- optimizeMap(
  map_base,
  number_of_dimensions = 2,
  number_of_optimizations = 100,
  minimum_column_basis = "none",
  options = my_optimizer_options(
    optimize_colbases = FALSE  # Original behavior
  )
)

stress_fixed <- mapStress(map_fixed)
colbases_fixed <- colBases(map_fixed)

cat("\nStress (fixed colbases):", stress_fixed, "\n")
cat("Column bases:\n")
print(colbases_fixed)

# ============================================
# Test 2: Optimized column bases (new feature)
# ============================================

cat("\n========================================\n")
cat("Test 2: Optimized Column Bases (New)\n")
cat("========================================\n")

set.seed(42)  # Same seed for fair comparison

map_optimized <- optimizeMap(
  map_base,
  number_of_dimensions = 2,
  number_of_optimizations = 100,
  minimum_column_basis = "none",
  options = my_optimizer_options(
    optimize_colbases = TRUE  # New behavior
  )
)

stress_optimized <- mapStress(map_optimized)
colbases_optimized <- colBases(map_optimized)

cat("\nStress (optimized colbases):", stress_optimized, "\n")
cat("Column bases:\n")
print(colbases_optimized)

# ============================================
# Comparison
# ============================================

cat("\n========================================\n")
cat("Comparison\n")
cat("========================================\n")

cat("\nStress reduction:", stress_fixed - stress_optimized, "\n")
if (stress_fixed > 0) {
  cat("Stress % reduction:", 100 * (stress_fixed - stress_optimized) / stress_fixed, "%\n")
}

cat("\nColumn base differences (optimized - fixed):\n")
print(colbases_optimized - colbases_fixed)

# Check that stress is reduced (or equal)
if (stress_optimized <= stress_fixed) {
  cat("\n✓ SUCCESS: Optimized stress <= Fixed stress\n")
} else {
  cat("\n✗ WARNING: Optimized stress > Fixed stress (unexpected)\n")
}

# ============================================
# Visualize (if viewer available)
# ============================================

cat("\n========================================\n")
cat("Visualization\n")
cat("========================================\n")

# Try to plot if the viewer is available
tryCatch({
  par(mfrow = c(1, 2))
  plot(map_fixed, main = "Fixed Column Bases")
  plot(map_optimized, main = "Optimized Column Bases")
}, error = function(e) {
  cat("Could not create plot:", e$message, "\n")
})

# ============================================
# Export results
# ============================================

results <- list(
  stress_fixed = stress_fixed,
  stress_optimized = stress_optimized,
  colbases_fixed = colbases_fixed,
  colbases_optimized = colbases_optimized,
  stress_reduction = stress_fixed - stress_optimized
)

cat("\nResults saved in 'results' variable\n")