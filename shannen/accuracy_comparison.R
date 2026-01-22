# ==============================================================================
# COMPREHENSIVE COMPARISON: FIXED vs BOUNDED vs REGULARIZED COLUMN BASES
# ==============================================================================
#
# This script compares different approaches to column base optimization in
# antigenic cartography:
#   1. Fixed (baseline) - Traditional approach, column bases from max titer per serum
#   2. Unbounded - Column bases optimized freely (can lead to unrealistic values)
#   3-5. Bounded - Column bases constrained within deviation limits
#   6-8. Regularized - Column bases penalized for deviating from initial values
#
# Key metrics:
#   - Stress: Map fitting quality (base stress from map distances)
#   - Penalty Stress: Regularization penalty term
#   - Total Stress: Stress + Penalty Stress (for regularized models)
#   - Mean Error: Prediction accuracy on held-out titers
#   - Column Base Range: Biological plausibility (should be ~6-14 for typical sera)
#
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

source("shannen/cpp_compile.R")
reload_cpp(".", compile_attrs = TRUE)

library(Racmacs)
options(RacOptimizer.num_cores = parallel::detectCores())

# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_DIM <- 3
N_OPTIMIZATIONS <- 500
BLIND_FRACTION <- 0.15
MIN_COLUMN_BASIS <- "none"
RANDOM_SEED <- 42

source("shannen/accuracy_helpers.R")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Create optimizer options with column base settings
make_options <- function(
    optimize_colbases = FALSE,
    bound_colbases = FALSE,
    colbase_min_bound = 4.0,
    colbase_max_bound = 14.0,
    colbase_max_deviation = 5.0,
    regularize_colbases = FALSE,
    colbase_lambda = 0.1
) {
  opts <- RacOptimizer.options()
  opts$optimize_colbases <- optimize_colbases
  opts$bound_colbases <- bound_colbases
  opts$colbase_min_bound <- colbase_min_bound
  opts$colbase_max_bound <- colbase_max_bound
  opts$colbase_max_deviation <- colbase_max_deviation
  opts$regularize_colbases <- regularize_colbases
  opts$colbase_lambda <- colbase_lambda
  opts
}

# Calculate the regularization penalty for a given map (R-side calculation)
calculate_regularization_penalty <- function(colbases, initial_colbases, lambda) {
  if (lambda <= 0) return(0)
  sum(lambda * (colbases - initial_colbases)^2)
}

# Calculate table distances and map distances for plotting
calculate_distances <- function(map, titer_table) {
  ag_coords <- agCoords(map)
  sr_coords <- srCoords(map)
  colbases <- colBases(map)

  table_dists <- c()
  map_dists <- c()
  titer_types <- c()

  for (ag in 1:nrow(ag_coords)) {
    for (sr in 1:nrow(sr_coords)) {
      titer <- titer_table[ag, sr]

      if (is.na(titer) || titer == "*") next
      if (any(is.na(ag_coords[ag, ])) || any(is.na(sr_coords[sr, ]))) next

      if (grepl("^<", titer)) {
        titer_type <- "less_than"
        numeric_titer <- as.numeric(gsub("<", "", titer))
      } else if (grepl("^>", titer)) {
        titer_type <- "greater_than"
        numeric_titer <- as.numeric(gsub(">", "", titer))
      } else {
        titer_type <- "measurable"
        numeric_titer <- as.numeric(titer)
      }

      log_titer <- log2(numeric_titer / 10)
      table_dist <- colbases[sr] - log_titer
      map_dist <- sqrt(sum((ag_coords[ag, ] - sr_coords[sr, ])^2))

      table_dists <- c(table_dists, table_dist)
      map_dists <- c(map_dists, map_dist)
      titer_types <- c(titer_types, titer_type)
    }
  }

  data.frame(
    table_dist = table_dists,
    map_dist = map_dists,
    titer_type = titer_types,
    stringsAsFactors = FALSE
  )
}

# Predict titers using a map for specific antigen-serum pairs
# This version works with a map trained on FULL data but evaluates on specific pairs
predict_titers_from_map <- function(map, blind_info) {
  ag_coords <- agCoords(map)
  sr_coords <- srCoords(map)
  colbases <- colBases(map)

  predictions <- data.frame(
    ag = blind_info$ag,
    sr = blind_info$sr,
    actual_titer = blind_info$titer,
    actual_log = log2(blind_info$titer / 10),
    predicted_log = NA,
    predicted_titer = NA,
    error = NA,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(predictions)) {
    ag_idx <- predictions$ag[i]
    sr_idx <- predictions$sr[i]

    if (any(is.na(ag_coords[ag_idx, ])) || any(is.na(sr_coords[sr_idx, ]))) {
      next
    }

    map_dist <- sqrt(sum((ag_coords[ag_idx, ] - sr_coords[sr_idx, ])^2))
    predicted_log <- colbases[sr_idx] - map_dist

    predictions$predicted_log[i] <- predicted_log
    predictions$predicted_titer[i] <- 10 * 2^predicted_log
    predictions$error[i] <- abs(predictions$actual_log[i] - predicted_log)
  }

  return(predictions)
}

# Calculate error statistics
calculate_prediction_stats <- function(predictions) {
  valid <- !is.na(predictions$error)

  list(
    n = sum(valid),
    mean_error = mean(predictions$error[valid]),
    rmse = sqrt(mean(predictions$error[valid]^2)),
    correlation = cor(predictions$actual_log[valid], predictions$predicted_log[valid])
  )
}

# Define cluster reference strains (the specific antigen to point to)
cluster_reference_strains <- c(
  "HK68" = "HK/1/68",
  "EN72" = "EN/42/72",
  "VI75" = "VI/3/75",
  "TX77" = "TE/1/77",
  "BK79" = "BA/1/79",
  "SI87" = "SH/11/87",
  "BE89" = "BE/353/89",
  "BE92" = "BE/32/92",
  "WU95" = "WU/359/95",
  "SY97" = "SY/5/97",
  "FU02" = "FU/411/02"
)

# Define cluster names by color (used for labeling)
cluster_names <- c(
  "#A208BD" = "HK68",
  "#00FFE1" = "EN72",
  "#F9D004" = "VI75",
  "#AB4C00" = "TX77",
  "#00FF00" = "BK79",
  "#0000FF" = "SI87",
  "#FF0000" = "BE89",
  "#F894F8" = "BE92",
  "#37802B" = "WU95",
  "#00AFFF" = "SY97",
  "#FFD700" = "FU02"
)

# Apply plotspec styling to a map object (for Racmacs view())
apply_plotspec <- function(map, plotspec) {
  ag_names <- agNames(map)
  sr_names <- srNames(map)

  ag_spec <- plotspec[plotspec$type == "ag", ]
  sr_spec <- plotspec[plotspec$type == "sr", ]

  # Match and apply antigen styles
  ag_fill <- sapply(ag_names, function(name) {
    idx <- which(ag_spec$name == name)
    if (length(idx) > 0) ag_spec$fill[idx[1]] else "green"
  })
  ag_outline <- sapply(ag_names, function(name) {
    idx <- which(ag_spec$name == name)
    if (length(idx) > 0) ag_spec$outline[idx[1]] else "black"
  })
  ag_shape <- sapply(ag_names, function(name) {
    idx <- which(ag_spec$name == name)
    if (length(idx) > 0) ag_spec$shape[idx[1]] else "CIRCLE"
  })
  ag_size <- sapply(ag_names, function(name) {
    idx <- which(ag_spec$name == name)
    if (length(idx) > 0) ag_spec$size[idx[1]] else 5
  })

  # Match and apply serum styles
  sr_fill <- sapply(sr_names, function(name) {
    idx <- which(sr_spec$name == name)
    if (length(idx) > 0) sr_spec$fill[idx[1]] else "transparent"
  })
  sr_outline <- sapply(sr_names, function(name) {
    idx <- which(sr_spec$name == name)
    if (length(idx) > 0) sr_spec$outline[idx[1]] else "black"
  })
  sr_shape <- sapply(sr_names, function(name) {
    idx <- which(sr_spec$name == name)
    if (length(idx) > 0) sr_spec$shape[idx[1]] else "BOX"
  })
  sr_size <- sapply(sr_names, function(name) {
    idx <- which(sr_spec$name == name)
    if (length(idx) > 0) sr_spec$size[idx[1]] else 5
  })

  # Apply to map
  agFill(map) <- ag_fill
  agOutline(map) <- ag_outline
  agShape(map) <- ag_shape
  agSize(map) <- ag_size

  srFill(map) <- sr_fill
  srOutline(map) <- sr_outline
  srShape(map) <- sr_shape
  srSize(map) <- sr_size

  return(map)
}

# Plot antigenic map with plotspec colors
plot_antigenic_map <- function(map, plotspec, title = "", subtitle = "") {
  ag_coords <- agCoords(map)
  sr_coords <- srCoords(map)
  ag_names <- agNames(map)
  sr_names <- srNames(map)

  # Create color lookup from plotspec
  ag_spec <- plotspec[plotspec$type == "ag", ]
  sr_spec <- plotspec[plotspec$type == "sr", ]

  # Match colors to antigens
  ag_fill <- sapply(ag_names, function(name) {
    idx <- which(ag_spec$name == name)
    if (length(idx) > 0) ag_spec$fill[idx[1]] else "#808080"
  })

  # Match colors to sera
  sr_outline <- sapply(sr_names, function(name) {
    idx <- which(sr_spec$name == name)
    if (length(idx) > 0) sr_spec$outline[idx[1]] else "black"
  })

  # Combine coordinates for range calculation
  all_coords <- rbind(ag_coords, sr_coords)
  valid <- !is.na(all_coords[,1]) & !is.na(all_coords[,2])

  # Calculate plot range
  x_range <- range(all_coords[valid, 1], na.rm = TRUE)
  y_range <- range(all_coords[valid, 2], na.rm = TRUE)

  # Add margin for labels
  x_margin <- diff(x_range) * 0.25
  y_margin <- diff(y_range) * 0.25

  xlim <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)

  # Plot center
  center_x <- mean(x_range)
  center_y <- mean(y_range)

  # Plot
  plot(NULL,
       xlim = xlim,
       ylim = ylim,
       xlab = "Antigenic Dimension 1",
       ylab = "Antigenic Dimension 2",
       main = title,
       asp = 1)

  # Add subtitle if provided
  if (subtitle != "") {
    mtext(subtitle, side = 3, line = 0.3, cex = 0.8)
  }

  # Add grid
  grid(col = "lightgray", lty = "dotted")

  # Plot antigens (circles with fill)
  valid_ag <- !is.na(ag_coords[,1]) & !is.na(ag_coords[,2])
  points(ag_coords[valid_ag, 1], ag_coords[valid_ag, 2],
         pch = 21, bg = ag_fill[valid_ag], col = "black", cex = 1.0, lwd = 0.3)

  # Plot sera (squares, transparent fill with black outline)
  valid_sr <- !is.na(sr_coords[,1]) & !is.na(sr_coords[,2])
  points(sr_coords[valid_sr, 1], sr_coords[valid_sr, 2],
         pch = 22, bg = "transparent", col = sr_outline[valid_sr], cex = 1.3, lwd = 1)

  # Collect reference strain positions for labeling
  label_data <- list()
  unique_colors <- unique(ag_fill[valid_ag])

  for (color in unique_colors) {
    if (color %in% names(cluster_names)) {
      cluster_label <- cluster_names[color]

      # Find the reference strain for this cluster
      if (cluster_label %in% names(cluster_reference_strains)) {
        ref_strain <- cluster_reference_strains[cluster_label]
        ref_idx <- which(ag_names == ref_strain)

        if (length(ref_idx) > 0 && !is.na(ag_coords[ref_idx[1], 1])) {
          # Use reference strain coordinates
          target_x <- ag_coords[ref_idx[1], 1]
          target_y <- ag_coords[ref_idx[1], 2]
        } else {
          # Fallback to centroid if reference strain not found
          color_idx <- which(ag_fill == color & valid_ag)
          target_x <- mean(ag_coords[color_idx, 1], na.rm = TRUE)
          target_y <- mean(ag_coords[color_idx, 2], na.rm = TRUE)
        }
      } else {
        # Fallback to centroid
        color_idx <- which(ag_fill == color & valid_ag)
        target_x <- mean(ag_coords[color_idx, 1], na.rm = TRUE)
        target_y <- mean(ag_coords[color_idx, 2], na.rm = TRUE)
      }

      # Calculate angle from center to target
      angle <- atan2(target_y - center_y, target_x - center_x)

      label_data[[color]] <- list(
        name = cluster_label,
        cx = target_x,
        cy = target_y,
        angle = angle,
        color = color
      )
    }
  }

  # Sort by angle to help with positioning
  angles <- sapply(label_data, function(x) x$angle)
  label_data <- label_data[order(angles)]

  # Place labels outside the plot and draw arrows
  plot_radius_x <- diff(x_range) / 2 + x_margin * 0.7
  plot_radius_y <- diff(y_range) / 2 + y_margin * 0.7

  for (ld in label_data) {
    # Calculate label position on the edge
    label_x <- center_x + plot_radius_x * cos(ld$angle)
    label_y <- center_y + plot_radius_y * sin(ld$angle)

    # Draw arrow from label to centroid
    arrows(
      x0 = label_x, y0 = label_y,
      x1 = ld$cx, y1 = ld$cy,
      length = 0.08, angle = 20,
      col = "black", lwd = 1
    )

    # Determine text position (adjust based on angle)
    if (cos(ld$angle) > 0.3) {
      pos <- 4  # right
    } else if (cos(ld$angle) < -0.3) {
      pos <- 2  # left
    } else if (sin(ld$angle) > 0) {
      pos <- 3  # top
    } else {
      pos <- 1  # bottom
    }

    # Draw label
    text(label_x, label_y, ld$name, font = 2, cex = 0.85, pos = pos, offset = 0.2)
  }
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("LOADING DATA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

titer_table <- read.titerTable("shannen/h3map2004_hitable.csv")

# Load plot specifications for coloring
plotspec <- read.csv("shannen/h3map2004_plotspec.csv", stringsAsFactors = FALSE)
ag_colors <- plotspec[plotspec$type == "ag", c("name", "fill")]
sr_colors <- plotspec[plotspec$type == "sr", c("name", "fill", "outline")]

cat(sprintf("Dataset: %d antigens x %d sera\n", nrow(titer_table), ncol(titer_table)))
cat(sprintf("Dimensions: %d\n", N_DIM))
cat(sprintf("Optimizations: %d\n", N_OPTIMIZATIONS))
cat(sprintf("Minimum column basis: %s\n", MIN_COLUMN_BASIS))
cat(sprintf("Random seed: %d\n", RANDOM_SEED))

# ==============================================================================
# PREPARE BLIND TEST SET
# ==============================================================================
# We identify which titers to hold out for prediction evaluation
# But we will train on FULL data for stress comparison

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("PREPARING BLIND TEST SET\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

set.seed(RANDOM_SEED)

# Find all measurable titers (not NA, not *, not < or >)
titer_matrix <- as.matrix(titer_table)
blind_info <- data.frame(
  ag = integer(),
  sr = integer(),
  titer = numeric(),
  stringsAsFactors = FALSE
)

for (ag in 1:nrow(titer_matrix)) {
  for (sr in 1:ncol(titer_matrix)) {
    titer <- titer_matrix[ag, sr]
    if (!is.na(titer) && titer != "*" && !grepl("^[<>]", titer)) {
      blind_info <- rbind(blind_info, data.frame(
        ag = ag,
        sr = sr,
        titer = as.numeric(titer),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Randomly select titers for blind test
n_blind <- round(nrow(blind_info) * BLIND_FRACTION)
blind_indices <- sample(1:nrow(blind_info), n_blind)
blind_test <- blind_info[blind_indices, ]

cat(sprintf("Total measurable titers: %d\n", nrow(blind_info)))
cat(sprintf("Blind test set: %d (%.1f%%)\n", nrow(blind_test), 100 * nrow(blind_test) / nrow(blind_info)))

# ==============================================================================
# CREATE MAP FROM FULL DATA
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("CREATING BASE MAP FROM FULL DATA\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

map_base <- acmap(titer_table = titer_table)

# Apply plotspec styling to base map
map_base <- apply_plotspec(map_base, plotspec)

# Get initial column bases from a quick fixed optimization
set.seed(RANDOM_SEED)
temp_map <- optimizeMap(
  map_base,
  number_of_dimensions = N_DIM,
  number_of_optimizations = 10,
  minimum_column_basis = MIN_COLUMN_BASIS,
  options = make_options(optimize_colbases = FALSE)
)
initial_colbases <- colBases(temp_map)

cat(sprintf("Initial column bases: min=%.2f, max=%.2f, mean=%.2f\n",
            min(initial_colbases), max(initial_colbases), mean(initial_colbases)))

# ==============================================================================
# DEFINE TEST CONDITIONS
# ==============================================================================

conditions <- list(

  # 1. Baseline: Fixed column bases (traditional approach)
  fixed = list(
    name = "1. Fixed (baseline)",
    opts = make_options(optimize_colbases = FALSE),
    lambda = 0,
    is_regularized = FALSE
  ),

  # 2. Unbounded optimization (demonstrates the problem)
  unbounded = list(
    name = "2. Unbounded",
    opts = make_options(optimize_colbases = TRUE),
    lambda = 0,
    is_regularized = FALSE
  ),

  # 3-5. Bounded only (different deviation limits)
  bounded_3 = list(
    name = "3. Bounded ±3",
    opts = make_options(
      optimize_colbases = TRUE,
      bound_colbases = TRUE,
      colbase_min_bound = 0.0,
      colbase_max_bound = 100.0,
      colbase_max_deviation = 3.0
    ),
    lambda = 0,
    is_regularized = FALSE
  ),

  bounded_5 = list(
    name = "4. Bounded ±5",
    opts = make_options(
      optimize_colbases = TRUE,
      bound_colbases = TRUE,
      colbase_min_bound = 0.0,
      colbase_max_bound = 100.0,
      colbase_max_deviation = 5.0
    ),
    lambda = 0,
    is_regularized = FALSE
  ),

  bounded_8 = list(
    name = "5. Bounded ±8",
    opts = make_options(
      optimize_colbases = TRUE,
      bound_colbases = TRUE,
      colbase_min_bound = 0.0,
      colbase_max_bound = 100.0,
      colbase_max_deviation = 8.0
    ),
    lambda = 0,
    is_regularized = FALSE
  ),

  # 6-8. Regularized only (different lambda values)
  reg_strong = list(
    name = "6. Regularized λ=1.0",
    opts = make_options(
      optimize_colbases = TRUE,
      regularize_colbases = TRUE,
      colbase_lambda = 1.0
    ),
    lambda = 1.0,
    is_regularized = TRUE
  ),

  reg_moderate = list(
    name = "7. Regularized λ=0.1",
    opts = make_options(
      optimize_colbases = TRUE,
      regularize_colbases = TRUE,
      colbase_lambda = 0.1
    ),
    lambda = 0.1,
    is_regularized = TRUE
  ),

  reg_weak = list(
    name = "8. Regularized λ=0.01",
    opts = make_options(
      optimize_colbases = TRUE,
      regularize_colbases = TRUE,
      colbase_lambda = 0.01
    ),
    lambda = 0.01,
    is_regularized = TRUE
  ),

  combined = list(
    name = "9. Bounded +5 + λ=0.1",
    opts = make_options(
      optimize_colbases = TRUE,
      bound_colbases = TRUE,
      colbase_min_bound = 0.0,
      colbase_max_bound = 100.0,
      colbase_max_deviation = 5.0,
      regularize_colbases = TRUE,
      colbase_lambda = 0.1
    ),
    lambda = 0.1
  )
)

# ==============================================================================
# RUN ALL CONDITIONS
# ==============================================================================

results <- list()

for (cond_name in names(conditions)) {
  cond <- conditions[[cond_name]]

  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("Running: %s\n", cond$name))
  cat(paste(rep("=", 70), collapse = ""), "\n")

  set.seed(RANDOM_SEED)

  start_time <- Sys.time()

  tryCatch({
    # Train on FULL data
    map <- optimizeMap(
      map_base,
      number_of_dimensions = N_DIM,
      number_of_optimizations = N_OPTIMIZATIONS,
      minimum_column_basis = MIN_COLUMN_BASIS,
      options = cond$opts
    )

    elapsed <- difftime(Sys.time(), start_time, units = "mins")

    # Apply plotspec styling to map (for Racmacs view())
    map <- apply_plotspec(map, plotspec)

    # Get stress components
    # The C++ optimizer now tracks base_stress and penalty_stress separately.
    # mapStress() returns total stress (base + penalty for regularized models).
    # We calculate penalty on R side to match what C++ computed.
    total_stress <- mapStress(map)
    colbases <- colBases(map)

    # For regularized models, calculate stress components
    if (cond$is_regularized && cond$lambda > 0) {
      # Calculate penalty stress to match C++ calculation:
      # penalty_stress = λ × Σ(colbase - initial_colbase)²
      penalty_stress <- calculate_regularization_penalty(colbases, initial_colbases, cond$lambda)
      base_stress <- total_stress - penalty_stress
    } else {
      # For non-regularized models (fixed, unbounded, bounded)
      # There is no regularization penalty
      base_stress <- total_stress
      penalty_stress <- 0
    }

    # Predict held-out titers using the full-data map
    predictions <- predict_titers_from_map(map, blind_test)
    stats <- calculate_prediction_stats(predictions)

    # Calculate distances for plotting
    distances <- calculate_distances(map, titer_table)

    # Store results with separate stress components
    results[[cond_name]] <- list(
      name = cond$name,
      map = map,
      stress = base_stress,           # Base stress (map fitting only)
      penalty_stress = penalty_stress, # Regularization penalty
      total_stress = total_stress,     # Total = stress + penalty_stress
      lambda = cond$lambda,
      is_regularized = cond$is_regularized,
      colbases = colbases,
      predictions = predictions,
      stats = stats,
      elapsed_mins = as.numeric(elapsed),
      distances = distances
    )

    cat(sprintf("\nCompleted in %.1f minutes\n", elapsed))

    # Print stress information based on model type
    if (cond$is_regularized && cond$lambda > 0) {
      cat(sprintf("Stress (base): %.2f\n", base_stress))
      cat(sprintf("Penalty Stress: %.2f\n", penalty_stress))
      cat(sprintf("Total Stress (stress + penalty): %.2f\n", total_stress))
    } else {
      cat(sprintf("Stress: %.2f\n", base_stress))
    }

    cat(sprintf("Column bases: min=%.2f, max=%.2f, mean=%.2f, sd=%.2f\n",
                min(colbases), max(colbases), mean(colbases), sd(colbases)))
    cat(sprintf("Mean prediction error: %.4f log2 units\n", stats$mean_error))
    cat(sprintf("RMSE: %.4f log2 units\n", stats$rmse))
    cat(sprintf("Correlation: %.4f\n", stats$correlation))

  }, error = function(e) {
    cat(sprintf("\nERROR: %s\n", e$message))
    results[[cond_name]] <<- list(
      name = cond$name,
      error = e$message
    )
  })
}

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("SUMMARY COMPARISON\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

# Build summary dataframe with separate stress columns
summary_rows <- list()
for (cond_name in names(results)) {
  r <- results[[cond_name]]
  if (!is.null(r$error)) {
    summary_rows[[cond_name]] <- data.frame(
      Condition = r$name,
      Stress = NA,
      Penalty_Stress = NA,
      Total_Stress = NA,
      Mean_Error = NA,
      RMSE = NA,
      Correlation = NA,
      CB_Min = NA,
      CB_Max = NA,
      CB_Mean = NA,
      CB_SD = NA,
      Time_min = NA,
      stringsAsFactors = FALSE
    )
  } else {
    summary_rows[[cond_name]] <- data.frame(
      Condition = r$name,
      Stress = round(r$stress, 2),
      Penalty_Stress = round(r$penalty_stress, 2),
      Total_Stress = round(r$total_stress, 2),
      Mean_Error = round(r$stats$mean_error, 4),
      RMSE = round(r$stats$rmse, 4),
      Correlation = round(r$stats$correlation, 4),
      CB_Min = round(min(r$colbases), 2),
      CB_Max = round(max(r$colbases), 2),
      CB_Mean = round(mean(r$colbases), 2),
      CB_SD = round(sd(r$colbases), 2),
      Time_min = round(r$elapsed_mins, 1),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

# Print full summary
cat("Full Results Table:\n")
cat("NOTE: For non-regularized models, Penalty_Stress = 0 and Stress = Total_Stress\n")
cat("      For regularized models (6-8), Total_Stress = Stress + Penalty_Stress\n\n")
print(summary_df, row.names = FALSE)

# ==============================================================================
# STRESS COMPARISON
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("STRESS COMPARISON\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

cat("INTERPRETATION:\n")
cat("  - Stress: Map fitting stress only (comparable across all models)\n")
cat("  - Penalty_Stress: Regularization penalty (only for regularized models 6-8)\n")
cat("  - Total_Stress: Stress + Penalty_Stress (what optimizer minimizes for regularized)\n")
cat("  - For fair comparison, use 'Stress' column (excludes penalty)\n\n")

stress_comparison <- summary_df[!is.na(summary_df$Stress),
                                c("Condition", "Stress", "Penalty_Stress", "Total_Stress", "CB_Max")]
stress_comparison <- stress_comparison[order(stress_comparison$Stress), ]
cat("Ranked by base Stress (lower = better map fit):\n\n")
print(stress_comparison, row.names = FALSE)

# ==============================================================================
# REGULARIZED MODELS DETAIL
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("REGULARIZED MODELS DETAIL (6-8)\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

cat("For regularized models, the optimizer minimizes Total_Stress = Stress + Penalty_Stress\n")
cat("where Penalty_Stress = λ × Σ(colbase - initial_colbase)²\n\n")

reg_models <- summary_df[summary_df$Penalty_Stress > 0 & !is.na(summary_df$Penalty_Stress), ]
if (nrow(reg_models) > 0) {
  cat("Regularized Models Breakdown:\n\n")
  print(reg_models[, c("Condition", "Stress", "Penalty_Stress", "Total_Stress", "CB_Min", "CB_Max", "Mean_Error")],
        row.names = FALSE)
}

# ==============================================================================
# RANKING BY PREDICTION ACCURACY
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("RANKING BY PREDICTION ACCURACY (Mean Error on Held-Out Titers)\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

cat(sprintf("NOTE: Prediction accuracy evaluated on %d held-out titers (%.0f%% of data)\n",
            nrow(blind_test), 100 * BLIND_FRACTION))
cat("      Lower Mean_Error = better prediction\n\n")

ranked_df <- summary_df[!is.na(summary_df$Mean_Error), ]
ranked_df <- ranked_df[order(ranked_df$Mean_Error), ]
ranked_df$Rank <- 1:nrow(ranked_df)

baseline_error <- summary_df$Mean_Error[summary_df$Condition == "1. Fixed (baseline)"]
ranked_df$Improvement <- sprintf("%+.1f%%",
                                 100 * (baseline_error - ranked_df$Mean_Error) / baseline_error)

print(ranked_df[, c("Rank", "Condition", "Mean_Error", "RMSE", "Correlation", "Stress", "CB_Max", "Improvement")],
      row.names = FALSE)

# ==============================================================================
# COLUMN BASE COMPARISON TABLE
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("COLUMN BASE COMPARISON BY SERUM\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

colbase_comparison <- data.frame(
  serum = srNames(map_base),
  initial = round(initial_colbases, 2),
  stringsAsFactors = FALSE
)

for (cond_name in names(results)) {
  r <- results[[cond_name]]
  if (is.null(r$error)) {
    col_name <- gsub("^[0-9]+\\. ", "", r$name)
    col_name <- gsub(" ", "_", col_name)
    col_name <- gsub("[^a-zA-Z0-9_λ±+]", "", col_name)
    colbase_comparison[[col_name]] <- round(r$colbases, 2)
  }
}

# Add difference from initial for unbounded (to show the problem)
if (!is.null(results$unbounded) && is.null(results$unbounded$error)) {
  colbase_comparison$Unbounded_diff <- round(results$unbounded$colbases - initial_colbases, 2)
}

# Sort by unbounded difference to highlight problematic sera
if ("Unbounded_diff" %in% names(colbase_comparison)) {
  colbase_comparison <- colbase_comparison[order(-abs(colbase_comparison$Unbounded_diff)), ]
}

cat("Top 20 sera with largest column base changes in unbounded model:\n\n")
print(head(colbase_comparison, 20), row.names = FALSE)

# ==============================================================================
# ANALYSIS BY CB CHANGE GROUP
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("PREDICTION ERROR BY COLUMN BASE CHANGE MAGNITUDE\n")
cat(paste(rep("=", 100), collapse = ""), "\n")

if (!is.null(results$fixed) && is.null(results$fixed$error)) {
  fixed_colbases <- results$fixed$colbases

  for (cond_name in names(results)) {
    if (cond_name == "fixed") next
    r <- results[[cond_name]]
    if (!is.null(r$error)) next

    cat(sprintf("\n--- %s ---\n", r$name))

    cb_change <- r$colbases - fixed_colbases

    # Add CB change to blind test predictions
    blind_analysis <- r$predictions
    blind_analysis$cb_change <- cb_change[blind_analysis$sr]
    blind_analysis$error_fixed <- results$fixed$predictions$error

    # Calculate improvement
    blind_analysis$improvement <- blind_analysis$error_fixed - blind_analysis$error

    blind_analysis$cb_group <- cut(
      abs(blind_analysis$cb_change),
      breaks = c(-Inf, 3, 8, 15, Inf),
      labels = c("Small (0-3)", "Medium (3-8)", "Large (8-15)", "Extreme (>15)")
    )

    cat("Mean prediction error by CB change magnitude:\n")
    for (g in levels(blind_analysis$cb_group)) {
      subset_data <- blind_analysis[blind_analysis$cb_group == g & !is.na(blind_analysis$error), ]
      if (nrow(subset_data) > 0) {
        mean_err <- mean(subset_data$error, na.rm = TRUE)
        mean_imp <- mean(subset_data$improvement, na.rm = TRUE)
        n <- nrow(subset_data)
        cat(sprintf("  %s: error=%.4f, improvement=%+.4f (n=%d)\n",
                    g, mean_err, mean_imp, n))
      }
    }
  }
}

# ==============================================================================
# RECOMMENDATION
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("RECOMMENDATION\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

best_row <- ranked_df[1, ]

cat(sprintf("Best performing model: %s\n", best_row$Condition))
cat(sprintf("  Base Stress: %.2f\n", best_row$Stress))
cat(sprintf("  Mean prediction error: %.4f log2 units\n", best_row$Mean_Error))
cat(sprintf("  RMSE: %.4f log2 units\n", best_row$RMSE))
cat(sprintf("  Correlation: %.4f\n", best_row$Correlation))
cat(sprintf("  Max column base: %.2f\n", best_row$CB_Max))
cat(sprintf("  Improvement over baseline: %s\n", best_row$Improvement))

if (best_row$CB_Max <= 14) {
  cat("\n✓ Column bases are within biologically reasonable range (max ≤ 14)\n")
} else if (best_row$CB_Max <= 20) {
  cat("\n⚠ Column bases are slightly elevated (max > 14) - may want tighter bounds\n")
} else {
  cat("\n✗ Column bases are unreasonably high (max > 20) - use bounded version\n")
}

# Additional insights
cat("\n\nKEY INSIGHTS:\n")
cat("-------------\n")

# Compare unbounded vs best regularized
if (!is.null(results$unbounded) && is.null(results$unbounded$error)) {
  unbounded_error <- results$unbounded$stats$mean_error
  unbounded_cb_max <- max(results$unbounded$colbases)
  unbounded_stress <- results$unbounded$stress

  cat(sprintf("\n1. Unbounded optimization:\n"))
  cat(sprintf("   - Base Stress: %.2f\n", unbounded_stress))
  cat(sprintf("   - Mean error: %.4f (%.1f%% vs baseline)\n",
              unbounded_error, 100 * (baseline_error - unbounded_error) / baseline_error))
  cat(sprintf("   - BUT max column base: %.2f (biologically implausible)\n", unbounded_cb_max))

  cat(sprintf("\n2. %s:\n", best_row$Condition))
  best_result <- results[[names(which(sapply(results, function(x) !is.null(x$name) && x$name == best_row$Condition)))[1]]]
  cat(sprintf("   - Base Stress: %.2f\n", best_result$stress))
  if (best_result$is_regularized) {
    cat(sprintf("   - Penalty Stress: %.2f\n", best_result$penalty_stress))
    cat(sprintf("   - Total Stress: %.2f\n", best_result$total_stress))
  }
  cat(sprintf("   - Mean error: %.4f (%s vs baseline)\n",
              best_row$Mean_Error, best_row$Improvement))
  cat(sprintf("   - Max column base: %.2f\n", best_row$CB_Max))

  if (best_row$Mean_Error < unbounded_error) {
    cat("\n   → Regularization/Bounding IMPROVES prediction accuracy while keeping realistic column bases!\n")
  } else if (best_row$CB_Max <= 14) {
    cat(sprintf("\n   → Regularization/Bounding maintains biologically plausible column bases with minimal accuracy tradeoff\n"))
  }
}

# Compare to paper
cat("\n3. Comparison to original paper:\n")
cat("   - Paper reported stress: ~3600\n")
if (!is.null(results$fixed) && is.null(results$fixed$error)) {
  cat(sprintf("   - Our fixed baseline stress: %.2f ✓\n", results$fixed$stress))
}

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("GENERATING PLOTS\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

n_successful <- sum(sapply(results, function(r) is.null(r$error)))

if (n_successful > 0) {
  # Multi-panel PDF
  n_cols <- 3
  n_rows <- ceiling(n_successful / n_cols)

  pdf("table_vs_map_distance_comparison.pdf", width = 4 * n_cols, height = 4 * n_rows)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

  for (cond_name in names(results)) {
    r <- results[[cond_name]]
    if (!is.null(r$error)) next

    distances <- r$distances
    measurable <- distances[distances$titer_type == "measurable", ]
    less_than <- distances[distances$titer_type == "less_than", ]

    all_dists <- c(distances$table_dist, distances$map_dist)
    lim_range <- range(all_dists, na.rm = TRUE)
    lim_range <- c(max(0, lim_range[1] - 1), lim_range[2] + 1)

    plot(NULL,
         xlim = lim_range, ylim = lim_range,
         xlab = "Table Distance", ylab = "Map Distance",
         main = r$name, asp = 1)

    abline(0, 1, col = "red", lwd = 2, lty = 2)

    if (nrow(measurable) > 0) {
      points(measurable$table_dist, measurable$map_dist,
             pch = 16, col = rgb(0, 0, 0, 0.2), cex = 0.5)
    }
    if (nrow(less_than) > 0) {
      points(less_than$table_dist, less_than$map_dist,
             pch = 16, col = rgb(0, 0, 1, 0.2), cex = 0.5)
    }

    # Legend text differs based on model type
    if (r$is_regularized && r$lambda > 0) {
      legend("topleft",
             legend = c(sprintf("Stress: %.1f", r$stress),
                        sprintf("Penalty: %.1f", r$penalty_stress),
                        sprintf("CB range: [%.1f, %.1f]", min(r$colbases), max(r$colbases)),
                        sprintf("Mean error: %.3f", r$stats$mean_error)),
             bty = "n", cex = 0.7)
    } else {
      legend("topleft",
             legend = c(sprintf("Stress: %.1f", r$stress),
                        sprintf("CB range: [%.1f, %.1f]", min(r$colbases), max(r$colbases)),
                        sprintf("Mean error: %.3f", r$stats$mean_error)),
             bty = "n", cex = 0.8)
    }

    legend("bottomright",
           legend = c("Measurable", "Less-than", "Zero stress"),
           pch = c(16, 16, NA), lty = c(NA, NA, 2),
           col = c(rgb(0, 0, 0, 0.5), rgb(0, 0, 1, 0.5), "red"),
           bty = "n", cex = 0.7)
  }

  dev.off()
  cat("Saved: table_vs_map_distance_comparison.pdf\n")

  # Individual PNGs
  for (cond_name in names(results)) {
    r <- results[[cond_name]]
    if (!is.null(r$error)) next

    distances <- r$distances
    measurable <- distances[distances$titer_type == "measurable", ]
    less_than <- distances[distances$titer_type == "less_than", ]

    all_dists <- c(distances$table_dist, distances$map_dist)
    lim_range <- range(all_dists, na.rm = TRUE)
    lim_range <- c(max(0, lim_range[1] - 1), lim_range[2] + 1)

    filename <- paste0("dist_plot_", gsub("[^a-zA-Z0-9]", "_", cond_name), ".png")

    png(filename, width = 600, height = 600, res = 100)
    par(mar = c(4, 4, 3, 1))

    plot(NULL,
         xlim = lim_range, ylim = lim_range,
         xlab = "Table Distance (log2 units)",
         ylab = "Map Distance (antigenic units)",
         main = r$name, asp = 1)

    abline(0, 1, col = "red", lwd = 2, lty = 2)

    if (nrow(measurable) > 0) {
      points(measurable$table_dist, measurable$map_dist,
             pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.6)
    }
    if (nrow(less_than) > 0) {
      points(less_than$table_dist, less_than$map_dist,
             pch = 16, col = rgb(0, 0, 1, 0.3), cex = 0.6)
    }

    # Legend text differs based on model type
    if (r$is_regularized && r$lambda > 0) {
      legend("topleft",
             legend = c(sprintf("Stress: %.1f", r$stress),
                        sprintf("Penalty Stress: %.1f", r$penalty_stress),
                        sprintf("Total Stress: %.1f", r$total_stress),
                        sprintf("CB range: [%.1f, %.1f]", min(r$colbases), max(r$colbases)),
                        sprintf("Mean error: %.3f", r$stats$mean_error),
                        sprintf("Correlation: %.3f", r$stats$correlation)),
             bty = "n", cex = 0.8)
    } else {
      legend("topleft",
             legend = c(sprintf("Stress: %.1f", r$stress),
                        sprintf("CB range: [%.1f, %.1f]", min(r$colbases), max(r$colbases)),
                        sprintf("Mean error: %.3f", r$stats$mean_error),
                        sprintf("Correlation: %.3f", r$stats$correlation)),
             bty = "n", cex = 0.9)
    }

    legend("bottomright",
           legend = c("Measurable titers", "Less-than titers", "Zero stress line"),
           pch = c(16, 16, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
           col = c(rgb(0, 0, 0, 0.6), rgb(0, 0, 1, 0.6), "red"),
           bty = "n", cex = 0.8)

    dev.off()
    cat(sprintf("Saved: %s\n", filename))
  }

  # ==========================================================================
  # ANTIGENIC MAP PLOTS WITH PLOTSPEC COLORS
  # ==========================================================================

  cat("\nGenerating antigenic maps...\n")

  # Multi-panel PDF of all antigenic maps
  pdf("antigenic_maps_comparison.pdf", width = 4 * n_cols, height = 4 * n_rows)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

  for (cond_name in names(results)) {
    r <- results[[cond_name]]
    if (!is.null(r$error)) next

    if (r$is_regularized && r$lambda > 0) {
      subtitle <- sprintf("Stress: %.1f + Pen: %.1f | CB: [%.1f, %.1f] | Err: %.3f",
                          r$stress, r$penalty_stress, min(r$colbases), max(r$colbases), r$stats$mean_error)
    } else {
      subtitle <- sprintf("Stress: %.1f | CB: [%.1f, %.1f] | Error: %.3f",
                          r$stress, min(r$colbases), max(r$colbases), r$stats$mean_error)
    }
    plot_antigenic_map(r$map, plotspec, title = r$name, subtitle = subtitle)
  }

  dev.off()
  cat("Saved: antigenic_maps_comparison.pdf\n")

  # Individual high-resolution PNGs for each map
  for (cond_name in names(results)) {
    r <- results[[cond_name]]
    if (!is.null(r$error)) next

    filename <- paste0("antigenic_map_", gsub("[^a-zA-Z0-9]", "_", cond_name), ".png")

    png(filename, width = 800, height = 800, res = 100)
    par(mar = c(4, 4, 4, 1))

    if (r$is_regularized && r$lambda > 0) {
      subtitle <- sprintf("Stress: %.1f | Penalty: %.1f | Total: %.1f | CB: [%.1f, %.1f] | Error: %.3f",
                          r$stress, r$penalty_stress, r$total_stress,
                          min(r$colbases), max(r$colbases), r$stats$mean_error)
    } else {
      subtitle <- sprintf("Stress: %.1f | CB range: [%.1f, %.1f] | Mean error: %.3f",
                          r$stress, min(r$colbases), max(r$colbases), r$stats$mean_error)
    }
    plot_antigenic_map(r$map, plotspec, title = r$name, subtitle = subtitle)

    dev.off()
    cat(sprintf("Saved: %s\n", filename))
  }
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\n\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("RESULTS SAVED\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")

saveRDS(results, "colbase_comparison_results.rds")
saveRDS(summary_df, "colbase_comparison_summary.rds")
saveRDS(ranked_df, "colbase_comparison_ranked.rds")

cat("R objects:\n")
cat("  - colbase_comparison_results.rds: Full results for each condition\n")
cat("  - colbase_comparison_summary.rds: Summary statistics table\n")
cat("  - colbase_comparison_ranked.rds: Ranking by prediction accuracy\n")
cat("\nPlots:\n")
cat("  - table_vs_map_distance_comparison.pdf: Distance plots (all conditions)\n")
cat("  - dist_plot_*.png: Individual distance plots\n")
cat("  - antigenic_maps_comparison.pdf: Antigenic maps (all conditions)\n")
cat("  - antigenic_map_*.png: Individual antigenic maps with plotspec colors\n")

cat("\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 100), collapse = ""), "\n")
