# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

make_options <- function(optimize_colbases = FALSE) {
  opts <- RacOptimizer.options()
  opts$optimize_colbases <- optimize_colbases
  opts
}

#' Extract measurable titers from a titer table
extract_measurable_titers <- function(titer_matrix) {
  measurable <- list()
  idx <- 1
  
  for (i in 1:nrow(titer_matrix)) {
    for (j in 1:ncol(titer_matrix)) {
      titer <- titer_matrix[i, j]
      # Check if titer is numeric (not *, <, >, etc.)
      if (!is.na(titer) && grepl("^[0-9]+$", as.character(titer))) {
        measurable[[idx]] <- list(ag = i, sr = j, titer = titer)
        idx <- idx + 1
      }
    }
  }
  
  return(measurable)
}

#' Prepare blind test data by masking a fraction of measurable titers
prepare_blind_test <- function(titer_matrix, blind_fraction = 0.2, min_blind = 10) {
  cat("\n=== PREPARING BLIND TEST DATA ===\n")
  
  # Extract measurable titers
  measurable_titers <- extract_measurable_titers(titer_matrix)
  n_total <- length(measurable_titers)
  
  if (n_total == 0) {
    stop("No measurable titers found in the dataset")
  }
  
  # Calculate blind test size
  n_blind <- max(min_blind, floor(n_total * blind_fraction))
  n_blind <- min(n_blind, floor(n_total * 0.5))  # Cap at 50%
  
  # Select random titers for blind testing
  blind_indices <- sample(n_total, n_blind)
  
  cat(sprintf("Data split:\n"))
  cat(sprintf("  Total measurable titers: %d\n", n_total))
  cat(sprintf("  Training: %d (%.1f%%)\n", 
              n_total - n_blind, 
              100 * (n_total - n_blind) / n_total))
  cat(sprintf("  Blind test: %d (%.1f%%)\n", 
              n_blind, 
              100 * n_blind / n_total))
  
  # Create training data with masked titers
  train_titers <- titer_matrix
  blind_data <- data.frame()
  
  for (idx in blind_indices) {
    pos <- measurable_titers[[idx]]
    train_titers[pos$ag, pos$sr] <- "*"
    
    numeric_titer <- as.numeric(pos$titer)
    blind_data <- rbind(blind_data, 
                        data.frame(ag = pos$ag, 
                                   sr = pos$sr, 
                                   actual = numeric_titer))
  }
  
  return(list(
    train = train_titers, 
    blind = blind_data,
    n_train = n_total - nrow(blind_data),
    n_blind = nrow(blind_data)
  ))
}

#' Predict titers using map coordinates and column bases
predict_titers <- function(map, blind_data, col_bases = NULL, max_distance = 10) {
  ag_coords <- agCoords(map)
  sr_coords <- srCoords(map)
  
  # Use provided column bases or get from map
  if (is.null(col_bases)) {
    col_bases <- colBases(map)
  }
  
  n <- nrow(blind_data)
  predictions <- numeric(n)
  errors <- numeric(n)
  valid <- logical(n)
  distances <- numeric(n)
  
  for (i in 1:n) {
    ag_idx <- blind_data$ag[i]
    sr_idx <- blind_data$sr[i]
    
    ag_pos <- ag_coords[ag_idx, ]
    sr_pos <- sr_coords[sr_idx, ]
    
    # Check if positions are valid
    if (any(is.na(ag_pos)) || any(is.na(sr_pos))) {
      predictions[i] <- NA
      errors[i] <- NA
      valid[i] <- FALSE
      distances[i] <- NA
    } else {
      # Calculate map distance
      distance <- sqrt(sum((ag_pos - sr_pos)^2))
      distances[i] <- distance
      
      # Skip if distance is too large
      if (distance > max_distance) {
        predictions[i] <- NA
        errors[i] <- NA
        valid[i] <- FALSE
      } else {
        # Predict titer: predicted_log2 = column_base - distance
        predicted_log2 <- col_bases[sr_idx] - distance
        predicted_log2 <- max(0, min(predicted_log2, 15))  # Reasonable bounds
        
        predictions[i] <- 10 * 2^predicted_log2
        
        # Calculate error in log2 units
        actual_log2 <- log2(blind_data$actual[i] / 10)
        errors[i] <- abs(predicted_log2 - actual_log2)
        valid[i] <- TRUE
      }
    }
  }
  
  return(list(
    predictions = predictions, 
    errors = errors, 
    valid = valid,
    distances = distances
  ))
}

#' Calculate error statistics from predictions
calculate_error_stats <- function(pred_results, blind_data) {
  valid <- pred_results$valid
  
  if (sum(valid) == 0) {
    return(list(
      mean_error = NA,
      median_error = NA,
      sd_error = NA,
      rmse = NA,
      correlation = NA,
      n_valid = 0,
      n_invalid = sum(!valid)
    ))
  }
  
  errors <- pred_results$errors[valid]
  pred_log2 <- log2(pred_results$predictions[valid] / 10)
  actual_log2 <- log2(blind_data$actual[valid] / 10)
  
  # Calculate correlation if there's variance
  if (sd(pred_log2) > 0 && sd(actual_log2) > 0) {
    correlation <- cor(pred_log2, actual_log2)
  } else {
    correlation <- NA
  }
  
  return(list(
    mean_error = mean(errors),
    median_error = median(errors),
    sd_error = sd(errors),
    rmse = sqrt(mean(errors^2)),
    correlation = correlation,
    n_valid = sum(valid),
    n_invalid = sum(!valid)
  ))
}