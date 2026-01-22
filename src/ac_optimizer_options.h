
#ifndef Racmacs__ac_optimizer_options__h
#define Racmacs__ac_optimizer_options__h

// Optimizer options
struct AcOptimizerOptions {

  bool dim_annealing;
  std::string method;
  int maxit;
  int num_basis;
  double armijo_constant;
  double wolfe;
  double min_gradient_norm;
  double factr;
  int max_line_search_trials;
  double min_step;
  double max_step;
  int num_cores;
  bool report_progress;
  int progress_bar_length;

  // Column base optimization options
  bool optimize_colbases = false;
  
  // NEW: Bounded column bases
  bool bound_colbases = false;           // Enable hard bounds on column bases
  double colbase_min_bound = 4.0;        // Minimum column base (log2 scale, 4 = titer 160)
  double colbase_max_bound = 14.0;       // Maximum column base (log2 scale, 14 = titer 163840)
  double colbase_max_deviation = 5.0;    // Max deviation from initial value (0 = no limit)
  
  // NEW: Regularized column bases  
  bool regularize_colbases = false;      // Enable regularization penalty
  double colbase_lambda = 0.1;           // Regularization strength (0 = no penalty)
};

#endif