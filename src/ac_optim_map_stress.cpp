// src/ac_optim_map_stress.cpp
// CORRECTED VERSION - Fixed forward declarations and removed duplicate functions

#include <math.h>
#include <RcppArmadillo.h>
#include <RcppEnsmallen.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]

#include "utils.h"
#include "utils_error.h"
#include "utils_progress.h"
#include "acmap_map.h"
#include "ac_stress.h"
#include "ac_optim_map_stress.h"
#include "ac_optimization.h"
#include "ac_optimizer_options.h"
#include "acmap_optimization.h"
#include "acmap_titers.h"


// ============================================================================
// FORWARD DECLARATIONS (must come before they are used)
// ============================================================================

arma::mat create_initial_pars(
    const arma::mat &ag_coords,
    const arma::mat &sr_coords,
    const arma::uvec &moveable_ags,
    const arma::uvec &moveable_sr,
    const arma::vec &colbases,
    const arma::uvec &moveable_colbases,
    bool optimize_colbases
);

void ac_relaxOptimizations(
  std::vector<AcOptimization>& optimizations,
  arma::uword num_dims,
  const arma::mat &tabledist_matrix,
  const arma::imat &titertype_matrix,
  const AcOptimizerOptions &options,
  const arma::mat &titer_weights,
  const double &dilution_stepsize
);

void ac_relaxOptimizations_with_colbases(
  std::vector<AcOptimization>& optimizations,
  arma::uword num_dims,
  const arma::mat &logtiter_matrix,
  const arma::imat &titertype_matrix,
  const arma::vec &initial_colbases,
  const AcOptimizerOptions &options,
  const arma::mat &titer_weights,
  const double &dilution_stepsize
);

void sort_optimizations_by_stress(
    std::vector<AcOptimization>& optimizations
);

void align_optimizations(
  std::vector<AcOptimization>& optimizations
);


// ============================================================================
// MapOptimizer CLASS
// ============================================================================

class MapOptimizer {

  public:

    // > ATTRIBUTES

    // - Coordinate matrices
    arma::mat ag_coords;
    arma::mat sr_coords;

    // - Distance matrices
    arma::mat mapdist_matrix;      // Distances from coordinates (computed)
    arma::mat tabledist_matrix;    // Distances from titers (now computed dynamically)

    // - NEW: Store log titers separately (constant throughout optimization)
    arma::mat logtiter_matrix;     // FIXED - the log titers (constant)

    // - NEW: Column bases (optimizable when optimize_colbases = true)
    arma::vec colbases;            // OPTIMIZABLE - column bases
    arma::vec colbase_gradients;
    arma::uvec moveable_colbases;  // Indices of colbases to optimize

    // NEW: Flag to control column base optimization
    bool optimize_colbases;

    // Titer information
    arma::imat titertype_matrix;
    
    // Dimensions
    arma::uword num_dims;
    arma::uword num_ags;
    arma::uword num_sr;
    
    // Point tracking
    arma::uvec moveable_ags;
    arma::uvec moveable_sr;
    arma::uvec included_ags;
    arma::uvec included_srs;
    arma::uvec::iterator agi;
    arma::uvec::iterator agi_end;
    arma::uvec::iterator sri;
    arma::uvec::iterator sri_end;
    
    // Weights and gradients
    arma::mat titer_weights;
    arma::mat ag_gradients;
    arma::mat sr_gradients;
    
    // Other parameters
    double dilution_stepsize;
    double gradient;
    double stress;

    // ========================================================================
    // CONSTRUCTOR 1: Basic constructor (backward compatible, no colbase optimization)
    // ========================================================================
    MapOptimizer(
      arma::mat ag_start_coords,
      arma::mat sr_start_coords,
      arma::mat tabledist,           // Pre-computed table distances
      arma::imat titertype,
      arma::uword dims,
      double dilution_stepsize_in
    )
      :ag_coords(ag_start_coords),
       sr_coords(sr_start_coords),
       tabledist_matrix(tabledist),
       titertype_matrix(titertype),
       num_dims(dims),
       num_ags(tabledist.n_rows),
       num_sr(tabledist.n_cols),
       dilution_stepsize(dilution_stepsize_in),
       optimize_colbases(false)
    {
      // Set default moveable antigens and sera to all
      moveable_ags = arma::regspace<arma::uvec>(0, num_ags - 1);
      moveable_sr = arma::regspace<arma::uvec>(0, num_sr - 1);

      // Set included antigens and sera
      included_ags = arma::find_finite(ag_start_coords.col(0));
      included_srs = arma::find_finite(sr_start_coords.col(0));
      agi_end = included_ags.end();
      sri_end = included_srs.end();

      // Set default weights to 1
      titer_weights.ones(num_ags, num_sr);

      // Setup map dist matrices
      mapdist_matrix = arma::mat(num_ags, num_sr, arma::fill::zeros);

      // Setup the gradient vectors
      ag_gradients.zeros(num_ags, num_dims);
      sr_gradients.zeros(num_sr, num_dims);

      // Update the map distance matrix according to coordinates
      update_map_dist_matrix();
    }

    // ========================================================================
    // CONSTRUCTOR 2: With fixed points (backward compatible)
    // ========================================================================
    MapOptimizer(
      arma::mat ag_start_coords,
      arma::mat sr_start_coords,
      arma::mat tabledist,
      arma::imat titertype,
      arma::uword dims,
      arma::uvec ag_fixed,
      arma::uvec sr_fixed,
      arma::mat titer_weights_in,
      double dilution_stepsize_in
    )
      :ag_coords(ag_start_coords),
       sr_coords(sr_start_coords),
       tabledist_matrix(tabledist),
       titertype_matrix(titertype),
       num_dims(dims),
       num_ags(tabledist.n_rows),
       num_sr(tabledist.n_cols),
       dilution_stepsize(dilution_stepsize_in),
       optimize_colbases(false)
    {
      // Set default weights to 1 if missing
      if (titer_weights_in.n_elem == 0) titer_weights.ones(num_ags, num_sr);
      else                              titer_weights = titer_weights_in;

      // Set moveable antigens
      moveable_ags = arma::find(ag_fixed == 0);
      moveable_sr = arma::find(sr_fixed == 0);

      // Set included antigens and sera
      included_ags = arma::find_finite(ag_start_coords.col(0));
      included_srs = arma::find_finite(sr_start_coords.col(0));
      agi_end = included_ags.end();
      sri_end = included_srs.end();

      // Setup map dist matrices
      mapdist_matrix = arma::mat(num_ags, num_sr, arma::fill::zeros);

      // Setup the gradient vectors
      ag_gradients.zeros(num_ags, num_dims);
      sr_gradients.zeros(num_sr, num_dims);

      // Update the map distance matrix according to coordinates
      update_map_dist_matrix();
    }

    // ========================================================================
    // CONSTRUCTOR 3: NEW - With column base optimization support
    // ========================================================================
    MapOptimizer(
      arma::mat ag_start_coords,
      arma::mat sr_start_coords,
      arma::mat logtiter_matrix_in,
      arma::imat titertype,
      arma::vec initial_colbases,
      arma::uword dims,
      arma::uvec ag_fixed,
      arma::uvec sr_fixed,
      arma::mat titer_weights_in,
      double dilution_stepsize_in,
      bool optimize_colbases_in
    )
      :ag_coords(ag_start_coords),
       sr_coords(sr_start_coords),
       logtiter_matrix(logtiter_matrix_in),
       colbases(initial_colbases),
       titertype_matrix(titertype),
       num_dims(dims),
       num_ags(logtiter_matrix_in.n_rows),
       num_sr(logtiter_matrix_in.n_cols),
       dilution_stepsize(dilution_stepsize_in),
       optimize_colbases(optimize_colbases_in)
    {
      // Set default weights
      if (titer_weights_in.n_elem == 0) titer_weights.ones(num_ags, num_sr);
      else                              titer_weights = titer_weights_in;

      // Set moveable antigens and sera
      moveable_ags = arma::find(ag_fixed == 0);
      moveable_sr = arma::find(sr_fixed == 0);

      // Set included antigens and sera
      included_ags = arma::find_finite(ag_start_coords.col(0));
      included_srs = arma::find_finite(sr_start_coords.col(0));
      agi_end = included_ags.end();
      sri_end = included_srs.end();

      // All column bases are moveable when optimize_colbases is true
      if (optimize_colbases) {
        moveable_colbases = arma::regspace<arma::uvec>(0, num_sr - 1);
        colbase_gradients.zeros(num_sr);
      }

      // Setup distance matrices
      mapdist_matrix = arma::mat(num_ags, num_sr, arma::fill::zeros);
      tabledist_matrix = arma::mat(num_ags, num_sr, arma::fill::zeros);

      // Setup the gradient vectors
      ag_gradients.zeros(num_ags, num_dims);
      sr_gradients.zeros(num_sr, num_dims);

      // Initialize table distances from log titers and initial colbases
      update_table_dist_matrix();

      // Update the map distance matrix according to coordinates
      update_map_dist_matrix();
    }

    // ========================================================================
    // EVALUATE OBJECTIVE FUNCTION
    // ========================================================================
    double Evaluate(const arma::mat &pars) {
      update_parameters(pars);
      if (optimize_colbases) {
        update_table_dist_matrix();
      }
      update_map_dist_matrix();
      return calculate_stress();
    }

    // ========================================================================
    // EVALUATE WITH GRADIENT
    // ========================================================================
    double EvaluateWithGradient(const arma::mat &pars, arma::mat &grad) {
      update_parameters(pars);

      if (optimize_colbases) {
        update_table_dist_matrix();
      }

      update_map_dist_matrix();
      update_gradients();
      
      if (optimize_colbases) {
        update_colbase_gradients();
      }

      arma::uword coord_rows = moveable_ags.n_elem + moveable_sr.n_elem;
      arma::uword total_rows = coord_rows;
      if (optimize_colbases) {
        total_rows += moveable_colbases.n_elem;
      }
      grad.set_size(total_rows, num_dims);
      
      for (arma::uword i = 0; i < moveable_ags.n_elem; ++i) {
        grad.row(i) = ag_gradients.row(moveable_ags(i));
      }
      for (arma::uword i = 0; i < moveable_sr.n_elem; ++i) {
        grad.row(i + moveable_ags.n_elem) = sr_gradients.row(moveable_sr(i));
      }

      if (optimize_colbases) {
        for (arma::uword i = 0; i < moveable_colbases.n_elem; ++i) {
          grad.at(coord_rows + i, 0) = colbase_gradients(moveable_colbases(i));
          for (arma::uword j = 1; j < num_dims; ++j) {
            grad.at(coord_rows + i, j) = 0.0;
          }
        }
      }

      return calculate_stress();
    }

    // ========================================================================
    // UPDATE GRADIENTS (for coordinates)
    // ========================================================================
    void update_gradients() {
      ag_gradients.zeros();
      sr_gradients.zeros();

      for(sri = included_srs.begin(); sri != sri_end; ++sri) {
        for(agi = included_ags.begin(); agi != agi_end; ++agi) {
          if(titertype_matrix.at(*agi, *sri) <= 0) continue;

          double ibase = titer_weights.at(*agi,*sri) * inc_base(
            mapdist_matrix.at(*agi, *sri),
            tabledist_matrix.at(*agi, *sri),
            titertype_matrix.at(*agi, *sri),
            dilution_stepsize
          );

          for(arma::uword i = 0; i < num_dims; ++i) {
            gradient = ibase*(ag_coords.at(*agi, i) - sr_coords.at(*sri, i));
            ag_gradients.at(*agi, i) -= gradient;
            sr_gradients.at(*sri, i) += gradient;
          }
        }
      }
    }

    // ========================================================================
    // UPDATE COLUMN BASE GRADIENTS
    // ========================================================================
    void update_colbase_gradients() {
      colbase_gradients.zeros();

      for (arma::uword sr = 0; sr < num_sr; ++sr) {
        for (arma::uword ag = 0; ag < num_ags; ++ag) {
          arma::sword titer_type = titertype_matrix.at(ag, sr);
          if (titer_type <= 0) continue;

          double table_dist = tabledist_matrix(ag, sr);
          double map_dist = mapdist_matrix(ag, sr);
          double weight = titer_weights(ag, sr);

          if (titer_type == 1) {
            colbase_gradients(sr) += weight * 2.0 * (table_dist - map_dist);
          }
          else if (titer_type == 2) {
            double x = table_dist - map_dist + dilution_stepsize;
            double sig = 1.0 / (1.0 + std::exp(-10.0 * x));
            double d_sig = sig * (1.0 - sig) * 10.0;
            double grad = 2.0 * x * sig + x * x * d_sig;
            colbase_gradients(sr) += weight * grad;
          }
        }
      }
    }

    // ========================================================================
    // UPDATE TABLE DISTANCE MATRIX
    // ========================================================================
    void update_table_dist_matrix() {
      for (arma::uword sr = 0; sr < num_sr; ++sr) {
        for (arma::uword ag = 0; ag < num_ags; ++ag) {
          if (titertype_matrix.at(ag, sr) > 0) {
            tabledist_matrix(ag, sr) = colbases(sr) - logtiter_matrix(ag, sr);
          } else {
            tabledist_matrix(ag, sr) = arma::datum::nan;
          }
        }
      }
    }

    // ========================================================================
    // UPDATE PARAMETERS (coordinates and optionally column bases)
    // ========================================================================
    void update_parameters(const arma::mat &pars) {
      for (arma::uword j = 0; j < num_dims; ++j) {
        for (arma::uword i = 0; i < moveable_ags.n_elem; ++i) {
          ag_coords.at(moveable_ags(i), j) = pars.at(i, j);
        }
      }

      for (arma::uword j = 0; j < num_dims; ++j) {
        for (arma::uword i = 0; i < moveable_sr.n_elem; ++i) {
          sr_coords.at(moveable_sr(i), j) = pars.at(i + moveable_ags.n_elem, j);
        }
      }

      if (optimize_colbases) {
        arma::uword coord_rows = moveable_ags.n_elem + moveable_sr.n_elem;
        for (arma::uword i = 0; i < moveable_colbases.n_elem; ++i) {
          colbases(moveable_colbases(i)) = pars.at(coord_rows + i, 0);
        }
      }
    }

    // ========================================================================
    // CALCULATE STRESS
    // ========================================================================
    double calculate_stress(){
      stress = 0;

      for(sri = included_srs.begin(); sri != sri_end; ++sri) {
        for(agi = included_ags.begin(); agi != agi_end; ++agi) {
          if(titertype_matrix.at(*agi,*sri) <= 0) continue;

          stress += titer_weights.at(*agi,*sri) * ac_ptStress(
            mapdist_matrix.at(*agi,*sri),
            tabledist_matrix.at(*agi,*sri),
            titertype_matrix.at(*agi,*sri),
            dilution_stepsize
          );
        }
      }

      return stress;
    }

    // ========================================================================
    // UPDATE MAP DISTANCE MATRIX
    // ========================================================================
    void update_map_dist_matrix(){
      for(sri = included_srs.begin(); sri != sri_end; ++sri) {
        for(agi = included_ags.begin(); agi != agi_end; ++agi) {
          if(titertype_matrix.at(*agi,*sri) <= 0) continue;

          mapdist_matrix.at(*agi,*sri) = sqrt(arma::accu(arma::square(
            ag_coords.row(*agi) - sr_coords.row(*sri)
          )));
        }
      }
    }

    // ========================================================================
    // GET COLUMN BASES
    // ========================================================================
    arma::vec get_colbases() const {
      return colbases;
    }

};


// ============================================================================
// HELPER FUNCTION: Create initial parameter matrix
// ============================================================================
arma::mat create_initial_pars(
    const arma::mat &ag_coords,
    const arma::mat &sr_coords,
    const arma::uvec &moveable_ags,
    const arma::uvec &moveable_sr,
    const arma::vec &colbases,
    const arma::uvec &moveable_colbases,
    bool optimize_colbases
) {
  arma::uword num_dims = ag_coords.n_cols;
  arma::uword coord_rows = moveable_ags.n_elem + moveable_sr.n_elem;
  arma::uword total_rows = coord_rows;
  
  if (optimize_colbases) {
    total_rows += moveable_colbases.n_elem;
  }
  
  arma::mat pars(total_rows, num_dims, arma::fill::zeros);
  
  for (arma::uword i = 0; i < moveable_ags.n_elem; ++i) {
    pars.row(i) = ag_coords.row(moveable_ags(i));
  }
  
  for (arma::uword i = 0; i < moveable_sr.n_elem; ++i) {
    pars.row(i + moveable_ags.n_elem) = sr_coords.row(moveable_sr(i));
  }
  
  if (optimize_colbases) {
    for (arma::uword i = 0; i < moveable_colbases.n_elem; ++i) {
      pars.at(coord_rows + i, 0) = colbases(moveable_colbases(i));
    }
  }
  
  return pars;
}


// ============================================================================
// EXPORTED FUNCTIONS
// ============================================================================

// [[Rcpp::export]]
double ac_coords_stress(
    const AcTiterTable &titers,
    const std::string &min_colbasis,
    const arma::vec &fixed_colbases,
    const arma::vec &ag_reactivity_adjustments,
    arma::mat &ag_coords,
    arma::mat &sr_coords,
    double dilution_stepsize
){
  int num_dims = ag_coords.n_cols;

  MapOptimizer map(
      ag_coords,
      sr_coords,
      titers.numeric_table_distances(
        min_colbasis,
        fixed_colbases,
        ag_reactivity_adjustments
      ),
      titers.get_titer_types(),
      num_dims,
      dilution_stepsize
  );

  return map.calculate_stress();
}

// [[Rcpp::export]]
arma::mat ac_point_stresses(
    AcTiterTable titer_table,
    std::string min_colbasis,
    arma::vec fixed_colbases,
    arma::vec ag_reactivity_adjustments,
    arma::mat map_dists,
    double dilution_stepsize
){
  arma::uword num_ags = map_dists.n_rows;
  arma::uword num_sr  = map_dists.n_cols;
  arma::mat numeric_table_dists = titer_table.numeric_table_distances(
    min_colbasis,
    fixed_colbases,
    ag_reactivity_adjustments
  );
  arma::imat titer_types = titer_table.get_titer_types();

  arma::mat stress_table(num_ags, num_sr);

  for (arma::uword ag = 0; ag < num_ags; ag++) {
    for (arma::uword sr = 0; sr < num_sr; sr++) {
      if (std::isnan(map_dists(ag, sr))) {
        stress_table(ag, sr) = arma::datum::nan;
      } else {
        stress_table(ag, sr) = ac_ptStress(
          map_dists(ag, sr),
          numeric_table_dists(ag, sr),
          titer_types(ag, sr),
          dilution_stepsize
        );
      }
    }
  }

  return(stress_table);
}


// [[Rcpp::export]]
arma::mat ac_point_residuals(
    const AcMap &map,
    const arma::uword &optimization_number
){
  arma::uword num_ags = map.antigens.size();
  arma::uword num_sr  = map.sera.size();
  arma::mat numeric_table_dists = map.optimizations.at(optimization_number).numeric_table_distances(
    map.titer_table_flat
  );
  arma::imat titer_types = map.titer_table_flat.get_titer_types();
  arma::mat map_dists = map.optimizations.at(optimization_number).distance_matrix();
  double dilution_stepsize = map.dilution_stepsize;

  arma::mat residual_table(num_ags, num_sr);

  for (arma::uword ag = 0; ag < num_ags; ag++) {
    for (arma::uword sr = 0; sr < num_sr; sr++) {
      if (std::isnan(map_dists(ag, sr))) {
        residual_table(ag, sr) = arma::datum::nan;
      } else {
        residual_table(ag, sr) = ac_ptResidual(
          map_dists(ag, sr),
          numeric_table_dists(ag, sr),
          titer_types(ag, sr),
          dilution_stepsize
        );
      }
    }
  }

  return(residual_table);
}


// [[Rcpp::export]]
double ac_relax_coords(
    const arma::mat &tabledist_matrix,
    const arma::imat &titertype_matrix,
    arma::mat &ag_coords,
    arma::mat &sr_coords,
    const AcOptimizerOptions &options,
    const arma::uvec &fixed_antigens,
    const arma::uvec &fixed_sera,
    const arma::mat &titer_weights,
    const double &dilution_stepsize
){
  arma::uvec ag_fixed(ag_coords.n_rows, arma::fill::zeros);
  arma::uvec sr_fixed(sr_coords.n_rows, arma::fill::zeros);
  ag_fixed.elem(fixed_antigens).ones();
  sr_fixed.elem(fixed_sera).ones();
  ag_fixed.elem(arma::find_nonfinite(ag_coords.col(0))).ones();
  sr_fixed.elem(arma::find_nonfinite(sr_coords.col(0))).ones();

  MapOptimizer map(
    ag_coords,
    sr_coords,
    tabledist_matrix,
    titertype_matrix,
    ag_coords.n_cols,
    ag_fixed,
    sr_fixed,
    titer_weights,
    dilution_stepsize
  );

  arma::mat pars = arma::join_cols(
    ag_coords.rows(arma::find(ag_fixed == 0)),
    sr_coords.rows(arma::find(sr_fixed == 0))
  );

  ens::L_BFGS lbfgs(
    options.num_basis,
    options.maxit,
    options.armijo_constant,
    options.wolfe,
    options.min_gradient_norm,
    options.factr,
    options.max_line_search_trials,
    options.min_step,
    options.max_step
  );

  lbfgs.Optimize(map, pars);

  ag_coords = map.ag_coords;
  sr_coords = map.sr_coords;
  return map.calculate_stress();
}


// [[Rcpp::export]]
Rcpp::List ac_relax_coords_with_colbases(
    const arma::mat &logtiter_matrix,
    const arma::imat &titertype_matrix,
    arma::mat &ag_coords,
    arma::mat &sr_coords,
    arma::vec &colbases,
    const AcOptimizerOptions &options,
    const arma::uvec &fixed_antigens,
    const arma::uvec &fixed_sera,
    const arma::mat &titer_weights,
    const double &dilution_stepsize
){
  arma::uword num_dims = ag_coords.n_cols;
  arma::uword num_ags = ag_coords.n_rows;
  arma::uword num_sr = sr_coords.n_rows;
  
  arma::uvec ag_fixed(num_ags, arma::fill::zeros);
  arma::uvec sr_fixed(num_sr, arma::fill::zeros);
  ag_fixed.elem(fixed_antigens).ones();
  sr_fixed.elem(fixed_sera).ones();
  ag_fixed.elem(arma::find_nonfinite(ag_coords.col(0))).ones();
  sr_fixed.elem(arma::find_nonfinite(sr_coords.col(0))).ones();

  MapOptimizer map(
    ag_coords,
    sr_coords,
    logtiter_matrix,
    titertype_matrix,
    colbases,
    num_dims,
    ag_fixed,
    sr_fixed,
    titer_weights,
    dilution_stepsize,
    true
  );

  arma::mat pars = create_initial_pars(
    ag_coords,
    sr_coords,
    map.moveable_ags,
    map.moveable_sr,
    colbases,
    map.moveable_colbases,
    true
  );

  ens::L_BFGS lbfgs(
    options.num_basis,
    options.maxit,
    options.armijo_constant,
    options.wolfe,
    options.min_gradient_norm,
    options.factr,
    options.max_line_search_trials,
    options.min_step,
    options.max_step
  );

  lbfgs.Optimize(map, pars);

  ag_coords = map.ag_coords;
  sr_coords = map.sr_coords;
  colbases = map.get_colbases();
  double final_stress = map.calculate_stress();

  return Rcpp::List::create(
    Rcpp::Named("stress") = final_stress,
    Rcpp::Named("ag_coords") = ag_coords,
    Rcpp::Named("sr_coords") = sr_coords,
    Rcpp::Named("colbases") = colbases
  );
}


// ============================================================================
// GENERATE OPTIMIZATIONS
// ============================================================================

std::vector<AcOptimization> ac_generateOptimizations(
    const arma::mat &tabledist_matrix,
    const arma::imat &titertype_matrix,
    const std::string &min_colbasis,
    const arma::vec &fixed_colbases,
    const arma::vec &ag_reactivity_adjustments,
    const int &num_dims,
    const int &num_optimizations,
    const AcOptimizerOptions &options,
    const double &dilution_stepsize
){
  int num_ags = tabledist_matrix.n_rows;
  int num_sr = tabledist_matrix.n_cols;

  AcOptimization initial_optim = AcOptimization(
    num_dims,
    num_ags,
    num_sr,
    min_colbasis,
    fixed_colbases,
    ag_reactivity_adjustments
  );

  initial_optim.randomizeCoords( tabledist_matrix.max() );
  initial_optim.relax_from_raw_matrices(
    tabledist_matrix,
    titertype_matrix,
    options,
    arma::uvec(),
    arma::uvec(),
    arma::mat(),
    dilution_stepsize
  );

  arma::mat distmat = initial_optim.distance_matrix();
  double coord_maxdist = distmat.max();
  double coord_boxsize = coord_maxdist*2;

  std::vector<AcOptimization> optimizations;
  for(int i=0; i<num_optimizations; i++){
    AcOptimization optimization(
        num_dims,
        num_ags,
        num_sr,
        min_colbasis,
        fixed_colbases,
        ag_reactivity_adjustments
    );
    optimization.randomizeCoords(coord_boxsize);
    optimizations.push_back(optimization);
  }

  return optimizations;
}


// ============================================================================
// RELAX OPTIMIZATIONS (original - fixed column bases)
// ============================================================================

void ac_relaxOptimizations(
  std::vector<AcOptimization>& optimizations,
  arma::uword num_dims,
  const arma::mat &tabledist_matrix,
  const arma::imat &titertype_matrix,
  const AcOptimizerOptions &options,
  const arma::mat &titer_weights,
  const double &dilution_stepsize
){
  int num_optimizations = optimizations.size();

  if(options.report_progress) REprintf("Performing %d optimizations\n", num_optimizations);
  AcProgressBar pb(options.progress_bar_length, options.report_progress);
  Progress p(num_optimizations, true, pb);

  arma::uvec dim_set { num_dims };
  if (options.dim_annealing) {
    dim_set.set_size(2);
    dim_set(0) = 5;
    dim_set(1) = num_dims;
  }

  #pragma omp parallel for schedule(dynamic) num_threads(options.num_cores)
  for (int i=0; i<num_optimizations; i++) {
    if (!p.check_abort()) {
      p.increment();

      for (arma::uword j=0; j<dim_set.n_elem; j++) {
        optimizations.at(i).relax_from_raw_matrices(
            tabledist_matrix,
            titertype_matrix,
            options,
            arma::uvec(),
            arma::uvec(),
            titer_weights,
            dilution_stepsize
        );

        if (dim_set(j) != num_dims) {
          optimizations.at(i).reduceDimensions(dim_set(j + 1));
        }
      }
    }
  }

  if (p.is_aborted()) {
    ac_error("Optimization runs interrupted");
  } else {
    pb.complete("Optimization runs complete");
  }
}


// ============================================================================
// RELAX OPTIMIZATIONS WITH COLUMN BASES
// ============================================================================

void ac_relaxOptimizations_with_colbases(
  std::vector<AcOptimization>& optimizations,
  arma::uword num_dims,
  const arma::mat &logtiter_matrix,
  const arma::imat &titertype_matrix,
  const arma::vec &initial_colbases,
  const AcOptimizerOptions &options,
  const arma::mat &titer_weights,
  const double &dilution_stepsize
){
  int num_optimizations = optimizations.size();

  if(options.report_progress) REprintf("Performing %d optimizations (with colbase optimization)\n", num_optimizations);
  AcProgressBar pb(options.progress_bar_length, options.report_progress);
  Progress p(num_optimizations, true, pb);

  arma::uvec dim_set { num_dims };
  if (options.dim_annealing) {
    dim_set.set_size(2);
    dim_set(0) = 5;
    dim_set(1) = num_dims;
  }

  #pragma omp parallel for schedule(dynamic) num_threads(options.num_cores)
  for (int i = 0; i < num_optimizations; i++) {
    if (!p.check_abort()) {
      p.increment();

      arma::vec opt_colbases = initial_colbases;

      for (arma::uword j = 0; j < dim_set.n_elem; j++) {
        arma::mat ag_coords = optimizations.at(i).agCoords();
        arma::mat sr_coords = optimizations.at(i).srCoords();

        Rcpp::List result = ac_relax_coords_with_colbases(
          logtiter_matrix,
          titertype_matrix,
          ag_coords,
          sr_coords,
          opt_colbases,
          options,
          arma::uvec(),
          arma::uvec(),
          titer_weights,
          dilution_stepsize
        );

        optimizations.at(i).set_ag_base_coords(Rcpp::as<arma::mat>(result["ag_coords"]));
        optimizations.at(i).set_sr_base_coords(Rcpp::as<arma::mat>(result["sr_coords"]));
        optimizations.at(i).set_stress(Rcpp::as<double>(result["stress"]));
        
        opt_colbases = Rcpp::as<arma::vec>(result["colbases"]);
        optimizations.at(i).set_fixed_column_bases(opt_colbases, false);

        if (dim_set(j) != num_dims) {
          optimizations.at(i).reduceDimensions(dim_set(j + 1));
        }
      }
    }
  }

  if (p.is_aborted()) {
    ac_error("Optimization runs interrupted");
  } else {
    pb.complete("Optimization runs complete");
  }
}


// ============================================================================
// RUN OPTIMIZATIONS (main entry point) - ONLY ONE VERSION
// ============================================================================

// [[Rcpp::export]]
std::vector<AcOptimization> ac_runOptimizations(
    const AcTiterTable &titertable,
    const std::string &minimum_col_basis,
    const arma::vec &fixed_colbases,
    const arma::vec &ag_reactivity_adjustments,
    const arma::uword &num_dims,
    const arma::uword &num_optimizations,
    const AcOptimizerOptions &options,
    const arma::mat &titer_weights,
    const double &dilution_stepsize
){
  arma::uword num_ags = titertable.nags();
  arma::uword num_sr = titertable.nsr();

  // Check if we should optimize column bases
  bool optimize_colbases = options.optimize_colbases;

  // Calculate log titer matrix
  arma::mat logtiter_matrix = arma::log2(titertable.get_numeric_titers() / 10.0);
  logtiter_matrix.each_col() += ag_reactivity_adjustments;
  
  arma::imat titertype_matrix = titertable.get_titer_types();

  // Set initial column bases
  arma::vec initial_colbases = titertable.calc_colbases(
    minimum_col_basis,
    fixed_colbases,
    ag_reactivity_adjustments
  );

  // Compute initial table distance matrix
  arma::mat tabledist_matrix(num_ags, num_sr);
  for (arma::uword sr = 0; sr < num_sr; ++sr) {
    for (arma::uword ag = 0; ag < num_ags; ++ag) {
      if (titertype_matrix(ag, sr) > 0) {
        tabledist_matrix(ag, sr) = initial_colbases(sr) - logtiter_matrix(ag, sr);
      } else {
        tabledist_matrix(ag, sr) = arma::datum::nan;
      }
    }
  }

  // Determine starting dimensions
  arma::uword start_dims;
  if (options.dim_annealing && num_dims < 5) {
    start_dims = 5;
  } else {
    start_dims = num_dims;
  }

  // Generate optimizations
  std::vector<AcOptimization> optimizations = ac_generateOptimizations(
    tabledist_matrix,
    titertype_matrix,
    minimum_col_basis,
    fixed_colbases,
    ag_reactivity_adjustments,
    start_dims,
    num_optimizations,
    options,
    dilution_stepsize
  );

  // Relax optimizations
  if (optimize_colbases) {
    ac_relaxOptimizations_with_colbases(
      optimizations,
      num_dims,
      logtiter_matrix,
      titertype_matrix,
      initial_colbases,
      options,
      titer_weights,
      dilution_stepsize
    );
  } else {
    ac_relaxOptimizations(
      optimizations,
      num_dims,
      tabledist_matrix,
      titertype_matrix,
      options,
      titer_weights,
      dilution_stepsize
    );
  }

  // Sort and align
  sort_optimizations_by_stress(optimizations);
  align_optimizations(optimizations);

  return optimizations;
}