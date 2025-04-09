// This is the code that defines the objective function, returning the negative log-likelihood value.
// This function is still a work in progress, as the R session crashes when calling from 
// prepare_TMB_objective_function2.R; the crash is produced by the definition of the 'calculate_growth'
// function in line 52, even if the obj. function (line 186) call to calculate_growth (and everything
// else) is commented.

#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type logit_l50, Type log_ratio_left, 
                              Type log_l50_right_offset, Type log_ratio_right,
                              Type log_catchability, vector<Type> blength, 
                              Type min_len, Type max_len) {
  
  Type l50 = min_len + (max_len - min_len) * invlogit(logit_l50);  
  Type l25 = l50 * (1 - exp(log_ratio_left)); 
  
  Type l50_right = l50 + exp(log_l50_right_offset); 
  Type l25_right = l50_right * (1 + exp(log_ratio_right));
  
  Type catchability = exp(log_catchability);
  
  vector<Type> F_mort(blength.size());
  
  Type sr = l50 - l25;
  Type s1 = l50 * log(Type(3.0)) / sr;
  Type s2 = s1 / l50;
  
  Type sr_right = l50_right - l25_right;
  Type s1_right = l50_right * log(Type(3.0)) / sr_right;
  Type s2_right = s1_right / l50_right;
  
  for (int i = 0; i < blength.size(); i++) {
    Type ilength = blength(i);
    Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * ilength))) * 
      (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * ilength)));
    F_mort(i) = catchability * sel;
  }
  
  return F_mort;
}


template<class Type>
vector<Type> calculate_mort(vector<Type> total_F_mort, Type M, Type d, vector<Type> weight)
{
  vector<Type> mort = M * pow(weight, d) + total_F_mort;
  return mort;
}


template<class Type>
vector<Type> calculate_growth(vector<Type> EReproAndGrowth, vector<Type> repro_prop,
                              Type w_mat, Type U, vector<Type> weight)
{
  Type c1 = Type(1.0);
  vector<Type> psi = repro_prop / (c1 + pow(weight / w_mat, -U));
  vector<Type> growth = EReproAndGrowth * (c1 - psi);
  return growth;
}


template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth,
                         Type biomass,
                         vector<Type> bin_widths,
                         vector<Type> weight)
{
  int size = bin_widths.size();
  vector<Type> N(size);
  N(0) = Type(1.0);
  for (int i = 1; i < size; ++i) {
    Type denominator = growth(i) + mort(i) * bin_widths(i);
    N(i) = N(i - 1) * growth(i - 1) / denominator;
  }
  Type total_biomass = Type(0.0);
  for (int i = 0; i < size; ++i) {
    total_biomass += N(i) * bin_widths(i) * weight(i);
  }
  N = N * biomass / total_biomass;
  
  return N;
}

template<class Type>
vector<Type> calculate_catch_per_bin(vector<Type> N, vector<Type> F_mort, vector<Type> bin_widths)
{

  vector<Type> densities = N * F_mort;
  int num_bins = bin_widths.size(); // Number of bins
  vector<Type> catch_per_bin(num_bins);
  for (int i = 0; i < num_bins; ++i) {
    catch_per_bin[i] = bin_widths[i] * densities[i];
  }
  return catch_per_bin;
}

template<class Type>
Type calculate_yield(vector<Type> catch_per_bin,
                     vector<Type> weight)
{
  Type model_yield = Type(0.0);
  for (int i = 0; i < catch_per_bin.size(); ++i) {
    model_yield += catch_per_bin[i] * weight[i];
  }
  return model_yield;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // **Data Section**
  DATA_MATRIX(counts);               // Counts per gear, dimensions: n_bins x n_g
  DATA_VECTOR(bin_widths);           // Width of each bin in grams
  DATA_VECTOR(bin_boundaries);       // Boundaries of each bin in grams
  DATA_VECTOR(weight);               // Mid of each bin in grams
  DATA_VECTOR(blength);              // Mid of each bin in grams
  DATA_VECTOR(bin_boundary_lengths); // Boundaries of each bin in cm
  DATA_VECTOR(yield);                // Observed yield
  DATA_SCALAR(biomass);              // Observed biomass
  DATA_VECTOR(EReproAndGrowth);      // The rate at which energy is available for growth and reproduction
  DATA_VECTOR(repro_prop);           // Proportion of energy allocated to reproduction
  DATA_SCALAR(w_mat);
  DATA_SCALAR(d);                    // Exponent of mortality power-law
  DATA_SCALAR(yield_lambda);         // controls the strength of the penalty for deviation from the observed yield.
  DATA_INTEGER(n_g);                 // Number of gears
  DATA_SCALAR(M);
  DATA_SCALAR(U);
  DATA_SCALAR(min_len);
  DATA_SCALAR(max_len);
  
  // **Parameter Section**
  PARAMETER_VECTOR(logit_l50);            // Length n_g
  PARAMETER_VECTOR(log_ratio_left);       // Length n_g
  PARAMETER_VECTOR(log_l50_right_offset); // Length n_g
  PARAMETER_VECTOR(log_ratio_right);      // Length n_g
  PARAMETER_VECTOR(log_catchability);     // Length n_g
  
  int n_bins = bin_widths.size();
  int n_bin_boundaries = bin_boundaries.size();
  
  vector<Type> total_F_mort(n_bins);
  total_F_mort.setZero(); // initialize to zero
  matrix<Type> F_mort_mat(n_bins, n_g);
  for (int g = 0; g < n_g; ++g) {
    vector<Type> F_mort_g = calculate_F_mort(
      logit_l50[g], log_ratio_left[g], log_l50_right_offset[g], log_ratio_right[g], log_catchability[g], blength, min_len, max_len);
    for (int i = 0; i < n_bins; ++i) {
      F_mort_mat(i, g) = F_mort_g(i); // Store each gear's F_mort in matrix for later use in 'calculate_catch_per_bin'
      total_F_mort(i) += F_mort_g(i); // Total F_mort is added when calculating overall mortality in 'calculate_mort'
    }
  }
  
  vector<Type> mort = calculate_mort(total_F_mort, M, d, weight);
  
  vector<Type> growth = calculate_growth(EReproAndGrowth, repro_prop, w_mat, U, weight);
  
  vector<Type> N = calculate_N(mort, growth, biomass, bin_widths, weight);
  
  matrix<Type> catch_per_bin_mat(n_bins, n_g);

  for (int g = 0; g < n_g; ++g) {
    
    Eigen::Matrix<Type, Eigen::Dynamic, 1> F_mort_column = F_mort_mat.col(g);
    
    vector<Type> F_mort_g(n_bins);
    for (int i = 0; i < n_bins; i++) {
      F_mort_g(i) = F_mort_column(i);
    }
    
    vector<Type> catch_per_bin_g = calculate_catch_per_bin(N, F_mort_g, bin_widths);
    catch_per_bin_mat.col(g) = catch_per_bin_g;
  }
  
  vector<Type> model_yield_g(n_g);

  for (int g = 0; g < n_g; ++g) {
    
    Eigen::Matrix<Type, Eigen::Dynamic, 1> catch_per_bin_column = catch_per_bin_mat.col(g);

    vector<Type> catch_per_bin_g(n_bins);
    for (int i = 0; i < n_bins; i++) {
      catch_per_bin_g(i) = catch_per_bin_column(i);
    }

    Type yield_g = calculate_yield(catch_per_bin_g, weight);
    model_yield_g(g) = yield_g;
  }
  
  Type nll = Type(0.0);
  for (int g = 0; g < n_g; ++g) {
    vector<Type> catch_per_bin_g = catch_per_bin_mat.col(g);
    vector<Type> counts_g = counts.col(g);
    
    vector<Type> probs_g = catch_per_bin_g + Type(1e-10);
    probs_g = probs_g / probs_g.sum();
    
    Type nll_g = -dmultinom(counts_g, probs_g, true);
    
    nll_g += yield_lambda * pow(log(model_yield_g(g) / yield(g)), Type(2));

    nll += nll_g;
  }
  
  TMBAD_ASSERT(nll >= 0);
  TMBAD_ASSERT(CppAD::isfinite(nll));
  if (!CppAD::isfinite(nll)) error("nll is not finite");
  
  return nll;
}
