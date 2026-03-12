#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type sel_type, Type logit_l50, Type log_ratio_left, Type log_l50_right_offset, Type log_ratio_right,
                              Type log_catchability, vector<Type> l, Type minl, Type maxl) 
{
  Type l50 = minl + (maxl - minl) * invlogit(logit_l50);  
  Type l25 = l50 * (1 - exp(log_ratio_left)); 
  
  Type sr = l50 - l25;
  Type s1 = l50 * log(Type(3.0)) / sr;
  Type s2 = s1 / l50;
  
  if(sr <= 0) sr = Type(1e-6);
  
  Type catchability = exp(log_catchability);
  
  vector<Type> F_mort(l.size());
  
  if(sel_type == 2) {
  
    Type l50_right = l50 + exp(log_l50_right_offset); 
    Type l25_right = l50_right * (1 + exp(log_ratio_right));
    
    Type sr_right = l50_right - l25_right;
    Type s1_right = l50_right * log(Type(3.0)) / sr_right;
    Type s2_right = s1_right / l50_right;
    
    if(sr_right <= 0) sr_right = Type(1e-6);
    
    for (int i = 0; i < l.size(); i++) {
      Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * l(i)))) * 
        (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * l(i))));
      F_mort(i) = catchability * sel;
    }
    
  } else if( sel_type == 1) {
    
    for (int i = 0; i < l.size(); i++) {
      Type sel = Type(1.0) / (Type(1.0) + exp(s1 - s2 * l(i)));
      F_mort(i) = catchability * sel;
    }
  }
  
  
  Rcout << "g=" << sel_type
        << " l50=" << l50
        << " l25=" << l25
        << " sr=" << sr
        << " catch=" << catchability
        << std::endl;
  
  
  
  return F_mort;
}

// Helper function: Steady-state numbers-at-size
template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth,
                         vector<Type> dw)
{
    int size = dw.size();
    vector<Type> N(size);
    N(0) = Type(1.0);
    for (int i = 1; i < size; ++i) {
        Type denominator = growth(i) + mort(i) * dw(i);
        N(i) = N(i - 1) * growth(i - 1) / denominator;
    }

    // Ensure all elements are finite and >= 0
    TMBAD_ASSERT((N.array().isFinite() && (N.array() >= 0)).all());
    
    for (int i = 1; i < size; ++i) {
      Type denominator = growth(i) + mort(i) * dw(i);
      if(!CppAD::isfinite(denominator) || denominator <= 0) {
        Rcout << "Problema en N(" << i << ") | denominator=" << denominator
              << " growth=" << growth(i)
              << " mort=" << mort(i)
              << " dw=" << dw(i) << std::endl;
      }
      N(i) = N(i - 1) * growth(i - 1) / denominator;
    }
    

    return N;
}

template<class Type>
vector<Type> calculate_catch_per_bin(vector<Type> N, vector<Type> F_mort, vector<Type> dw)
{
  
  vector<Type> densities = N * F_mort;
  int num_bins = dw.size();
  vector<Type> catch_per_bin(num_bins);
  for (int i = 0; i < num_bins; ++i) {
    catch_per_bin[i] = dw[i] * densities[i];
  }
  return catch_per_bin;
}

template<class Type>
Type calculate_yield(vector<Type> catch_per_bin, vector<Type> w)
{
  Type model_yield = Type(0.0);
  for (int i = 0; i < catch_per_bin.size(); ++i) {
    model_yield += catch_per_bin[i] * w[i];
  }
  return model_yield;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
    // **Data Section**
    // Introduce an integer flag to indicate whether 'counts' data are present.
    // Calling code can set 'use_counts = 1' if counts exist, or '0' if they do not.
    DATA_INTEGER(use_counts);
    // The next lines read counts data, but if use_counts = 0, we won't use them.
    DATA_MATRIX(counts);     // Observed count data for each bin
    DATA_VECTOR(sel_type);
    DATA_IVECTOR(bin_index); // Bin indices (for overlapping segments)
    DATA_IVECTOR(f_index);   // Function indices (for overlapping segments)
    DATA_VECTOR(coeff_fj);
    DATA_VECTOR(coeff_fj1);

    DATA_VECTOR(dw);
    DATA_VECTOR(w);
    DATA_VECTOR(l);                 // lengths corresponding to w
    DATA_SCALAR(minl);
    DATA_SCALAR(maxl);
    DATA_VECTOR(yield);             // Observed yield
    DATA_SCALAR(production);        // Observed production
    DATA_SCALAR(biomass);           // Observed biomass
    DATA_INTEGER(biomass_cutoff_idx);  // Index for biomass cutoff weight
    DATA_VECTOR(growth);            // Growth rate
    DATA_SCALAR(w_mat);             // Maturity size (weight)
    DATA_SCALAR(d);                 // Exponent for mortality power-law
    DATA_SCALAR(yield_lambda);      // Penalty strength for yield deviation
    DATA_SCALAR(production_lambda); // Penalty strength for production deviation

    // **Parameter Section**
    PARAMETER_VECTOR(logit_l50);
    PARAMETER_VECTOR(log_ratio_left);
    PARAMETER_VECTOR(log_l50_right_offset);
    PARAMETER_VECTOR(log_ratio_right);
    PARAMETER(mu_mat);  
    PARAMETER_VECTOR(log_catchability); 
    
    int n_bins = l.size();
    int n_g = sel_type.size();
    
    vector<Type> total_F_mort(n_bins);
    matrix<Type> F_mort_mat(n_bins, n_g);
    
    // **Calculate fishing mortality rate**
    total_F_mort.setZero();
    for (int g = 0; g < n_g; ++g) {
      vector<Type> F_mort_g = calculate_F_mort(sel_type[g], logit_l50[g], log_ratio_left[g], log_l50_right_offset[g], 
          log_ratio_right[g], log_catchability[g], l, minl, maxl);
      for (int i = 0; i < n_bins; ++i) {
        F_mort_mat(i, g) = F_mort_g(i);
        total_F_mort(i) += F_mort_g(i);
      }
    }

    // **Calculate total mortality rate**
    vector<Type> mort = mu_mat * pow(w / w_mat, d) + total_F_mort;

    // **Calculate steady-state number density**
    //   This is unscaled (N is not matched to 'biomass' yet).
    vector<Type> N = calculate_N(mort, growth, dw);

    // Rescale to match observed biomass
    Type unscaled_biomass = 0;
    for (int i = biomass_cutoff_idx; i < N.size(); ++i) {
        unscaled_biomass += N(i) * w(i) * dw(i);
    }
    N *= (biomass / unscaled_biomass);  // Now total biomass matches 'biomass'.

    matrix<Type> catch_per_bin_mat(n_bins, n_g);
    vector<Type> model_yield_g(n_g);
    
    for (int g = 0; g < n_g; ++g) {
      vector<Type> F_mort_g(n_bins);
      for (int i = 0; i < n_bins; i++) {
        F_mort_g(i) = F_mort_mat(i, g); 
      }
      
      vector<Type> catch_per_bin_g = calculate_catch_per_bin(N, F_mort_g, dw);
      catch_per_bin_mat.col(g) = catch_per_bin_g;
      model_yield_g(g) = calculate_yield(catch_per_bin_g, w);
    }
    
    // **Negative Log-Likelihood (NLL)**
    Type nll = Type(0.0);
    
    if (use_counts == 1) {
      for (int g = 0; g < n_g; ++g) {
        vector<Type> catch_per_bin_g = catch_per_bin_mat.col(g);
        vector<Type> counts_g = counts.col(g);
        
        int num_bins = counts_g.size();
        int num_segs = bin_index.size();
        vector<Type> probs(num_bins);
        probs.fill(Type(1e-10));
        
        // Accumulate the bin probabilities from overlapping segments
        for (int k = 0; k < num_segs; k++) {
          int i = bin_index(k);
          int j = f_index(k);
          probs(i) += coeff_fj(k)  * catch_per_bin_g(j)
            + coeff_fj1(k) * catch_per_bin_g(j+1);
        }
        
        // Normalize probabilities
        probs /= probs.sum();
        
        // Multinomial likelihood por gear
        nll -= dmultinom(counts_g, probs, true);
        
        // Penalization
        nll += yield_lambda * pow(log(model_yield_g(g) / yield(g)), Type(2));
        
        
        if(model_yield_g(g) <= 0 || yield(g) <= 0) {
          Rcout << "Yield problem: model=" << model_yield_g(g) << " obs=" << yield(g) << std::endl;
        }
        
        
      }
    } // else: skip all calculations involving 'counts'

    if (production_lambda > 0) {
        // **Calculate production**
        vector<Type> production_per_bin = N * growth * dw;
        Type model_production = production_per_bin.sum();
        REPORT(model_production);
        // **Add penalty for deviation from observed production**
        nll += production_lambda * pow(log(model_production / production), Type(2));
    }

    // Final sanity checks
    TMBAD_ASSERT(CppAD::isfinite(nll));
    if (!CppAD::isfinite(nll)) {
        error("nll is not finite");
    }

    // **Reporting**

    REPORT(N);
    REPORT(total_F_mort);
    REPORT(F_mort_mat);
    REPORT(mort);

    // Check final biomass again
    Type total_biomass = 0;
    for (int i = biomass_cutoff_idx; i < N.size(); ++i) {
        total_biomass += N(i) * w(i) * dw(i);
    }
    REPORT(total_biomass);

    return nll;
}
