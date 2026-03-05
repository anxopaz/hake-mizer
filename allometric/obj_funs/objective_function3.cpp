#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type sel_type, Type logit_l50, Type log_ratio_left, Type log_l50_right_offset, Type log_ratio_right,
                              Type log_catchability, vector<Type> l, Type minl, Type maxl) 
{
  Type l50 = minl + (maxl - minl) * invlogit(logit_l50);  
  Type l25 = l50 * (1 - exp(log_ratio_left)); 
  Type sr = l50 - l25;
  if(sr <= 0) sr = Type(1e-6);
  
  Type s1 = l50 * log(Type(3.0)) / sr;
  Type s2 = s1 / l50;
  Type catchability = exp(log_catchability);
  
  vector<Type> F_mort(l.size());
  if(sel_type == 2) {
    Type l50_right = l50 + exp(log_l50_right_offset); 
    Type l25_right = l50_right * (1 + exp(log_ratio_right));
    Type sr_right = l50_right - l25_right;
    if(sr_right <= 0) sr_right = Type(1e-6);
    Type s1_right = l50_right * log(Type(3.0)) / sr_right;
    Type s2_right = s1_right / l50_right;
    
    for (int i = 0; i < l.size(); i++) {
      Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * l(i)))) * 
        (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * l(i))));
      F_mort(i) = catchability * sel;
    }
  } else {
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

template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth, vector<Type> dw)
{
  int size = dw.size();
  vector<Type> N(size);
  N(0) = Type(1.0);
  for (int i = 1; i < size; ++i) {
    Type denom = growth(i) + mort(i) * dw(i);
    if (!CppAD::isfinite(denom) || denom <= Type(0)) {
      Rcout << "Problema en N(" << i << ") denom=" << denom
            << " growth=" << growth(i)
            << " mort=" << mort(i)
            << " dw=" << dw(i) << std::endl;
      denom = Type(1e-6);
    }
    N(i) = N(i - 1) * growth(i - 1) / denom;
    if (!CppAD::isfinite(N(i)) || N(i) < Type(0)) {
      Rcout << "N(" << i << ") no finito: " << N(i) << std::endl;
      N(i) = Type(1e-6);
    }
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
    if (!CppAD::isfinite(catch_per_bin[i])) {
      Rcout << "catch_per_bin[" << i << "] no finito: " << catch_per_bin[i] << std::endl;
      catch_per_bin[i] = Type(0);
    }
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
        vector<Type> probs(num_bins);
        probs.fill(Type(1e-10));
        
        // Construir probs
        for (int k = 0; k < bin_index.size(); k++) {
          int bi = bin_index(k);
          int fj = f_index(k);
          if (bi < 0 || bi >= num_bins || fj < 0 || fj+1 >= catch_per_bin_g.size()) {
            Rcout << "Index out of range: bi=" << bi << " fj=" << fj << std::endl;
            continue;
          }
          Type add = coeff_fj(k)  * catch_per_bin_g(fj)
            + coeff_fj1(k) * catch_per_bin_g(fj+1);
          probs(bi) += add;
        }
        
        // Floor y normalización
        for (int i = 0; i < num_bins; i++) {
          if (!CppAD::isfinite(probs(i)) || probs(i) < Type(0)) {
            Rcout << "Prob no válida: gear=" << g << " bin=" << i << " prob=" << probs(i) << std::endl;
            probs(i) = Type(1e-12);
          }
        }
        Type s = probs.sum();
        if (!CppAD::isfinite(s) || s <= Type(0)) {
          Rcout << "Suma probs inválida: gear=" << g << " s=" << s << std::endl;
          s = Type(1e-6);
        }
        probs /= s;
        
        // Poisson por bin
        Type Ng = counts_g.sum();
        if (!CppAD::isfinite(Ng) || Ng <= Type(0)) {
          Rcout << "Ng inválido: gear=" << g << " Ng=" << Ng << std::endl;
          Ng = Type(1e-6);
        }
        for (int i = 0; i < num_bins; i++) {
          Type mu_i = Ng * probs(i);
          if (!CppAD::isfinite(mu_i) || mu_i < Type(1)) {
            Rcout << "mu_i inválido: gear=" << g << " bin=" << i << " mu=" << mu_i << std::endl;
            mu_i = Type(1);
          }
          Type ll_i = dpois(counts_g(i), mu_i, true);
          if (!CppAD::isfinite(ll_i)) {
            Rcout << "dpois no finito: gear=" << g << " bin=" << i
                  << " count=" << counts_g(i) << " mu=" << mu_i << std::endl;
            ll_i = Type(0.0); // penaliza suavemente
          }
          nll -= ll_i;
        }
        
        // Penalización yield
        Type ratioY = model_yield_g(g) / yield(g);
        if (!CppAD::isfinite(ratioY) || ratioY <= Type(0)) {
          Rcout << "Yield ratio inválido: gear=" << g
                << " model=" << model_yield_g(g)
                << " obs=" << yield(g) << std::endl;
          nll += Type(1e6);
        } else {
          nll += yield_lambda * pow(log(ratioY), Type(2));
        }
      }
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
