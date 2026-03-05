prepare_data <- function(params, species = 1, catch,
                         yield_lambda = 1, production_lambda = 1) {
  
  # Validate MizerParams object and extract data for the selected species ----
  params <- validParams(params)
  species <- valid_species_arg(params, species, error_on_empty = TRUE)
  if (length(species) > 1) {
    stop("Only one species can be updated at a time.") 
  }
  sp <- species_params(params)
  sp_select <- sp$species == species
  sps <- sp[sp_select, ]
  gp <- params@gear_params
  gp_select <- gp$species == species
  gps <- gp[gp_select, ]
  # if (nrow(gps) > 1) {
  #   stop("The code currently assumes that there is only a single gear for each species.")
  # }
  # catch <- valid_catch(catch, species)              # Not true now
  
  # Validate catch data frame and extract data for the selected species ----
  if (nrow(catch) == 0) {
    use_counts <- 0
    counts <- numeric(0)
    w_min <- 1
    w_max <- sps$w_max
  }
  else {
    use_counts <- 1
    max_length <- max(catch$length + catch$dl)
    max_weight <- sps$a * max_length^sps$b
    if (max_weight > sps$w_max) {
      stop("For ", species, " you have observed catches of larger weight than the `w_max` that you specified.")
    }
    
    # Fill in missing zero counts ----
    
    # Sort bins
    # catch <- catch[order(catch$length), ]   # better as follows for different gears
    observed_bins <- catch |>
      mutate(bin_start = length, bin_end = length + dl, count = catch) |>
      select(bin_start, bin_end, count, gear) |>
      arrange(gear, bin_start)
    # if (any(observed_bins$bin_end[-nrow(observed_bins)] > 
    #         observed_bins$bin_start[-1])) {
    #   stop("Bins in the catch data must not overlap.")
    # }     # Not true now
    # Add empty bins at either end. This will have an effect only when the
    # catch data is very poor and would be matched by curves that are still
    # large at the end of the observation interval.
    min_length <- (sps$w_min / sps$a) ^ (1 / sps$b)
    observed_bins <- rbind(observed_bins, data.frame(bin_start = min_length,
                                                     bin_end = min(catch$length), count = 0, gear = unique(catch$gear))) # add diff gears
    max_idx <- which.max(catch$length)
    max_length <- catch$length[max_idx] + catch$dl[max_idx]
    l_max <- (sps$w_max / sps$a) ^ (1 / sps$b)
    observed_bins <- rbind(observed_bins, data.frame(bin_start = max_length,
                                                     bin_end = l_max, count = 0, gear = unique(catch$gear)))  # add diff gears
    observed_bins <- observed_bins |> arrange(gear, bin_start)  # better for different gears
    
    # Create a comprehensive set of bin edges covering all observed bins
    bin_edges <- sort(unique(c(observed_bins$bin_start, observed_bins$bin_end)))
    bins <- data.frame( bin_start = head(bin_edges, -1), bin_end   = tail(bin_edges, -1))
    all_combos <- tidyr::expand_grid(gear=unique(observed_bins$gear),bins)
    
    # Define full bins covering the observed range
    full_bins <- all_combos |> 
      left_join(observed_bins, by = c("gear","bin_start","bin_end")) |>
      mutate(count = tidyr::replace_na(count, 0))  # different fill for missing data
    
    # Extract counts, bin boundaries and widths ----
    counts <- full_bins |> tidyr::pivot_wider(names_from = gear, values_from = count, values_fill = 0)
    counts <- as.matrix(counts)[,-c(1:2)]    # for different gears
    l_bin_boundaries <- unique(c(full_bins$bin_start, full_bins$bin_end))
    w_bin_boundaries <- sps$a * l_bin_boundaries^sps$b
    w_bin_widths <- diff(w_bin_boundaries)
    w <- w(params)
    w_min <- max(w[w <= min(w_bin_boundaries)], sps$w_min)
    w_max <- min(w[w >= max(w_bin_boundaries)], sps$w_max)
  }
  
  w <- w(params)
  w_select <- w >= w_min & w <= w_max
  w <- w[w_select]
  dw <- dw(params)[w_select]
  l <- (w/sps$a)^(1/sps$b)
  
  N <- initialN(params)[sp_select, w_select]
  
  # Calculate biomass above cutoff
  biomass_cutoff <- sps$biomass_cutoff
  # Determine the C++ array index for the first weight bin to
  # be included in the biomass calculation.
  if (is.null(biomass_cutoff) || is.na(biomass_cutoff)) {
    biomass_cutoff_idx <- as.integer(0)
  } else {
    biomass_cutoff_idx <- as.integer(sum(w < biomass_cutoff))
  }
  # The cutoff index for R is one more than the C++ index
  biomass <- sum((N * w * dw)[(biomass_cutoff_idx + 1):length(w)])
  growth <- getEGrowth(params)[sp_select, w_select]
  if (use_counts) {
    # Precompute weights for interpolation
    weight_list <- precompute_weights(w_bin_boundaries, w)
  } else {
    weight_list <- list(bin_index = integer(0), f_index = integer(0),
                        coeff_fj = numeric(0), coeff_fj1 = numeric(0))
  }
  
  # The w_mat relevant for calculating mortality is the w just below it
  w_mat_idx <- sum(params@w < sps$w_mat)
  w_mat <- params@w[w_mat_idx]
  
  # If production is not observed
  if (is.null(sps$production_observed) || is.na(sps$production_observed)) {
    production_lambda <- 0
    production <- 0
    if (!use_counts) {
      # Not enough data
      return(NULL)
    }
  } else {
    production <- sps$production_observed
  }
  
  # If yield is not observed
  if (is.null(gps$yield_observed) || any(is.na(gps$yield_observed)) ||
      any(!(gps$yield_observed > 0))) {
    yield <- 0
    yield_lambda <- 0
  } else {
    yield <- gps$yield_observed
  }
  # Prepare data list for TMB ----
  data <- list(
    use_counts = use_counts,
    counts = counts,
    bin_index = weight_list$bin_index,
    f_index = weight_list$f_index,
    sel_func = ifelse(gp$sel_func=='double_sigmoid_length',1,2), # different selectivity functions
    coeff_fj = weight_list$coeff_fj,
    coeff_fj1 = weight_list$coeff_fj1,
    dw = dw,
    w = w,
    l = l,
    # minl = min(l),
    # maxl = max(l),   # already calculated in C++
    yield = yield,
    production = production,
    biomass = biomass,
    biomass_cutoff_idx = biomass_cutoff_idx,
    growth = growth,
    w_mat = w_mat,
    d = sps$d,
    yield_lambda = yield_lambda,
    production_lambda = production_lambda
  )
  
  return(data)
}

#' Precompute weights for integration of density over observed bins
#'
#' We have a set of weight bins with boundaries `w_bin_boundaries` and need
#' an efficient way to integrate a probability density to determine a
#' probability for each bin. The probability density is available at the set
#' of weights given in the vector `w`. We want to use linear interpolation
#' for the density between these values. Because the values in `w` do not
#' align with the values in `w_bin_boundaries` we need to split each bin into
#' segments. In this function we want to
#' precompute the weights with which we need to add up the density values to
#' approximate the integral over each bin.
#'
#' #' @return A list with vectors `bin_index`, `f_index`, `coeff_fj`, and `coeff_fj1`
#'   used for efficient linear interpolation of density values over bins.

precompute_weights <- function(w_bin_boundaries, w) {
  # Precompute overlaps and weights
  num_bins <- length(w_bin_boundaries) - 1
  num_w <- length(w)
  
  # Initialise lists to store precomputed data
  bin_index <- c()
  f_index <- c()
  coeff_fj <- c()
  coeff_fj1 <- c()
  
  for (i in 1:num_bins) { # Loop over bins
    bin_start <- w_bin_boundaries[i]
    bin_end <- w_bin_boundaries[i + 1]
    
    for (j in 1:(num_w - 1)) { # Loop over w
      x0 <- w[j]
      x1 <- w[j + 1]
      
      # Check for overlap
      segment_start <- max(x0, bin_start)
      segment_end <- min(x1, bin_end)
      
      if (segment_start < segment_end) {
        delta_x <- segment_end - segment_start
        dx_j <- x1 - x0
        
        # Calculate interpolation weights
        w0_k <- (segment_start - x0) / dx_j  # Weight at segment_start
        w1_k <- (segment_end - x0) / dx_j    # Weight at segment_end
        
        # Coefficients for f(j) and f(j+1)
        coeff_fj_k = delta_x * (1 - (w0_k + w1_k)/2)
        coeff_fj1_k = delta_x * (w0_k + w1_k)/2
        
        # Store the precomputed data
        bin_index <- c(bin_index, i - 1)  # Zero-based indexing for C++
        f_index <- c(f_index, j - 1)      # Zero-based indexing for C++
        coeff_fj <- c(coeff_fj, coeff_fj_k)
        coeff_fj1 <- c(coeff_fj1, coeff_fj1_k)
      }
    }
  }
  
  return(list(
    bin_index = bin_index,
    f_index = f_index,
    coeff_fj = coeff_fj,
    coeff_fj1 = coeff_fj1
  ))
}

#' Validate and extract catch data for a single species
#'
#' This helper function ensures that the catch data frame contains the required columns,
#' extracts only the rows for the specified species, and checks that a single gear is used
#' (current assumption of the code which may subsequently be overhauled).

#' @export
valid_catch <- function(catch, species) {
  # Allow "catch" as an alternative name to "count"
  if ("catch" %in% names(catch)) {
    catch$count <- catch$catch
  }
  if (!all(c('length', 'dl', 'count') %in% names(catch))) {
    stop("Data frame 'catch' must contain columns 'length', 'dl', and 'count'.")
  }
  # If this contains data for several species, extract the desired species
  if ("species" %in% names(catch)) {
    catch <- catch[catch$species == species, ]
  }
  # if ("gear" %in% names(catch) && length(unique(catch$gear)) > 1) {
  #   stop("The code currently assumes that there is only a single gear for each species.")
  # }   # not true now
  
  return(catch)
}


#' Match the observed catch and yield
#'
#' This function adjusts various model parameters for the selected species so
#' that the model in steady state reproduces the observed catch size
#' distribution, the observed yield and the observed production, if available.
#'
#' Currently this function is implemented only for the case where there is a
#' single gear catching each species.
#'
#' The function sets new values for the following parameters:
#' * `l50`: The size at which the gear selectivity is 50%.
#' * `l25`: The size at which the gear selectivity is 25%.
#' * `l50_right`: The size at which the gear selectivity is 50% (right side for double_sigmoid_length selectivity).
#' * `l25_right`: The size at which the gear selectivity is 25% (right side for double_sigmoid_length selectivity).
#' * `catchability`: The catchability of the gear.
#' * `mu_mat`: The external mortality at maturity.
#'
#' Only the parameters of the selected species are adjusted. The function then
#' recalculates the corresponding rate arrays in the params object. It sets the
#' initial size spectrum to the steady state size spectrum. The total biomass of
#' each species remains unchanged.
#'
#' The function estimates these parameters by minimizing an objective function.
#' The objective function is the negative log likelihood of the observed catch
#' size distribution given the probabilities predicted by the model plus the sum
#' of squares difference between the log of the observed yield and the log of
#' the predicted yield, multiplied by `yield_lambda`, as well as the sum of
#' squares difference between the log of the observed production and the log of
#' the predicted production, multiplied by `production_lambda`.
#'
#' The function deals with missing data in the following way, for each species
#' individually:
#'
#' -  If the observed yield is not available, the function will only match the
#' observed catch size distribution and the observed production.
#'
#' -  If the observed production is not available, the function will only match the
#' observed catch size distribution and the observed yield.
#'
#' -  If the observed catch size distribution is not available, the function will
#' only match the observed yield and the observed production.
#'
#' -  If neither the observed yield nor the observed production are available, the
#' function raises an error.
#'
#' The catch predicted by the model is calculated by integrating the
#' catchability and the gear selectivity over the size distribution of the
#' species. The size distribution itself is shaped by the model through the
#' interplay between growth and mortality. This is why the gear selectivity and
#' catchability need to be adjusted together with the other parameters that
#' shape the size distribution.
#'
#' The objective function is coded in C++ and the TMB package is used to compile
#' the objective function and to create functions for automatically calculating
#' the gradients. These are then passed to the `nlminb` function to minimize the
#' objective function.
#'
#' @param params A MizerParams object
#' @param species The species for which to match the catch. Optional. By default
#'   all target species are selected. A vector of species names, or a numeric
#'   vector with the species indices, or a logical vector indicating for each
#'   species whether it is to be selected (TRUE) or not.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param lambda The slope of the community spectrum. Default is 2.05.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#' @param production_lambda A parameter that controls the strength of the penalty
#'  for deviation from the observed production.
#'
#' @return A MizerParams object with the adjusted external mortality, gear
#'   selectivity, catchability and steady-state spectrum for the selected
#'   species.
#' @family match functions
#' @examples
#' params <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
#' plot_catch(params, species = "Hake", catch = celtic_catch)
#' # The function leaves the biomass of the species unchanged
#' all.equal(getBiomass(params, use_cutoff = TRUE),
#'           getBiomass(celtic_params, use_cutoff = TRUE),
#'           tol = 1e-4)
#' # It also leaves the energy available to an individual for reproduction
#' # and growth unchanged
#' all.equal(getEReproAndGrowth(params),
#'           getEReproAndGrowth(celtic_params))
#' # The initial size spectrum is set to the steady state size spectrum
#' params_steady <- mizerEcopath::steadySingleSpecies(params)
#' all.equal(initialN(params), initialN(params_steady))
#' @export
#' 
matchCatch <- function(params, species = NULL, catch, lambda = 2.05,
                        yield_lambda = 1, production_lambda = 1) 
{
  species <- valid_species_arg(params, species = species, error_on_empty = TRUE)
  params <- validParams(params)
  if (length(species) > 1) {
    for (s in species) {
      params <- matchCatch(params, species = s, catch = catch, 
                           yield_lambda = yield_lambda, production_lambda = production_lambda)
    }
    return(params)
  }
  
  data <- prepare_data(params, species = species, catch,
                       yield_lambda = yield_lambda,
                       production_lambda = production_lambda)
  if (is.null(data)) {
    warning(species, " can not be matched because neither catches nor production are given.")
    return(params)
  }
  
  sp <- species_params(params)
  gp <- gear_params(params)
  sp_select <- sp$species == species
  sps <- sp[sp_select, ]
  gps <- gp[gp$species == species, ]
  
  mat_idx <- sum(params@w < sps$w_mat)
  w_mat <- params@w[mat_idx]
  if (!"mu_mat" %in% names(sps) || is.na(sps$mu_mat)) {
    # determine external mortality at maturity
    mu_mat <- ext_mort(params)[sp_select, mat_idx]
  } else {
    mu_mat <- sps$mu_mat
  }
  
  # Initial parameter estimates
  # differenet parametrization for selectivity (note that C++ fails with NAs)
  initial_params <- list(
    logit_l50 = qlogis((gps$l50 - min(data$l))/(max(data$l) - min(data$l))),
    log_ratio_left = qlogis((gps$l50 - gps$l25)/gps$l50),
    log_l50_right_offset = ifelse( gps$sel_func == 'double_sigmoid_length', 
                                   log(pmax(1e-3, gps$l50_right - gps$l50)), 1),
    log_ratio_right = ifelse( gps$sel_func == 'double_sigmoid_length', 
                              log((gps$l25_right - gps$l50_right)/gps$l50_right), 1),
    mu_mat = mu_mat,
    log_catchability = log(ifelse(gps$catchability <= 0, 1e-8, gps$catchability)))
  
  # Set parameter bounds
  # Mortality is bounded by the requirement that the juvenile spectrum of
  # each species must be less steep than the community spectrum.
  # With g(w) = g w^n and mu(w) = mu w^{n-1}, the exponent of the juvenile
  # spectrum is -mu/g-n. The exponent of the community spectrum is -lambda.
  g_mat <- getEReproAndGrowth(params)[sp_select, mat_idx]
  mu_mat_max <- g_mat / w_mat * (lambda - sps$n)
  lower_bounds <- upper_bounds <- NULL
  lower_bounds <- c(
    rep(-10, length(data$sel_func))*4,    # l50, ratio_left, l50_right_offset, ratio_right
    0.2, rep(-10, length(data$sel_func))) # mu_mat and log_catchability
  
  upper_bounds <- c( rep(10, length(data$sel_func)*4), mu_mat_max, rep(10, length(data$sel_func)))
  
  map <- list()
  if (!data$use_counts) {
      map$l50 <- factor(NA)
      map$ratio <- factor(NA)
      lower_bounds <- lower_bounds[-(1:2)]
      upper_bounds <- upper_bounds[-(1:2)]
  }
  if (data$yield_lambda == 0) {
      map$catchability <- factor(NA)
      lower_bounds <- lower_bounds[-4]
      upper_bounds <- upper_bounds[-4]
  }
  # Prepare the objective function.
  
  # TMB::compile("./allometric/new_objective_function.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
  
  # load('./allometric/prefit.RData')
  
  dyn.load( TMB::dynlib("./allometric/new_objective_function"))
  
  obj <- TMB::MakeADFun(data = data, parameters = initial_params, DLL = "new_objective_function", 
                        silent = FALSE, debug = TRUE)
  
  # Perform the optimization.
  optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                         lower = lower_bounds, upper = upper_bounds,
                         control = list(trace = 0))
  
  # Set model to use the optimal parameters
  w_select <- w(params) %in% data$w
  optimal_params <- update_params(params, species, optim_result$par, data)
  
  return(optimal_params)
}


#' Prepare a TMB Objective Function for Optimising Model Parameters
#'
#' This function returns a list with the data to be passed to the TMB objective
#' function. The data includes the observed catch data, the model parameters,
#' and some precomputed values that are used in the likelihood calculation.
#' The main preprocessing makes sure that we have a comprehensive set of bins
#' that cover the entire size range, even though there will not be observations
#' at all sizes. Missing observations should be interpreted as a 0 count.
#'
#' @param params A MizerParams object
#' @param species The species for which the data is to be prepared. By default
#'   the first species in the model.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#'
#' @return A list with the data to be passed to the TMB objective function. If
#'   there is no catch data for the species, the function returns NULL.
#' @export
update_params <- function(params, species = 1, pars, data) {
  params <- validParams(params)
  sp <- species_params(params)
  
  # Need to add mu_mat column if it does not exist
  # TODO: change once we have introduced a standard species param for this
  mat_idx <- colSums(outer(params@w, sp$w_mat, "<"))
  mu_mat <- ext_mort(params)[cbind(seq_len(nrow(sp)), mat_idx)]
  params <- set_species_param_default(params, "mu_mat", mu_mat)
  
  species <- valid_species_arg(params, species, error_on_empty = TRUE)
  if (length(species) > 1) {
    stop("Only one species can be updated at a time.")
  }
  
  sp_select <- sp$species == species
  sps <- sp[sp_select, ]
  
  gp <- params@gear_params
  gp_select <- gp$species == species
  gps <- gp[gp_select, ]
  # if (nrow(gps) > 1) {
  #   stop("The code currently assumes that there is only a single gear for each species.")
  # }   # not true now
  
  # Update the gear parameters
  gplist <- list()
  gpnames <- c( 'logit_l50', 'log_ratio_left', 'log_l50_right_offset', 'log_ratio_right',
                'log_catchability')
  
  for (i in gpnames) gplist[[i]] <- as.numeric(pars[grep(i, names(pars))])
  
  l50 <- min(data$l) + (max(data$l) - min(data$l)) * plogis(gplist$logit_l50)
  l25 <- l50 * (1 - plogis(gplist$log_ratio_left))
  l50_right <- l50 + exp(gplist$log_l50_right_offset)
  l25_right <- l50_right * (1 + exp(gplist$log_ratio_right))
  catchability <- exp(gplist$log_catchability)
  
  gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)
  
  gps[,'l50'] <- gp_res$l50
  gps[,'l25'] <- gp_res$l25
  gps[,'l50_right'] <- ifelse(gps$sel_func=='double_sigmoid_length', gp_res$l50_right, NA)
  gps[,'l25_right'] <- ifelse(gps$sel_func=='double_sigmoid_length', gp_res$l25_right, NA)
  gps[,'catchability'] <- gp_res$catchability
  
  gear_params(params)[gp_select, ] <- gps
  
  # recalculate the power-law mortality rate
  sps$mu_mat <- pars["mu_mat"]
  # Note that `mu_mat` is the mortality at the w just below w_mat
  mat_idx <- sum(params@w < sps$w_mat)
  w_mat <- params@w[mat_idx]
  ext_mort(params)[sp_select, ] <-
    pars["mu_mat"] * (params@w / w_mat)^sps$d
  
  params@species_params[sp_select, ] <- sps
  params <- setReproduction(params)
  
  # Calculate the new steady state ----
  params <- mizerEcopath::steadySingleSpecies(params, species = species)
  # Rescale it to get the observed biomass
  params <- matchBiomasses(params, species = species)
  # Set the reproduction level to zero
  rl <- numeric(length(species))
  names(rl) <- species
  params <- setBevertonHolt(params, reproduction_level = rl)
  
  return(params)
}


plotYieldVsSizeByGear <- function( model, catch, return_df = FALSE){
  
  params <- validParams(model)
  
  s <- params@species_params['species']
  a <- as.numeric(params@species_params["a"])
  b <- as.numeric(params@species_params["b"])
  w <- params@w
  l <- wlf(w,a,b)
  gears <- unique(catch$gear)
  glength <- length(gears)
  
  df <- NULL
  
  for( i in 1:glength){
    
    igear <- gears[i]
    
    f_mort <- getFMortGear(params)[i,1,]
    
    catch_w <- f_mort * params@initial_n[1,]
    catch_w <- catch_w/sum(catch_w * params@dw)
    catch_l <- catch_w * b * w/l
    
    df <- rbind(df, data.frame( Length=l, Catch_l=catch_l, 
                                Gear=igear, Type="Estimated"))
    
    cind <- which(catch$gear==igear)
    
    len <- catch$length[cind]
    catch_l <- catch$catch[cind]
    catch_l <- catch_l/sum(catch_l)
    
    df <- rbind(df, data.frame(Length=len, Catch_l=catch_l, 
                               Gear=igear, Type = "Observed"))
    
  }
  
  pl <-   ggplot( df, aes(x = Length, y = Catch_l, color = Type, fill = Type)) + 
    geom_bar(data = subset(df, Type != 'Estimated'), stat = "identity", position = "dodge", alpha = 0.6) + 
    geom_line(data = subset(df, Type == 'Estimated'), linewidth = 1) + theme_bw() +
    facet_wrap( ~Gear) + 
    labs( x = "Size [cm]", y = "Normalised number density [1/cm]", color = NULL, fill = NULL)
  
  print(pl)
  
  if(return_df) return(df)
  
}

