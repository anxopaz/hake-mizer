prepare_data <- function (params, species = 1, catch, yield_lambda = 1, production_lambda = 1) 
{
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
  if (nrow(gps) > 1) {
    stop("The code currently assumes that there is only a single gear for each species.")
  }
  catch <- valid_catch(catch, species)
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
    catch <- catch[order(catch$length), ]
    observed_bins <- data.frame(bin_start = catch$length, 
                                bin_end = catch$length + catch$dl, count = catch$count)
    if (any(observed_bins$bin_end[-nrow(observed_bins)] > 
            observed_bins$bin_start[-1])) {
      stop("Bins in the catch data must not overlap.")
    }
    min_length <- (sps$w_min/sps$a)^(1/sps$b)
    observed_bins <- rbind(observed_bins, data.frame(bin_start = min_length, 
                                                     bin_end = min(catch$length), count = 0))
    max_idx <- which.max(catch$length)
    max_length <- catch$length[max_idx] + catch$dl[max_idx]
    l_max <- (sps$w_max/sps$a)^(1/sps$b)
    observed_bins <- rbind(observed_bins, data.frame(bin_start = max_length, 
                                                     bin_end = l_max, count = 0))
    bin_edges <- sort(unique(c(observed_bins$bin_start, observed_bins$bin_end)))
    full_bins <- data.frame(bin_start = bin_edges[-length(bin_edges)], 
                            bin_end = bin_edges[-1], count = 0)
    bin_key <- paste0(full_bins$bin_start, "_", full_bins$bin_end)
    observed_bin_key <- paste0(observed_bins$bin_start, "_", 
                               observed_bins$bin_end)
    matched_indices <- match(observed_bin_key, bin_key)
    full_bins$count[matched_indices] <- observed_bins$count
    counts <- full_bins$count
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
  biomass_cutoff <- sps$biomass_cutoff
  if (is.null(biomass_cutoff) || is.na(biomass_cutoff)) {
    biomass_cutoff_idx <- as.integer(0)
  }
  else {
    biomass_cutoff_idx <- as.integer(sum(w < biomass_cutoff))
  }
  biomass <- sum((N * w * dw)[(biomass_cutoff_idx + 1):length(w)])
  growth <- getEGrowth(params)[sp_select, w_select]
  if (use_counts) {
    weight_list <- precompute_weights(w_bin_boundaries, w)
  }
  else {
    weight_list <- list(bin_index = integer(0), f_index = integer(0), 
                        coeff_fj = numeric(0), coeff_fj1 = numeric(0))
  }
  w_mat_idx <- sum(params@w < sps$w_mat)
  w_mat <- params@w[w_mat_idx]
  if (is.null(sps$production_observed) || is.na(sps$production_observed)) {
    production_lambda <- 0
    production <- 0
    if (!use_counts) {
      return(NULL)
    }
  }
  else {
    production <- sps$production_observed
  }
  if (is.null(gps$yield_observed) || is.na(gps$yield_observed) || 
      !(gps$yield_observed > 0)) {
    yield <- 0
    yield_lambda <- 0
  }
  else {
    yield <- gps$yield_observed
  }
  data <- list(use_counts = use_counts, counts = counts, bin_index = weight_list$bin_index, 
               f_index = weight_list$f_index, coeff_fj = weight_list$coeff_fj, 
               coeff_fj1 = weight_list$coeff_fj1, dw = dw, w = w, l = l, 
               yield = yield, production = production, biomass = biomass, 
               biomass_cutoff_idx = biomass_cutoff_idx, growth = growth, 
               w_mat = w_mat, d = sps$d, yield_lambda = yield_lambda, 
               production_lambda = production_lambda)
  return(data)
}


update_data <- function (params, species = 1, pars, data) 
{
  params <- validParams(params)
  sp <- species_params(params)
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
  if (nrow(gps) > 1) {
    stop("The code currently assumes that there is only a single gear for each species.")
  }
  if (data$use_counts) {
    gps$l50 <- pars["l50"]
    gps$l25 <- pars["ratio"] * pars["l50"]
  }
  if (data$yield_lambda > 0) {
    gps$catchability <- pars["catchability"]
  }
  gear_params(params)[gp_select, ] <- gps
  sps$mu_mat <- pars["mu_mat"]
  mat_idx <- sum(params@w < sps$w_mat)
  w_mat <- params@w[mat_idx]
  ext_mort(params)[sp_select, ] <- pars["mu_mat"] * (params@w/w_mat)^sps$d
  params@species_params[sp_select, ] <- sps
  params <- setReproduction(params)
  params <- mizerEcopath::steadySingleSpecies(params, species = species)
  params <- matchBiomasses(params, species = species)
  rl <- numeric(length(species))
  names(rl) <- species
  params <- setBevertonHolt(params, reproduction_level = rl)
  return(params)
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
matchCatch <- function(params, species = NULL, catch, lambda = 2.05,
                       yield_lambda = 1, production_lambda = 1) {
  species <- valid_species_arg(params, species = species,
                               error_on_empty = TRUE)
  params <- validParams(params)
  if (length(species) > 1) {
    for (s in species) {
      params <- matchCatch(params, species = s, catch = catch,
                           yield_lambda = yield_lambda,
                           production_lambda = production_lambda)
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
  initial_params <- c(l50 = gps$l50, ratio = gps$l25 / gps$l50,
                      mu_mat = mu_mat,
                      # we need non-zero catchability to match catch
                      catchability = max(gps$catchability, 1e-8))
  
  # Set parameter bounds
  # Mortality is bounded by the requirement that the juvenile spectrum of
  # each species must be less steep than the community spectrum.
  # With g(w) = g w^n and mu(w) = mu w^{n-1}, the exponent of the juvenile
  # spectrum is -mu/g-n. The exponent of the community spectrum is -lambda.
  g_mat <- getEReproAndGrowth(params)[sp_select, mat_idx]
  mu_mat_max <- g_mat / w_mat * (lambda - sps$n)
  lower_bounds <- c(l50 = 5, ratio = 0.1, mu_mat = 0.2,
                    catchability = 1e-8)
  upper_bounds <- c(l50 = Inf, ratio = 0.99, mu_mat = mu_mat_max,
                    catchability = Inf)
  
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
  obj <- MakeADFun(data = data,
                   parameters = initial_params,
                   map = map,
                   DLL = "mizerEcopath",
                   silent = TRUE)
  
  # Perform the optimization.
  optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                         lower = lower_bounds, upper = upper_bounds,
                         control = list(trace = 0))
  
  # Set model to use the optimal parameters
  w_select <- w(params) %in% data$w
  optimal_params <- update_params(params, species, optim_result$par, data)
  
  return(optimal_params)
}