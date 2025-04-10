### Von Bertalanffy

alf <- function( a, linf, K, al0) { linf * ( 1 - exp( -K * (a - al0)))}
laf <- function( l, linf, K, al0) { al0 - log(1 - l / linf) / K}


### Weight - Length

lwf <- function( l, a, b) { a*l^b}
wlf <- function( w, a, b) { (w/a)^(1/b)}


### Logistic function: f(x) = 1 / (1 + exp(-b (x - c)))

logf <- function(x, s, a50) { 1 / (1 + exp( -s * (x - a50)))}
logf_shift <- function(x, s, a50, shift) { 1 / (1 + exp( -s * (x - a50 - shift)))}


### Line equation

linef <- function( x, L0, L1, M0, M1) { M0 + (M1-M0)*(x-L0)/(L1-L0)}


### Power Low function

powerlow <- lwf


# Double Normal Length

# 'params' vector elements:
# 1: location of ascending peak
# 2: logistic scale width of plateau
# 3: ascending slope log scale value (slope = exp(p[3]))
# 4: descending slope log scale value (slope = exp(p[4]))
# 5: sel value at initial bin;
# 6: sel value at final bin;

double_normal_length <- function( w, params, species_params) {
  # Preemptive checks
  assert_that(is.numeric(w) && is.numeric(params))
  assert_that(params[1] > 0)
  assert_that(params[5]>=0, params[5]<=1, params[6]>=0, params[6]<=1)
  
  # Weight-length conversion
  a <- species_params[["a"]]
  b <- species_params[["b"]]
  if (is.null(a) || is.null(b)) {
    stop("The selectivity function needs the weight-length parameters ", 
         "`a` and `b` to be provided in the species_params data frame.")
  }
  l <- (w/a)^(1/b)
  
  # Extract parameters
  peak1 <- params[1]
  upselex <- exp(params[3])
  downselex <- exp(params[4])
  
  # Compute derived parameters and points
  peak2 <- peak1 + (0.99 * max(l) - peak1) / (1 + exp(-params[2]))
  point1 <- ifelse(params[5] > 0, params[5], NA)
  point2 <- ifelse(params[6] > 0, params[6], NA)
  
  # Precompute initial scaling factors if points are defined
  t1min <- if (!is.na(point1)) exp(-((min(l) - peak1)^2) / upselex) else NA
  t2min <- if (!is.na(point2)) exp(-((max(l) - peak2)^2) / downselex) else NA
  
  # Vectorized computation for asc and dsc selectivity across x
  t1 <- l - peak1
  t2 <- l - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  
  # Ascending and descending selectivity calculations
  asc <- exp(-t1^2 / upselex)
  dsc <- exp(-t2^2 / downselex)
  
  # Scale asc and dsc selectivity if points are defined
  asc_scl <- if (!is.na(point1)) point1 + (1 - point1) * (asc - t1min) / (1 - t1min) else asc
  dsc_scl <- if (!is.na(point2)) 1 + (point2 - 1) * (dsc - 1) / (t2min - 1) else dsc
  
  # Compute final selectivity using vectorized operations
  sel <- asc_scl * (1 - join1) + join1 * (1 - join2 + dsc_scl * join2)
  
  # Plot and return
  plot(l, sel, col = "red", type = "l", ylab = "Selectivity", xlab = "Length")
  return(sel)
}


plot_lfd <- function( params, catch, dl=1) {
  
  heights <- aggregate(number ~ length, data = catch, FUN = sum)$number / sum(catch$number)
  
  plot( NULL, xlim = c(min(catch$length) - dl, max(catch$length) + dl),
        ylim = c(0, max(heights)), xlab = "Length [cm]", ylab = "Density",
        main = "Histogram with Areas Representing Counts")
  
  rect( catch$length - dl/2, 0, catch$length + dl/2, heights, col = "blue", border = "blue")
  # rect( catch$length, 0, catch$length + dl, heights, col = "blue", border = "blue")
  
  lengths <- (params@w / params@species_params$a)^(1/params@species_params$b)
  
  model_catch <- params@initial_n * getFMort(params)
  model_catch <- model_catch  / sum(model_catch * params@dw)
  model_catch <- model_catch * params@species_params$b * params@w / lengths
  
  lines(lengths, model_catch, col = 'red', lwd = 2)
  legend('topright', legend = c('Observed', 'Modelled'), col = c('blue', 'red'), lwd = 2)

}


plot_lfd_gear <- function( model, catch, maxlim = max( catch$number/sum(catch$number))*unieu(catch$fleet)){

  plist <- list()
  gear_names <- unique(catch$fleet)

  for( i in gear_names){ 
    
    iLFD <- catch %>% filter( fleet == i)
    iLFD$catch <- iLFD$number; iLFD$dl <- 1; iLFD$species <- 'Hake'; iLFD$gear <- iLFD$fleet
    
    plist[[i]] <- plotYieldVsSize( model, species="Hake", gear=i, catch=iLFD, x_var="Length", return_data=FALSE) + 
      theme_bw() + theme(legend.position = "none",axis.title.x=element_blank()) + coord_cartesian(ylim = c(0, maxlim))
    
  }
  
  do.call( gridExtra::grid.arrange, c(plist, ncol = 3))
  
}


prefit <- function( model, catch, dl = 1, yield_lambda = 1e7) {
  
  params <- validParams(model)
  sp <- params@species_params
  
  gears <- unique(catch$fleet)
  n_g <- length(gears)
  
  lengths <- unique(catch$length)
  weights <- sp$a * lengths^sp$b
  
  bins <- data.frame( bin_start = lengths-dl/2, bin_end = lengths+dl/2)
  # bins <- data.frame( bin_start = lengths, bin_end = lengths+dl)
  
  counts <- catch %>% pivot_wider( names_from = fleet, values_from = number, values_fill = 0)
  counts <- as.matrix(counts[,gears])
  
  l_bins <- c( bins[,1], bins[nrow(bins),2]) 
  w_bins <- sp$a * l_bins^sp$b
  
  w_bin_widths <- diff(w_bins)
  
  EReproAndGrowth <- approx(w(params), getEReproAndGrowth(params), xout = weights, rule = 2)$y
  repro_prop <- approx(w(params), repro_prop(params), xout = weights, rule = 2)$y
  repro_prop <- repro_prop / max(repro_prop)
  
  data_list <- list( counts = counts, 
                     bin_widths = w_bin_widths, 
                     bin_boundaries = w_bins, 
                     bin_boundary_lengths = l_bins,
                     weight = weights, 
                     blength = lengths, 
                     yield = params@gear_params$yield_observed, 
                     biomass = sp$biomass_observed, 
                     EReproAndGrowth = EReproAndGrowth, 
                     repro_prop = repro_prop, 
                     w_mat = sp$w_mat, 
                     d = sp$d, 
                     yield_lambda = yield_lambda, 
                     n_g = n_g, 
                     M = sp$M, 
                     U = sp$U,
                     min_len = min(l_bins), 
                     max_len = max(l_bins) )
  
  gp <- gear_params(model)
  
  pars <- list(
    logit_l50 = qlogis((gp$l50 - min(l_bins))/(max(l_bins) - min(l_bins))),
    log_ratio_left = log((gp$l50 - gp$l25)/gp$l50),
    log_l50_right_offset = log(pmax(1e-3, gp$l50_right - gp$l50)),
    log_ratio_right = log((gp$l25_right - gp$l50_right)/gp$l50_right),
    log_catchability = log(gp$catchability))
  
  return(list(data_list = data_list, pars = pars))
  
}


update_params <- function( model, pars, lmin, lmax) {
  
  sp <- model@species_params
  gp <- model@gear_params
  
  logit_l50 <- as.numeric(pars[grep('logit_l50', names(pars))])
  log_ratio_left <- as.numeric(pars[grep('log_ratio_left', names(pars))])
  log_l50_right_offset <- as.numeric(pars[grep('log_l50_right_offset', names(pars))])
  log_ratio_right <- as.numeric(pars[grep('log_ratio_right', names(pars))])
  log_catchability <- as.numeric(pars[grep('log_catchability', names(pars))])
  
  l50 <- lmin + (lmax - lmin) * plogis(logit_l50)
  l25 <- l50 * (1 - exp(log_ratio_left))
  l50_right <- l50 + exp(log_l50_right_offset)
  l25_right <- l50_right * (1 + exp(log_ratio_right))
  catchability <- exp(log_catchability)
  
  gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)
  
  gp[,'l50'] <- gp_res$l50
  gp[,'l25'] <- gp_res$l25
  gp[,'l50_right'] <- gp_res$l50_right
  gp[,'l25_right'] <- gp_res$l25_right
  gp[,'catchability'] <- gp_res$catchability
  
  gear_params(model) <- gp
  
  # recalculate the power-law mortality rate
  # sp$M <- pars["M"]
  # ext_mort(params)[] <- sp$M * params@w^sp$d
  
  # Update the steepness of the maturity ogive
  # sp$w_mat25 <- sp$w_mat / 3^(1 / U)
  # model@species_params <- sp
  # model <- setReproduction(model)
  
  # Calculate the new steady state ----
  model <- steadySingleSpecies(model)
  
  # Rescale it to get the observed biomass
  total <- sum(model@initial_n * model@w * model@dw)
  factor <- sp$biomass_observed / total
  model@initial_n <- model@initial_n * factor
  
  return(model)
}
