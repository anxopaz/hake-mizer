#' Turn on uniform interaction so that preference is solely size-based
#' 
#' Set all entries in the interaction matrix and all resource interactions to 1
#' Then set the search volume so that predation explains as much of the
#' encounter and mortality as possible. The external encounter and mortality
#' rates are reduced so that the total rates are not changed by this function.
#' 
#' @params params MizerParams object
#' @return An updated MizerParams object
setUniformInteraction <- function(params) {
    if (!isNonInteracting(params)) {
        stop("This function should be called with a non-interacting model.")
    }
    if (!isAllometric(params)) {
        stop("This function requires power law encounter and mortality rates.")
    }
    old_encounter <- getEncounter(params)
    old_mort <- getMort(params)
    
    params@interaction[] <- 1
    params@species_params$interaction_resource <- 1
    
    # To get the encounter and mortality from this interaction set external
    # rates to 0 temporarily
    temp_params <- params
    temp_params@ext_encounter[] <- 0
    temp_params@mu_b[] <- 0
    temp_params@initial_effort[] <- 0
    encounter <- getEncounter(temp_params)
    mort <- getMort(temp_params)
    
    # Determine maximum ratio between new rates and external rates for each
    # species so that we can then divide the search volumes by those ratios.
    # We do not have to worry about division by zero because we know that the
    # external rates are given by power laws, hence never zero.
    encounter_ratio <- encounter / params@ext_encounter
    mort_ratio <- mort / params@mu_b
    # Calculate the maximum ratios for each species
    max_encounter_ratio <- apply(encounter_ratio, 1, max)
    max_mort_ratio <- apply(mort_ratio, 1, max)
    max_ratio <- pmax(max_encounter_ratio, max_mort_ratio)
    # Rescale gamma and search volume
    params@species_params$gamma <- 
        params@species_params$gamma / max_ratio
    params@search_vol <- params@search_vol / max_ratio
    # This rescales also the predation encounter and mortality
    encounter <- encounter / max_ratio
    mort <- mort / max_ratio
    
    new_encounter <- getEncounter(params)
    new_mort <- getMort(params)
    
    
    # Reduce external rates accordingly
    params@ext_encounter <- params@ext_encounter - new_encounter + old_encounter
    params@mu_b <- params@mu_b - new_mort + old_mort
    # Rounding errors could lead to negative values
    params@ext_encounter[params@ext_encounter < 0] <- 0
    params@mu_b[params@mu_b < 0] <- 0
    
    return(params)
}

isNonInteracting <- function(params) {
    all(params@interaction == 0) && 
    all(params@species_params$interaction_resource ==0)
}