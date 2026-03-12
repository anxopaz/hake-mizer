
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~   Hake's MIZER model  ~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rm(list=ls())

library(dplyr) 
library(tidyr)
library(ggplot2)
library(gridExtra)
library(plotly)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(TMB)

source( './scripts/aux_functions.R')


# Biological parameters ----------------------------

load( './input/Bio_Pars.RData')     # './1-Bio_Pars.R' with Biological Parameters


# Fishing Mortality --------------------------

load( './input/Catch.RData')   # './2-Catch.R' with Catch and LFD data


# SSB / Bio ----------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

quantity <- 'biomass' # 'SSB'

ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)

obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm

species_params(bio_pars)$biomass_observed <- obs_q
species_params(bio_pars)$biomass_cutoff <- b_min

bio_pars <- setBevertonHolt( bio_pars,  # Rdd = Rdi * (Rmax/(Rdi+Rmax))
   reproduction_level = 0.001)          # rep_level = Rdd/Rmax (density dependance degree)



# MIZER model --------------------------------

library(mizerEcopath)


## Species params

sp <- bio_pars@species_params |>
  select(species, w_mat, age_mat, w_max, a, b, n, pred_kernel_type, beta, sigma)

# Biomass
quantity <- 'biomass' # 'SSB'
ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)
obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm
sp$biomass_observed <- obs_q
sp$biomass_cutoff <- b_min

# Max size (from landings)
l_max <- 1.001 * max(corLFD$length + 1, na.rm = TRUE)
sp$w_max = lwf(l_max, sp$a, sp$b)

# Mortality exponent (from SS)
sp$d <- -0.1217


## Allometric params

# We now set up a Hake model with power law encounter rate and
# power law mortality rate.
hake_model <- newAllometricParams(sp)
# This is not yet calibrated


## Fishing (1 gear) -------------------------

# We work with only one "Total" gear
gp <- data.frame(
  gear = "Total", species = "Hake", catchability = 1,
  sel_func = "sigmoid_length",
  l50 = 30,
  l25 = 28,
  yield_observed = sum(corLFDs$catch)
)
gear_params(hake_model) <- gp
initial_effort(hake_model) <- 1

# Set initial catchability to get the observed yield
# This will be calibrated properly in `matchCatch()` below.
yield <- getYield(hake_model)
gp$catchability <- gp$yield_observed / yield
gear_params(hake_model) <- gp

# Transform landings data to required shape, adding landings from
# all gears to give total landings
catch_onefleet <- corLFD |> group_by(length) |> summarise(catch = sum(number)) |>
  mutate( dl = 1, species = 'Hake', gear = "Total")


### MatchCatch

# Calibrate to total yield and landings size distribution
hake_model_onefleet <- matchCatch(hake_model, catch = catch_onefleet)
plotYieldVsSize(hake_model_onefleet, x_var = "Length", catch = catch_onefleet)

source('./allometric/new_funs.R')

hake_model_newfun <- matchCatch(hake_model, catch = catch_onefleet)
plotYieldVsSize(hake_model_newfun, x_var = "Length", catch = catch_onefleet)

plotSpectra2( hake_model_onefleet, hake_model_newfun)
hake_model_onefleet@gear_params; hake_model_newfun@gear_params


## Fishing (different gears) -------------------------
# Same for different gears
catch_allfleets <- corLFD |>
  mutate( catch = number, gear = fleet, dl = 1, species = 'Hake') |>
  select( length, catch, dl, species, gear)

gear_names <- unique( catch_allfleets$gear)

gp <- data.frame(
  gear = gear_names, 
  species = "Hake",
  sel_func = ifelse( gear_names %in% c('palangre','vol'), 'sigmoid_length', 'double_sigmoid_length'),
  l50 =       c( 28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
  l25 =       c( 23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
  l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, NA,   58.0, 54.4, NA  ),
  l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, NA,   67.9, 60.8, NA  ),
  yield_observed = corLFDs$catch,
  catchability = corLFDs$catch/sum(corLFDs$catch))

gear_params(hake_model) <- gp

initial_effort(hake_model) <- 1

### MatchCatch

hake_model_fleets <- matchCatch(hake_model, catch = catch_allfleets)
plotYieldVsSizeByGear(hake_model_fleets, catch = catch_allfleets)


## Surveys ---------------------

surveys_sum <- LFD_sum_surv$aver_y
surveys_data <- LFD_list_surv$aver_y


### All in one survey 

onesurvey <- surveys_data |> group_by(length) |> summarise(catch = sum(number)) |>
  mutate( dl = 1, species = 'Hake', gear = "Surveys")

max_length <- onesurvey %>% filter(catch == max(catch)) %>% ungroup()

gear_names <- unique( onesurvey$gear)

gp <- data.frame(
  gear = gear_names, 
  species = "Hake",
  sel_func = rep( 'sigmoid_length', length(gear_names)),
  l50 = max_length$length*0.8,
  l25 = max_length$length*0.7,
  l50_right = rep( NA, length(gear_names)),
  l25_right = rep( NA, length(gear_names)),
  yield_observed = sum(surveys_sum$catch),
  catchability = 1)

gear_params(hake_model) <- gp

initial_effort(hake_model) <- 1


### MatchCatch

hake_model_onesurvey <- matchCatch(hake_model, catch = onesurvey)
plotYieldVsSizeByGear(hake_model_onesurvey, catch = onesurvey)


### Separated surveys

surveys <- surveys_data |>
  mutate( catch = number, gear = fleet, dl = 1, species = 'Hake') |>
  select( length, catch, dl, species, gear)

max_length <- surveys %>% group_by(gear) %>%
  filter(catch == max(catch)) %>% ungroup()

gear_names <- unique( surveys$gear)

gp <- data.frame(
  gear = gear_names, 
  species = "Hake",
  sel_func = rep( 'sigmoid_length', length(gear_names)),
  l50 = max_length$length*0.8,
  l25 = max_length$length*0.7,
  l50_right = rep( NA, length(gear_names)),
  l25_right = rep( NA, length(gear_names)),
  yield_observed = surveys_sum$catch,
  catchability = surveys_sum$catch/sum(surveys_sum$catch))

gear_params(hake_model) <- gp

initial_effort(hake_model) <- 1


### MatchCatch

hake_model_surveys <- matchCatch(hake_model, catch = surveys)
plotYieldVsSizeByGear(hake_model_surveys, catch = surveys)



# Metabolic loss rate and density dependencies -----------

metab_and_dens <- function( model, rep_level = 0.6, feed_level = 0.6, res_level = 1/2){
  
  # Set metabolic loss rate
  
  # So far we ran the model without metabolic loss. If we now introduce this loss, we have to increase the encounter rate to make up for this.
  species_params(model)$ks <- bio_pars@species_params$ks
  ext_encounter(model) <- ext_encounter(model) + metab(model) / species_params(model)$alpha
  
  # Add density dependencies
  
  # We now have a model whose steady state matches the landings data, but we still need to calibrate its sensitivity to changes away from the steady state, for example its sensitivity to changes in fishing. We don't have a good way to choose the following three parameters yet,  so we'll just make up values for now.
  
  # Set reproduction level
  model <- setBevertonHolt(model, reproduction_level = 0.6)
  
  # Set feeding level
  model <- setFeedingLevel(model, 0.6)
  
  # Set resource level
  model <- alignResource(model) |>
    setResourceInteraction(resource_dynamics = "resource_semichemostat")
  resource_level(model) <- 1/2
  
  return(model)
  
}

hake_model_onefleet <- metab_and_dens( hake_model_onefleet, bio_pars)

hake_model_newfun <- metab_and_dens( hake_model_newfun, bio_pars)

hake_model_fleets <- metab_and_dens( hake_model_fleets, bio_pars)

hake_model_onesurvey <- metab_and_dens( hake_model_onesurvey, bio_pars)

hake_model_surveys <- metab_and_dens( hake_model_surveys, bio_pars)



# Cannibalism ----------------------

hake_model <- hake_model_onefleet

# I'll assume below that 14% of the total diet comes from cannibalism. You will get a warning if you try to increase this and we should discuss this.

load( './data/Diet.RData')
pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% aver_y)]); pcann

diet_matrix <- matrix( c(pcann, 1-pcann), ncol = 2,
    dimnames = list( predator = "Hake", prey = c("Hake", "other")))

cannibal_hake_model <- matchDiet(hake_model, diet_matrix)
plotDietX(cannibal_hake_model)


# Steady state ------------------------

sim <- project(hake_model, t_max = 8)
plotBiomass(sim)

sim <- project(hake_model_newfun, t_max = 8)
plotBiomass(sim)

sim <- project(hake_model_fleets, t_max = 8)
plotBiomass(sim)

sim <- project(hake_model_onesurvey, t_max = 8)
plotBiomass(sim)

sim <- project(hake_model_surveys, t_max = 8)
plotBiomass(sim)

sim <- project(cannibal_hake_model, t_max = 8)
plotBiomass(sim)


# Simulations with increased fishing effort ---------------
sim12 <- project(hake_model, effort = 1.2, t_max = 12, t_save = 0.2)
sim_cannibal12 <- project(cannibal_hake_model, effort = 1.2, t_max = 12, t_save = 0.2)


# Compare biomasses --------------------------

bio12 <- melt(getBiomass(sim12))
bio12$Cannibalism <- "Off"
bio_cannibal12 <- melt(getBiomass(sim_cannibal12))
bio_cannibal12$Cannibalism <- "On"
bio <- rbind(bio12, bio_cannibal12)
ggplot(bio) + geom_line(aes(x = time, y = value, colour = Cannibalism)) 



# MSY curve --------------------

yieldf <- plotYieldVsF(hake_model, species = 'Hake', F_range = seq(0, 1.2, 0.02)); yieldf
max(yieldf$data$yield); yieldf$data$F[which.max(yieldf$data$yield)]

yieldf2 <- plotYieldVsF(cannibal_hake_model, species = 'Hake', F_range = seq(0, 1.2, 0.02)); yieldf2
max(yieldf2$data$yield); yieldf2$data$F[which.max(yieldf2$data$yield)]


# Growth vs NatMort (k/M) -------------------------

## Growth ---------------

mizergrowth <- plotGrowthCurves(hake_model)$data
gage <- mizergrowth$Age

ssgrowthf <- lwf( alf(gage, Linf_f, Kvb, al0_f), a, b)
ssgrowthm <- lwf( alf(gage, Linf_m, Kvb, al0_m), a, b)
ssgrowth <- lwf( alf(gage, Linf, Kvb, al0), a, b)

grdf <- data.frame( Age = gage, Mizer = mizergrowth$value, SS_f = ssgrowthf, ss_m = ssgrowthm, ss_mean = ssgrowth)

grdf %>% pivot_longer( cols = Mizer:ss_mean) %>% 
  ggplot( aes(x = Age, y = value, color = name)) + geom_line() + theme_bw()


## NatMort ------------------------

mizerextmort <- getExtMort(hake_model)
sizes <- as.numeric(colnames(mizerextmort))

testnm <- hake_model@species_params$mu_mat*sizes^hake_model@species_params$d
testnm2 <- hake_model@species_params$mu_mat*(sizes/w50)^hake_model@species_params$d

ssfitextmort <- mu0*sizes^(d)

M <- replist$Natural_Mortality_Bmark
M_female <- as.numeric(subset( M, Seas==1 & Settlement==1 & Sex==1)[-(1:4)])
M_male <- as.numeric(subset( M, Seas==1 & Settlement==1 & Sex==2)[-(1:4)])
NatM <- ( M_female + M_male)/2    # Mean of both sexes

weights <- lwf( c(4, alf( 0:15, Linf, Kvb, al0)[-1]),a,b)

mdf <- rbind(
  data.frame( Weight = weights, Sex = 'SS Female', Mortality = M_female),
  data.frame( Weight = weights, Sex = 'SS Male', Mortality = M_male),
  data.frame( Weight = weights, Sex = 'SS Average', Mortality = NatM),
  data.frame( Weight = sizes, Sex = 'SS refit', Mortality = ssfitextmort),
  data.frame( Weight = sizes, Sex = 'MIZER', Mortality = as.numeric(mizerextmort)),
  data.frame( Weight = sizes, Sex = 'MIZER test', Mortality = testnm),
  data.frame( Weight = sizes, Sex = 'MIZER test 2', Mortality = testnm2))

ggplot( mdf, aes( x = Weight, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point() + xlim(NA, max(weights)) + ylim(NA, 2) + 
  theme(legend.title=element_blank())


# Save ----------------------

hake_mizer <- hake_model

save.image( './output/hake_model.RData')







# Time periods --------------------------------

# Loop for time periods

tpdir <- paste0( getwd(), '/output/time_periods/')
dir.create( path = tpdir, showWarnings = TRUE, recursive = TRUE)

for(i in 1:length(years)){
  names(years)[[i]] <- paste0( years[[i]][1],' - ', years[[i]][length(years[[i]])])}

onefleet_mods <- fleets_mods <- onesurvey_mods <- surveys_mods <- cannibal_mods <-
  tp_mods <- list()


for( i in 1:length(years)){

  ny <- names(years)[i]
  vy <- years[[i]]
  ychar <- paste0( vy[1],'-',vy[length(vy)])

  catch_sum <- corLFD_sum[[i]]
  catch_data <- corLFD_list[[i]]
  
  surveys_sum <- LFD_sum_surv[[i]]
  surveys_data <- LFD_list_surv[[i]]
  
  sp <- bio_pars@species_params |>
    select(species, w_mat, age_mat, w_max, a, b, n, pred_kernel_type, beta, sigma)
  
  quantity <- 'biomass' # 'SSB'
  ss_biomass <- sum(assessment[,quantity][assessment$Year %in% vy]*1e6)/length(vy)
  b_min <- lwf(4,a,b)
  sp$biomass_observed <- ss_biomass
  sp$biomass_cutoff <- b_min
  
  l_max <- 1.001 * max(catch_data$length + 1, na.rm = TRUE)
  sp$w_max = lwf(l_max, sp$a, sp$b)

  sp$d <- -0.1217

  ihake_model <- newAllometricParams(sp)
  
  gp <- data.frame( gear = "Total", species = "Hake", catchability = 1,
    sel_func = "sigmoid_length", l50 = 30, l25 = 28,
    yield_observed = sum(catch_sum$catch))
  
  gear_params(ihake_model) <- gp
  initial_effort(ihake_model) <- 1

  yield <- getYield(ihake_model)
  gp$catchability <- gp$yield_observed / yield
  gear_params(ihake_model) <- gp
  
  catch_onefleet <- catch_data |> group_by(length) |> summarise(catch = sum(number)) |>
    mutate( dl = 1, species = 'Hake', gear = "Total")
  
  ihake_model_onefleet <- matchCatch(ihake_model, catch = catch_onefleet)
  plotYieldVsSize(ihake_model_onefleet, x_var = "Length", catch = catch_onefleet)
  
  
  catch_allfleets <- catch_data |>
    mutate( catch = number, gear = fleet, dl = 1, species = 'Hake') |>
    select( length, catch, dl, species, gear)
  
  gear_names <- unique( catch_allfleets$gear)
  
  gp <- data.frame(
    gear = gear_names, 
    species = "Hake",
    sel_func = ifelse( gear_names %in% c('palangre','vol'), 'sigmoid_length', 'double_sigmoid_length'),
    l50 =       c( 28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
    l25 =       c( 23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
    l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, NA,   58.0, 54.4, NA  ),
    l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, NA,   67.9, 60.8, NA  ),
    yield_observed = catch_sum$catch,
    catchability = catch_sum$catch/sum(catch_sum$catch))
  
  gear_params(ihake_model) <- gp
  
  initial_effort(ihake_model) <- 1
  
  ihake_model_fleets <- matchCatch(ihake_model, catch = catch_allfleets)
  plotYieldVsSizeByGear(ihake_model_fleets, catch = catch_allfleets)
  
  onesurvey <- surveys_data |> group_by(length) |> summarise(catch = sum(number)) |>
    mutate( dl = 1, species = 'Hake', gear = "Surveys")
  
  max_length <- onesurvey %>% filter(catch == max(catch)) %>% ungroup()
  
  gear_names <- unique( onesurvey$gear)
  
  gp <- data.frame(
    gear = gear_names, 
    species = "Hake",
    sel_func = rep( 'sigmoid_length', length(gear_names)),
    l50 = max_length$length*0.8,
    l25 = max_length$length*0.7,
    l50_right = rep( NA, length(gear_names)),
    l25_right = rep( NA, length(gear_names)),
    yield_observed = sum(surveys_sum$catch),
    catchability = 1)
  
  gear_params(ihake_model) <- gp
  
  initial_effort(ihake_model) <- 1
  
  ihake_model_onesurvey <- matchCatch(ihake_model, catch = onesurvey)
  plotYieldVsSizeByGear(ihake_model_onesurvey, catch = onesurvey)
  
  surveys <- surveys_data |>
    mutate( catch = number, gear = fleet, dl = 1, species = 'Hake') |>
    select( length, catch, dl, species, gear)
  
  max_length <- surveys %>% group_by(gear) %>%
    filter(catch == max(catch)) %>% ungroup()
  
  gear_names <- unique( surveys$gear)
  
  gp <- data.frame(
    gear = gear_names, 
    species = "Hake",
    sel_func = rep( 'sigmoid_length', length(gear_names)),
    l50 = max_length$length*0.8,
    l25 = max_length$length*0.7,
    l50_right = rep( NA, length(gear_names)),
    l25_right = rep( NA, length(gear_names)),
    yield_observed = surveys_sum$catch,
    catchability = surveys_sum$catch/sum(surveys_sum$catch))
  
  gear_params(ihake_model) <- gp
  
  initial_effort(ihake_model) <- 1
  
  ihake_model_surveys <- matchCatch(ihake_model, catch = surveys)
  plotYieldVsSizeByGear(ihake_model_surveys, catch = surveys)

  ihake_model_onefleet <- metab_and_dens( ihake_model_onefleet)
  ihake_model_fleets <- metab_and_dens( ihake_model_fleets)
  ihake_model_onesurvey <- metab_and_dens( ihake_model_onesurvey)
  ihake_model_surveys <- metab_and_dens( ihake_model_surveys)
  
  ihake_model <- ihake_model_onefleet
  
  pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% vy)])

  diet_matrix <- matrix( c(pcann, 1-pcann), ncol = 2,
     dimnames = list( predator = "Hake", prey = c("Hake", "other")))
  
  cannibal_ihake_model <- matchDiet(ihake_model, diet_matrix)
  plotDietX(cannibal_ihake_model)
  
  sim <- project(ihake_model, t_max = 8); plotBiomass(sim)
  sim <- project(ihake_model_fleets, t_max = 8); plotBiomass(sim)
  sim <- project(ihake_model_onesurvey, t_max = 8); plotBiomass(sim)
  sim <- project(ihake_model_surveys, t_max = 8); plotBiomass(sim)
  sim <- project(cannibal_ihake_model, t_max = 8); plotBiomass(sim)
  
  tp_mods[[ny]] <- list( onefleet = ihake_model, fleets = ihake_model_fleets, 
    onesurvey = ihake_model_onesurvey, surveys = ihake_model_surveys, 
    cannibal = cannibal_ihake_model)
  
  onefleet_mods[[ny]] = ihake_model
  fleets_mods[[ny]] = ihake_model_fleets
  onesurvey_mods[[ny]] = ihake_model_onesurvey
  surveys_mods[[ny]] = ihake_model_surveys 
  cannibal_mods[[ny]] = cannibal_ihake_model
  
}


# tpdf <- NULL
# for(i in names(tp_mods)){ tpdf <- rbind( tpdf, spf( tp_mods[[i]], name = i))}
# 
# ggplot( tpdf %>% filter(w>10), aes( x = w, y = value, color = Model)) +
#   geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
#   labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Period')
# 
# ggsave( paste0( tpdir, 'spectra_comparison.jpg'), width = 12, height = 7)
# ggsave( paste0( tpdir, 'spectra_comparison_ppt.jpg'), width = 6, height = 4)
# 
# ggplot( tpdf %>% filter(w>10), aes( x = w, y = value2, color = Model)) +
#   geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
#   labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Period')
# 
# ggsave( paste0( tpdir, 'spectra_comparison2.jpg'), width = 12, height = 7)
# ggsave( paste0( tpdir, 'spectra_comparison2_ppt.jpg'), width = 6, height = 4)
# 
# 
# tp_table <- matrix( NA, nrow = length(tp_mods), ncol = length(all_pars),
#                     dimnames = list( names(tp_mods), all_pars))
# 
# for(i in names(tp_mods)){ tp_table[i,] <- parsf( tp_mods[[i]])}
# tp_table








# Save all ----------------------

save.image( './output/alldata.RData')




