
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
library(mizerMR)
library(TMB)

source( './scripts/aux_functions.R')


# Biological parameters ----------------------------

load( './input/Bio_Pars.RData')     # './1-Bio_Pars.R' with Biological Parameters


# Fishing Mortality --------------------------

load( './input/Catch.RData')   # './2-Catch.R' with Catch and LFD data


# SSB ----------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

aver_y

ss_biomass <- assessment$SSB[assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)

species_params(bio_pars)$biomass_observed <- sum(ss_biomass)/length(aver_y)
species_params(bio_pars)$biomass_cutoff <- lwf(4,a,b)   # SS smallest size is 4 cm

species_params(bio_pars)$biomass_observed/1e6
species_params(bio_pars)$biomass_cutoff

bio_pars <- setBevertonHolt( bio_pars,        # Rdd = Rdi * (Rmax/(Rdi+Rmax))
                reproduction_level = 0.001)   # rep_level = Rdd/Rmax (density dependance degree)



# MIZER model --------------------------------

## Match for biomass and growth ---------------

hake_model <- bio_pars |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 


## Catch data by fishing gear ----------------

gear_names <- Catch$fleet

### Double sigmoid selectivity initial parameters

gear_params( hake_model) <- data.frame(
  gear = gear_names, species = "Hake", catchability = 1,
  sel_func = "double_sigmoid_length",
  l50 = c(       28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
  l25 = c(       23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
  l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, 35.6, 58.0, 54.4, 27.3),
  l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, 45.9, 67.9, 60.8, 28.4))

### Catch by gear

gear_params( hake_model)$yield_observed <- corLFDs$catch   # == Catch$catch; != LFDs$catch

gear_params( hake_model)$yield_observed/1e6
sum(gear_params( hake_model)$yield_observed/1e6)


### Initial effort

initial_effort( hake_model) <- .25

hake_model <- matchYield( hake_model)
hake_model <- steady( hake_model)


# Fit ------------------- 

pre_obj <- prefit( model = hake_model, catch = LFD, dl = 1, yield_lambda = 1e7)
pre_obj_m <- prefit( model = hake_model, catch = corLFD, dl = 1, yield_lambda = 1e7)

data_list <- pre_obj$data_list
pars <- pre_obj$pars

data_list_m <- pre_obj_m$data_list
pars_m <- pre_obj_m$pars

# save( data_list, pars, file = "./data/prefit.RData")

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
dyn.load(dynlib("./TMB/fit"))

obj <- MakeADFun(data = data_list, parameters = pars, DLL = "fit")
obj_m <- MakeADFun(data = data_list_m, parameters = pars_m, DLL = "fit")


optim_result <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1, eval.max = 10000, iter.max = 10000))
optim_result_m <- nlminb(obj_m$par, obj_m$fn, obj_m$gr, control = list(trace = 1, eval.max = 10000, iter.max = 10000))

# optim_result <- optimx::optimx(obj$par, obj$fn, obj$gr, control = list(trace = 1, maxit = 10000))
# optim_result_m <- optimx::optimx(obj_m$par, obj_m$fn, obj_m$gr, control = list(trace = 1, maxit = 10000))
# hake_model_fitted <- update_params(  hake_model, optim_result[1,], data_list$min_len, data_list$max_len)
# hake_model_fitted_m <- update_params(  hake_model, optim_result_m[1,], data_list$min_len, data_list$max_len)


hake_model_fitted <- update_params( hake_model, optim_result$par, data_list$min_len, data_list$max_len)
hake_model_fitted_m <- update_params( hake_model, optim_result_m$par, data_list_m$min_len, data_list_m$max_len)



# Check ---------------------------------------

plot_lfd( hake_model, LFD)
plot_lfd( hake_model_fitted, LFD)
plot_lfd( hake_model_fitted_m, corLFD)

plot_lfd_gear( hake_model, LFD, 0.17)
plot_lfd_gear( hake_model_fitted, LFD, 0.12)
plot_lfd_gear( hake_model_fitted_m, corLFD, 0.12)

getYield( hake_model)
sum( hake_model_fitted@gear_params$yield_observed)
sum( hake_model_fitted_m@gear_params$yield_observed)
getYield( hake_model_fitted)
getYield( hake_model_fitted_m)

hake_model_fitted@species_params$biomass_observed
getBiomass( hake_model_fitted)
getBiomass( hake_model_fitted_m)



# Steady ------------------

# hake_mizer <- scaleDownBackground( hake_model_fitted, 1/8000000)
# 
# hake_mizer <- hake_mizer |>
#   calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> matchYield() |> steady() |>
#   calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()
# 
# plot_lfd( hake_mizer, LFD)
# plot_lfd_gear( hake_mizer, LFD, 0.17)
# sum( hake_mizer@gear_params$yield_observed)/ getYield( hake_mizer)
# hake_mizer@species_params$biomass_observed/ getBiomass( hake_mizer)



## ## ## ## ## ## ## ##

# catch_lengths<-data.frame(species="Hake",gear=rep(gear_names,each=129),
#                           length=rep(LFDc$length,9),dl=rep(1,9*129),
#                           weight=rep(w(hake_model),9),dw=rep(dw(hake_model),9),
#                           catch=c(LFD$number))
# 
# tuneParams( hake_mizer, catch = catch_lengths)
# 
# hake_model <- readParams("./output/hake_model.rds")
# hake_model <- scaleDownBackground( hake_model, 1/8000000)






# Save ----------------------

save.image( './output/hake_model.RData')

