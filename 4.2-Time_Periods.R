
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~   Hake's MIZER model by time periods  ~~~~~~ #
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


# SS Stock Assessment data --------------------------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results


# Biological parameters ----------------------------

load( './input/Bio_Pars.RData')


# Fishing Mortality --------------------------

load( './input/Catch.RData')   # './scripts/Natural_Mortality.R' with SS's catch information


# MIZER model -----------------

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
dyn.load( TMB::dynlib("./TMB/fit"))

gear_names <- Catch$fleet

tp_mods <- list()

for( i in 2:length(years)){
  
  ny <- names(years)[i]
  vy <- years[[i]]
  
  ssbio_tp <- sum(assessment$SSB[assessment$Year %in% vy]*1e6)/5
  
  tp_mods[[ny]] <- bio_pars
  
  species_params(tp_mods[[ny]])$biomass_observed <- ssbio_tp
  species_params(tp_mods[[ny]])$biomass_cutoff <- lwf(4,a,b)
  
  tp_mods[[ny]] <- setBevertonHolt( tp_mods[[ny]], reproduction_level = 0.001)
  
  tp_mods[[ny]] <- tp_mods[[ny]] |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 
  
  gear_params( tp_mods[[ny]]) <- data.frame(
    gear = gear_names, species = "Hake", catchability = 1,
    sel_func = "double_sigmoid_length",
    l50 = c(       28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
    l25 = c(       23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
    l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, 35.6, 58.0, 54.4, 27.3),
    l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, 45.9, 67.9, 60.8, 28.4))
  
  gear_params( tp_mods[[ny]])$yield_observed <- corLFD_sum[[ny]]$catch   # == Catch$catch; != LFDs$catch
  
  initial_effort( tp_mods[[ny]]) <- .25
  
  tp_mods[[ny]] <- matchYield( tp_mods[[ny]])
  tp_mods[[ny]] <- steady( tp_mods[[ny]])
  
  pre_obj <- prefit( model =  tp_mods[[ny]], catch = LFD_list[[ny]], dl = 1, yield_lambda = 1e7)
  
  data_list <- pre_obj$data_list
  pars <- pre_obj$pars
  
  obj <- MakeADFun(data = data_list, parameters = pars, DLL = "fit")
  
  
  optim_result <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1, eval.max = 10000, iter.max = 10000))
  
  tp_mods[[ny]] <- update_params(  tp_mods[[ny]], optim_result$par, data_list$min_len, data_list$max_len)

}



plot_lfd( tp_mods[[1]], LFD)

plot_lfd_gear( tp_mods[[1]], LFD, 0.17)

getYield( tp_mods[[1]])
sum( tp_mods[[1]]@gear_params$yield_observed)

tp_mods[[1]]@species_params$biomass_observed
getBiomass( tp_mods[[1]])




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




