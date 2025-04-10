
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

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
dyn.load( TMB::dynlib("./TMB/fit"))

# Load mizer model output

load( './output/hake_model.RData')   # Results from "3-MIZER.R"


# Loop for tim periods

tp_mods <- list()

for( i in 2:length(years)){
  
  ny <- names(years)[i]
  vy <- years[[i]]
  
  ssbio_tp <- sum(assessment$SSB[assessment$Year %in% vy]*1e6)/length(vy)
  
  tp_mods[[ny]] <- hake_model_fitted_m
  
  species_params(tp_mods[[ny]])$biomass_observed <- ssbio_tp
  species_params(tp_mods[[ny]])$biomass_cutoff <- lwf(4,a,b)
  
  tp_mods[[ny]] <- setBevertonHolt( tp_mods[[ny]], reproduction_level = 0.001)
  
  tp_mods[[ny]] <- tp_mods[[ny]] |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 
  
  gear_params( tp_mods[[ny]])$yield_observed <- corLFD_sum[[ny]]$catch   # == Catch$catch; != LFDs$catch
  
  initial_effort( tp_mods[[ny]]) <- .25
  
  # tp_mods[[ny]] <- matchYield( tp_mods[[ny]])
  # tp_mods[[ny]] <- steady( tp_mods[[ny]])
  
  ipre_obj <- prefit( model =  tp_mods[[ny]], catch = corLFD_list[[ny]], dl = 1, yield_lambda = 1e7)
  
  idata_list <- ipre_obj$data_list
  ipars <- ipre_obj$pars
  
  iobj <- MakeADFun(data = idata_list, parameters = ipars, DLL = "fit")
  
  ioptim_result <- nlminb(iobj$par, iobj$fn, iobj$gr, control = list(trace = 1, eval.max = 10000, iter.max = 10000))
  
  tp_mods[[ny]] <- update_params(  tp_mods[[ny]], ioptim_result$par, data_list$min_len, data_list$max_len)
  
  # ioptim_result <- optimx::optimx(iobj$par, iobj$fn, iobj$gr, control = list(trace = 1, maxit = 10000))
  # tp_mods[[ny]] <- update_params(  tp_mods[[ny]], ioptim_result[2,], data_list$min_len, data_list$max_len)
  
  tp_mods[[ny]] <- scaleDownBackground( tp_mods[[1]], 1/8000000)

  # tp_mods[[ny]] <- tp_mods[[ny]] |>
  #   calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
  #   calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()
  
  plot_lfd_gear( tp_mods[[ny]], corLFD_list[[ny]], 0.17)
  
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




