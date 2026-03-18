
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~   Multi-species model  ~~~~~~~~~~~~~~ #
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
library(mizerEcopath)
library(TMB)

source( './scripts/aux_functions.R')
source('./allometric/new_funs.R')

load( './output/hake_models.RData')
load( './output/other_spp.RData')


# Species parameters -----------

sp_list <- lapply(spp_mods, function(m) m@species_params)
sp_list[[length(sp_list) + 1]] <- hake_model@species_params

pars_spp <- bind_rows(sp_list) 
spp_pars <- pars_spp |> select(species, w_min, w_max, w_mat, beta, sigma, a, b, n, p, 
   R_max, erepro, age_mat, d, mu_mat)

spp_pars

n_spp <- nrow(spp_pars)

interaction <- matrix(0, n_spp, n_spp, dimnames = list(
  predator = spp_pars$species, prey = spp_pars$species))

interaction

resource_params <- resource_params(hake_model)


# Multi-species model -----------------

multi_sp <- newMultispeciesParams(
  species_params = spp_pars,
  interaction    = interaction,
  min_w  = min(spp_pars$w_min),
  max_w  = max(spp_pars$w_max),
  no_w = 200,
  w_pp_cutoff = 1
)



# Gear parameters --------------

gp_list <- lapply(spp_mods, function(m) m@gear_params)
gp_list[[length(gp_list) + 1]] <- hake_model@gear_params

gp_list

multi_gp <- bind_rows(gp_list)

gear_params(multi_sp) <- multi_gp

initial_effort(multi_sp) <- 1

multi_sp@species_params$interaction_resource <- 0


# Interaction -------------

setUniformInteraction( multi_sp)


# Check --------------

species_params(multi_sp)
gear_params(multi_sp)
interaction_matrix(multi_sp)

plotSpectra( multi_sp)


# Project biomass ------------------

sim <- project(multi_sp, t_max = 20)
plotBiomass(sim)




