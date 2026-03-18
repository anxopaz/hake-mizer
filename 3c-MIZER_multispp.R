
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
spp_pars <- pars_spp |> select(species, w_min, w_max, w_mat, beta, sigma, a, b, n,
   age_mat, d)

spp_pars

n_spp <- nrow(spp_pars)


# Multi-species model -----------------

multi_sp <- newAllometricParams(spp_pars)



# Gear parameters --------------

gp_list <- lapply(spp_mods, function(m) m@gear_params)
gp_list[[length(gp_list) + 1]] <- hake_model@gear_params

gp_list

multi_gp <- bind_rows(gp_list)

gear_params(multi_sp) <- multi_gp

initial_effort(multi_sp) <- 1
multi_sp <- steadySingleSpecies(multi_sp)
multi_sp <- setBevertonHolt(multi_sp, reproduction_level = 0.9)

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

# Interaction -------------
multi_sp <- alignResource(multi_sp, w_pp_cutoff = 1)

plotSpectra(multi_sp)

backup <- multi_sp
multi_sp <- backup

sim <- project(multi_sp, t_max = 20)
plotBiomass(sim)

old_mort <- getMort(multi_sp)
old_encounter <- getEncounter(multi_sp)

multi_sp <- setUniformInteraction( multi_sp)

new_mort <- getMort(multi_sp)
new_encounter <- getEncounter(multi_sp)

max((old_encounter - new_encounter)/new_encounter)
max((old_mort - new_mort)/new_mort)

compareParams(multi_sp, backup)

all(multi_sp@ext_encounter >= 0)
all(multi_sp@mu_b >= 0)

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

# Check --------------

species_params(multi_sp)
gear_params(multi_sp)
interaction_matrix(multi_sp)

plotSpectra( multi_sp)


# Project biomass ------------------

sim <- project(multi_sp, t_max = 20)
plotBiomass(sim)




