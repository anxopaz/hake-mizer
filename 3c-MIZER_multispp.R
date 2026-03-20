
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

load( './output/hake_model.RData')
load( './output/other_spp.RData')


# Species parameters -----------

for(i in 1:length(spp_mods)) print(getBiomass(spp_mods[[i]])/spp_mods[[i]]@species_params$biomass_observed)
getBiomass(hake_model)/hake_model@species_params$biomass_observed

for(i in 1:length(spp_mods)) print(getBiomass(spp_mods[[i]]))
getBiomass(hake_model)

sp_list <- lapply(spp_mods, function(m) m@species_params)
sp_list[[length(sp_list) + 1]] <- hake_model@species_params

pars_spp <- bind_rows(sp_list) 
spp_pars <- pars_spp |> select(species, w_min, w_max, w_mat, beta, sigma, a, b, n,
   age_mat, d)

spp_pars

n_spp <- nrow(spp_pars)


# Other option -----------

allcatch <- c( list( Hake = hake_catch), catch_mods)
allmods <- c( Hake = hake_model, spp_mods)

catchdf <- bind_rows( allcatch) 

save( allmods, file='./mods.RData')

# msmod <- bindParams( allmods)


# Multi-species model -----------------

multi_sp <- newAllometricParams(spp_pars)


## Gear parameters --------------

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


getBiomass(multi_sp)

for(i in names(allmods)){
  print(paste(i,round(allmods[[i]]@species_params$biomass_observed/10^6)))
}

for(i in names(allmods)){
  print(plotBiomass(project(allmods[[i]], t_max = 10)))
}

## Initial N and mu_b ---------

nsp <- nrow(multi_sp@species_params)
nw  <- ncol(multi_sp@initial_n)

mu_mat <- n_mat <- matrix(NA, nrow = nsp, ncol = nw,
  dimnames = list(multi_sp@species_params$species, NULL))

for (sp in 1:length(spp_mods)){
  n_mat[sp, ] <- as.numeric(spp_mods[[sp]]@initial_n)
  mu_mat[sp, ]  <- as.numeric(spp_mods[[sp]]@mu_b)
}

n_mat['Hake', ] <- as.numeric(hake_model@initial_n)
mu_mat['Hake', ]  <- as.numeric(hake_model@mu_b)


multi_sp@initial_n[] <- n_mat
multi_sp@mu_b[] <- mu_mat

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)



## Predation ---------------------

sp <- multi_sp@species_params

celpars <- mizerEcopath::celtic_params@species_params

anc_sar_kerpars <- 
  readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/hake-mizer/data/anchovy_sardine_kernel_params.rds")

trcols <- c("pred_kernel_type", "kernel_exp", "kernel_l_l", "kernel_u_l", "kernel_l_r", "kernel_u_r")

for (col in trcols) {if (!col %in% names(sp)) { sp[[col]] <- NA}}

for (i in seq_len(nrow(sp))) {
  sp_name <- sp$species[i]
  if (sp_name %in% celpars$species) {
    row_cel <- celpars[celpars$species == sp_name, trcols]
    sp[i, trcols] <- row_cel
  }
}

sp[sp$species == "Four-spot megrim", trcols] <- sp[sp$species == "Megrim", trcols]

sp[which(sp$species=='Hake'), 'pred_kernel_type'] <- 'lognormal'     # comment for using celtic values

anc_sar_kerpars$pred_kernel_type <- anc_sar_kerpars$kernel_type
sp[c('Anchovy','Pilchard'),trcols] <- anc_sar_kerpars[,trcols]

# species_params(multi_sp) <- sp




## Interaction -------------

multi_sp <- alignResource(multi_sp, w_pp_cutoff = 1)

plotSpectra(multi_sp)

backup <- multi_sp
multi_sp <- backup

multi_sp <- setUniformInteraction( multi_sp)

plotDiet(multi_sp)
plotDeath(multi_sp)

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

interaction_matrix(multi_sp)
sp <- species_params(multi_sp); sp




sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)


