
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

for(i in 1:length(spp_mods)) print(getBiomass(spp_mods[[i]], use_cutoff = TRUE)/spp_mods[[i]]@species_params$biomass_observed)
getBiomass(hake_model, use_cutoff = TRUE)/hake_model@species_params$biomass_observed


# Bind species in MS model -----------

allcatch <- c( list( Hake = hake_catch), catch_mods)
allmods <- c( Hake = hake_model, spp_mods)

catchdf <- bind_rows( allcatch) 

save( allmods, allcatch, catchdf, file='./input/mods.RData')

msm <- mizerEcopath::bindParams( allmods)

msm@linecolour[1:8] <- NS_params@linecolour[1:8]

sim <- project(msm, t_max = 10)
plotBiomass(sim)


## Predation ---------------------

# Taking predation kernel parameters from the Celtic Sea model for all species
# except Anchovy and Sardine that are determined from Mariella's data
# and Hake from own data

sp <- msm@species_params

celpars <- mizerEcopath::celtic_params@species_params

anc_sar_kerpars <- 
  readRDS("data/anchovy_sardine_kernel_params.rds")

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

msm@species_params <- sp


## Feeding Levels --------------

msm <- setFeedingLevels(msm, f = 0.6, f_c = 0.2)

sim <- project(msm, t_max = 10)
plotBiomass(sim)


## Interaction -------------

msm <- setUniformInteraction( msm)

plotDiet(msm)
plotDeath(msm)
getDietMatrix(msm)
plotSpectra(msm)

hake_diet <- getDietMatrix(msm, min_w_pred = 100)['Hake',]; hake_diet
round( hake_diet/ sum(hake_diet), 5)
round( hake_diet[-(9:10)]/ sum(hake_diet[-(9:10)]), 5)

sim <- project(msm, t_max = 10)
plotBiomass(sim)



## Resource ------------

backup <- msm

msm <- backup

msm <- alignResource(msm, w_pp_cutoff = 1)
# msm <- setResourceInteraction(msm, "resource_semichemostat")

resource_level(msm)
plotSpectra(msm)

sim <- project(msm, t_max = 10)
plotBiomass(sim)


sim <- setBevertonHolt(msm, reproduction_level = 0.5)

sim <- project(msm, t_max = 10)
plotBiomass(sim)

nsp <- nrow(species_params(msm))
init_eff <- msm@initial_effort

sim08 <- project( msm, effort = init_eff * 0.8, t_max = 10)
plotBiomass(sim08)

sim12 <- project( msm, effort = init_eff * 1.2, t_max = 10)
plotBiomass(sim12)
