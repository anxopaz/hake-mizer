
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

for(i in 1:length(spp_mods)) print(getBiomass(spp_mods[[i]]), use_cutoff = TRUE)
getBiomass(hake_model, use_cutoff = TRUE)


# Other option -----------

allcatch <- c( list( Hake = hake_catch), catch_mods)
allmods <- c( Hake = hake_model, spp_mods)

catchdf <- bind_rows( allcatch) 

save( allmods, file='./mods.RData')

multi_sp <- mizerEcopath::bindParams( allmods)


## Gear parameters --------------
gp <- multi_sp@gear_params

backup <- multi_sp
multi_sp <- backup

multi_sp <- setBevertonHolt(multi_sp, reproduction_level = 0)

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

## Predation ---------------------

# Taking predation kernel parameters from the Celtic Sea model for all species
# except Anchovy and Sardine that are determined from Mariella's data
# and Hake from own data

sp <- multi_sp@species_params

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

multi_sp@species_params <- sp

## Feeding Levels
# The next line is needed only due to bug in newAllometricParams()
multi_sp@species_params$p <- multi_sp@species_params$n

multi_sp <- setFeedingLevels(multi_sp, f = 0.6, f_c = 0.2)


sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

## Interaction -------------

backup <- multi_sp
multi_sp <- backup

multi_sp <- setUniformInteraction( multi_sp)
multi_sp@linecolour[1:8] <- NS_params@linecolour[1:8]
plotDiet(multi_sp)
plotDeath(multi_sp)
getDietMatrix(multi_sp)
plotSpectra(multi_sp)

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

## Resource

multi_sp <- alignResource(multi_sp, w_pp_cutoff = 1)
multi_sp <- setResourceInteraction(multi_sp, "resource_semichemostat")

resource_level(multi_sp)
plotSpectra(multi_sp)

backup2 <- multi_sp
multi_sp <- backup2

sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)

interaction_matrix(multi_sp)
sp <- species_params(multi_sp); sp




sim <- project(multi_sp, t_max = 10)
plotBiomass(sim)


