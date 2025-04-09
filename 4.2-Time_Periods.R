
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~   Hake's MIZER model by time periods  ~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #



rm(list=ls())

library(dplyr)
library(ggplot2)
library(gridExtra)
library(plotly)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(mizerMR)

source( './scripts/aux_functions.R')


# SS Stock Assessment data --------------------------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results


# Biological parameters ----------------------------

load( './input/Bio_Pars.RData')


# Fishing Mortality --------------------------

load( './input/Catch.RData')   # './scripts/Natural_Mortality.R' with SS's catch information


# MIZER model -----------------

## SSB from SS ----------------
# average over last 30 years !!

ssbio_tp <- c()
for( i in 1:nrow(tps)) ssbio_tp[paste0(tps[i,1],' - ',tps[i,2])] <-
    sum(assessment$SSB[assessment$Year %in% tps[i,1]:tps[i,2]]*1e6)/5

tp_names <- c()
tp_mods <- tp_cld <- list()

for( i in 1:nrow(tps)){
  
  tp_names[i] <- paste0(tps[i,1],' - ',tps[i,2])
  tp_mods[[tp_names[i]]] <- bio_pars
  
  species_params(tp_mods[[tp_names[i]]])$biomass_observed <- ssbio_tp[i]
  species_params(tp_mods[[tp_names[i]]])$biomass_cutoff <- lwf(4,a,b)
  
  tp_mods[[tp_names[i]]] <- setBevertonHolt( tp_mods[[tp_names[i]]], reproduction_level = 0.001)
  
  tp_mods[[tp_names[i]]] <- tp_mods[[tp_names[i]]] |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 

  tp_cld[[tp_names[i]]] <- data.frame(species="Hake",gear=rep(gear_names,each=129),
                              length=rep(ld_tp[[tp_names[i]]]$length,9),dl=rep(1,9*129),
                              weight=rep(w(bio_pars),9),dw=rep(dw(bio_pars),9),
                              catch=c(ld_catch_tp[[tp_names[i]]]$catch))
  
}


#- Catch data by fishing gear ----------------

## TP1 

itp <- tp_names[1]

yield_gear <- Catch_tp[[itp]]


### Double sigmoid selectivity -------------

ds_sel <- data.frame(
  gear = gear_names, species = "Hake", catchability = 1, 
  sel_func = "double_sigmoid_length",
  l50 = c( 26.4, 27.5, 27.4, 14.1, 27.0, 28.2, 47.1, 53.4, 14.3),
  l50_right = c( 32.4, 33.2, 38.3, 20.6, 24.6, 34.3, 57.9, 54.4, 26.7),
  l25 = c( 24.1, 26.2, 25.5, 12.6, 26.3, 26.7, 42.5, 48.5, 12.1),
  l25_right = c( 38, 43.5, 44.1, 27, 30.4, 41.4, 66.2, 60.8, 28.4))

gear_params( tp_mods[[itp]]) <- ds_sel
gear_params( tp_mods[[itp]])$yield_observed <- Catch_tp[[itp]]$avg_yield

initial_effort( tp_mods[[itp]]) <- .25

tp_mods[[itp]] <- matchYield( tp_mods[[itp]])
tp_mods[[itp]] <- steady( tp_mods[[itp]])

gear_params( tp_mods[[itp]]) <- ds_sel
gear_params( tp_mods[[itp]])$yield_observed <- Catch_tp[[itp]]$avg_yield

# tp_mods[[itp]] <- tuneParams( tp_mods[[itp]], catch = tp_cld[[itp]])


#- Resulting model after parameter tuning and calibration --------

tp_mods[[itp]] <- readParams(paste0("./output/tp_mods_",1,".rds"))
tp_mods[[itp]] <- scaleDownBackground( tp_mods[[itp]], 1/3000000)

catchabilities <- gear_params(tp_mods[[itp]])$catchability




# 
# 
# 
# # Yield vs size plots ------------------------------------------
# 
# plist <- plist_mix <- list()
# 
# for( i in gear_names){ 
#   
#   plist[[i]] <- plotYieldVsSize( hake_model, species="Hake", gear=i,
#       catch=catch_lengths, x_var="Length", return_data=FALSE) + theme_bw() +
#     theme(legend.position = "none",axis.title.x=element_blank()) + ylim(c(0,0.13))
#   
#   plist_mix[[i]] <- plotYieldVsSize( mix_hake_model, species="Hake", gear=i,
#       catch=catch_lengths, x_var="Length", return_data=FALSE) + theme_bw() +
#     theme(legend.position = "none",axis.title.x=element_blank()) + ylim(c(0,0.13))
# 
# }
# 
# 
# grid.arrange( plist[[1]], plist[[2]], plist[[3]], plist[[4]], plist[[5]],
#     plist[[6]], plist[[7]], plist[[8]], plist[[9]], ncol = 3)
# 
# grid.arrange( plist_mix[[1]], plist_mix[[2]], plist_mix[[3]], plist_mix[[4]], plist_mix[[5]],
#               plist_mix[[6]], plist_mix[[7]], plist_mix[[8]], plist_mix[[9]], ncol = 3)
# 
# 
# pdf("./plots/selectivity.pdf", width = 10, height = 6, onefile = TRUE)
# 
# grid.arrange( plist[[1]], plist[[2]], plist[[3]], plist[[4]], plist[[5]],
#     plist[[6]], plist[[7]], plist[[8]], plist[[9]], ncol = 3)
# 
# grid.arrange( plist_mix[[1]], plist_mix[[2]], plist_mix[[3]], plist_mix[[4]], plist_mix[[5]],
#     plist_mix[[6]], plist_mix[[7]], plist_mix[[8]], plist_mix[[9]], ncol = 3)
# 
# dev.off()
# 
# 
# Save ----------------------

save.image( './output/hake_model.RData')

