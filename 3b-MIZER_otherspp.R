
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~   Other spp MIZER model  ~~~~~~~~~~~~~~ #
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
library(mizerEcopath)

source( './scripts/aux_functions.R')
source('./allometric/new_funs.R')



# Multi-species model ------------------------

load('./input/Catch_spp.RData')
load('./input/other_spp_rmd.RData')
load( './output/hake_models.RData')
load( './input/Hake_SS_Data.RData')


## Parameters --------------

pars <- other_spp$short_fb
tibble( pars)

spps <- pars$red


## SSB --------------

ssb_long <- bind_rows(other_spp_ices, .id = "species") |> filter( Year >= 1992, Year <= 2023); ssb_long
ssb_spp <- ssb_long |> pivot_wider( names_from = species, values_from = SSB) |> arrange(Year); ssb_spp


# assessment <- assessment |> filter( Year %in% 2014:2023)
ratio_sb <- mean( assessment$biomass/assessment$SSB)

## LFD ------------------

LFD_data <- lapply(LFD_spp_tp, function(df) { split(df, df$spp)})
LFD_total <- lapply(LFD_spp_tp, function(df) { split(df, df$spp) |> lapply(function(x) sum(x$value, na.rm = TRUE))})


# MIZER fit (first time period) ---------------------

lfd <- LFD_data$aver_y; lfd
lfd_total <- LFD_total$aver_y; lfd_total

bio <- ssb_spp |> filter(Year %in% years$aver_y) |> 
  summarise(across(-Year, ~ mean(.x, na.rm = TRUE))) |> 
  pivot_longer(cols = everything(), names_to = "species", values_to = "SSB") |>
  mutate( Bio = SSB * ratio_sb)

bio


# Single fit ---------------

spp_mods <- catch_mods <- plots_mods <- list()

for(i in spps){
  
  ipars <- pars |> filter( red == i)
  issb <- bio |> filter( species == i) |> mutate( SSB = 10^6 * SSB, Bio = 10^6 * Bio)
  ilfd <- lfd[[i]] |> mutate(catch = value * 10^6)
  
  l_max <- 1.001 * max(ilfd$length + 5, na.rm = TRUE)
  
  isp <- data.frame( 
    species = ipars$common, 
    w_mat = lwf( ipars$l_mat, ipars$a, ipars$b),
    w_max = max( lwf( ipars$l_max, ipars$a, ipars$b), lwf( l_max, ipars$a, ipars$b)))
  
  isp$age_mat = laf( ipars$l_mat, ipars$l_inf, ipars$kvb, ipars$al0)
  isp$a = ipars$a
  isp$b = ipars$b 
  
  isp$biomass_observed <- issb$Bio
  isp$biomass_cutoff <- lwf(4,ipars$a,ipars$b)
  
  imodel <- newAllometricPars( isp, max_w = hake_model@species_params$w_max)
  
  igp <- data.frame(
    gear = "Demersales", 
    species = ipars$common, 
    catchability = 1,
    sel_func = "sigmoid_length",
    l50 = ipars$l_mat,
    l25 = ipars$l_mat*0.8,
    yield_observed = sum(ilfd$catch)
  )
  
  gear_params(imodel) <- igp
  initial_effort(imodel) <- 1
  
  yield <- getYield(imodel)
  igp$catchability <- igp$yield_observed / yield
  gear_params(imodel) <- igp
  
  icatch <- ilfd |> mutate( dl = 1, species = ipars$common, gear = "Demersales") |>
    select( length, catch, dl, species, gear)
  
  imodel <- matchCatch( imodel, catch = icatch)
  
  imodel <- steadySingleSpecies( imodel) 
  imodel <- setBevertonHolt( imodel, reproduction_level = 0.9)
  
  print(plotYieldVsSize(imodel, x_var = "Length", catch = icatch))
  print(plotBiomass(project(imodel,t_max=10)))
  
  spp_mods[[i]] <- imodel
  catch_mods[[i]] <- icatch
  plots_mods[[i]] <- plotYieldVsSize(imodel, x_var = "Length", catch = icatch)
  
}



save( spp_mods, catch_mods, plots_mods, file = './output/other_spp.RData')


