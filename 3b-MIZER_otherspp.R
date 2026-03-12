

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

# load( './output/alldata.RData')


# Multi-species model ------------------------

load('./input/Catch_spp.RData')
load('./input/other_spp_rmd.RData')


## Parameters --------------

tibble( other_spp$short_fb)
pars <- other_spp$short_fb

spps <- pars$red


## SSB --------------

ssb_long <- bind_rows(other_spp_ices, .id = "species") |> filter( Year >= 1992, Year <= 2023); ssb_long
ssb_spp <- ssb_long |> pivot_wider( names_from = species, values_from = SSB) |> arrange(Year); ssb_spp

load( './input/Hake_SS_Data.RData')
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


## Single spp fit for all spps -------------------

spp_mods <- catch_mods <- plots_mods <- list()

for(i in spps){
  
  ipars <- pars |> filter( red == i)
  issb <- bio |> filter( species == i) |> mutate( SSB = 10^6 * SSB, Bio = 10^6 * Bio)
  ilfd <- lfd[[i]] |> mutate(catch = value * 10^6)
  
  l_max <- 1.001 * max(ilfd$length + 5, na.rm = TRUE)
  
  isp <- newSingleSpeciesParams( 
    species_name = ipars$common, 
    w_mat = lwf( ipars$l_mat, ipars$a, ipars$b),
    w_max = max( lwf( ipars$l_max, ipars$a, ipars$b), lwf( l_max, ipars$a, ipars$b)), 
    n = 0.75, 
    # pred_kernel_type = 'lognormal', 
    beta = 11.33, 
    sigma = 0.46)
  
  isp@species_params$age_mat = laf( ipars$l_mat, ipars$l_inf, ipars$kvb, ipars$al0)
  isp@species_params$a = ipars$a
  isp@species_params$b = ipars$b 
 
  isp@species_params$biomass_observed <- issb$Bio
  isp@species_params$biomass_cutoff <- lwf(4,ipars$a,ipars$b)
 
  isp <- setBevertonHolt( isp, reproduction_level = 0.001)
   
  # ## Species params
  # # Max size (from landings)
  # l_max <- 1.001 * max(corLFD$length + 1, na.rm = TRUE)
  # sp$w_max = lwf(l_max, sp$a, sp$b)
  
  isp@species_params$d <- -0.1217
  
  isp <- isp@species_params
  
  imodel <- newAllometricParams(isp)
  
  
  igp <- data.frame(
    gear = "Total", 
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
  
  icatch <- ilfd |> mutate( dl = 1, species = ipars$common, gear = "Total")
  
  imodel <- matchCatch(imodel, catch = icatch)
  imodel <- metab_and_dens( imodel, imodel)
  
  spp_mods[[i]] <- imodel
  catch_mods[[i]] <- icatch
  plots_mods[[i]] <- plotYieldVsSize(imodel, x_var = "Length", catch = icatch)
  
}

plots_mods$Eng_encr
plots_mods$Lep_bosc
plots_mods$Lep_whif
plots_mods$Mic_pout
plots_mods$Sar_pilc
plots_mods$Sco_scom
plots_mods$Tra_trac


save( spp_mods, catch_mods, file = './output/other_spp.RData')


