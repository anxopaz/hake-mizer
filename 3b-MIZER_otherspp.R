
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
# source('./allometric/new_funs.R')



# Multi-species model ------------------------

load('./input/Catch_spp.RData')
load('./input/other_spp_rmd.RData')
load( './input/Hake_SS_Data.RData')


## Parameters --------------

pars <- other_spp$short_fb
tibble( pars)

spps <- pars$red


## SSB --------------

ssb_long <- bind_rows(other_spp_ices, .id = "species") |> filter( Year >= 1992, Year <= 2023); ssb_long
ssb_spp <- ssb_long |> pivot_wider( names_from = species, values_from = SSB) |> arrange(Year); ssb_spp


## LFD ------------------

LFD_data <- lapply(LFD_spp_tp, function(df) { split(df, df$spp)})
LFD_data$aver_y |> bind_rows() |> pivot_wider( names_from = spp, values_from = value) |> arrange(length)
LFD_total <- lapply(LFD_spp_tp, function(df) { split(df, df$spp) |> 
    lapply(function(x) sum(x$value, na.rm = TRUE))})


# MIZER fit (first time period) ---------------------

lfd <- LFD_data$aver_y; lfd
lfd_total <- LFD_total$aver_y; lfd_total

bio <- ssb_spp |> filter(Year %in% years$aver_y) |> 
  summarise(across(-Year, ~ mean(.x, na.rm = TRUE))) |> 
  pivot_longer(cols = everything(), names_to = "species", values_to = "SSB") 


## Using hake as a reference -----------

load( './output/hake_model.RData')


# # For Biomass/SSB ratio
# 
# ratio_sb <- mean( assessment$biomass/assessment$SSB)
# bio <- bio |> mutate( Bio = SSB * ratio_sb); bio


# For mortality

herg <- getEReproAndGrowth(hake_model)
hwmat <- hake_model@species_params$w_mat

wmatidx <- which( hake_model@w >= hwmat)[1]

hergmat <- approx( hake_model@w, y= herg, xout=hwmat)$y

plot( hake_model@w[1:(wmatidx+5)], herg[1:(wmatidx+5)], 
      type='l', xlab = 'Weight (g)', ylab = 'EReproAndGrowth')

abline( v = hwmat, col='red', lty = 'dashed')
abline( h = hergmat, col='red', lty = 'dashed')

hmu <- hake_model@species_params$mu_mat

xi <- (hergmat/hmu)/hwmat; xi



## Loop ---------------

spp_mods <- catch_mods <- plots_mods <- list()

for(i in spps){
  
  ipars <- pars |> filter( red == i)
  # issb <- bio |> filter( species == i) |> mutate( SSB = 10^6 * SSB, Bio = 10^6 * Bio)
  issb <- bio |> filter( species == i) |> mutate( SSB = 10^6 * SSB)
  ilfd <- lfd[[i]] |> mutate(catch = value * 10^6)
  icom <- ipars$common
  
  # ftype <- ifelse( icom %in% c( 'Mackerel', 'Four-spot megrim', 'Megrim'), 1, 2)
  ftype <- ifelse( icom %in% c( 'Mackerel', 'Four-spot megrim', 'Megrim'), 2, 2)
  
  l_max <- 1.001 * max(ilfd$length + 5, na.rm = TRUE)
  
  isp <- data.frame( 
    species = icom, 
    w_mat = lwf( ipars$l_mat, ipars$a, ipars$b),
    w_max = max( lwf( ipars$l_max, ipars$a, ipars$b), lwf( l_max, ipars$a, ipars$b)))
  
  isp$age_mat = laf( ipars$l_mat, ipars$l_inf, ipars$kvb, ipars$al0)
  isp$a = ipars$a
  isp$b = ipars$b 
  
  isp$biomass_observed <- issb$SSB   # issb$Bio for Biomass
  isp$biomass_cutoff <- isp$w_mat    # lwf(4,ipars$a,ipars$b) if Biomass
  
  imodel <- newAllometricParams( isp, max_w = hake_model@species_params$w_max)
  
  if( ftype == 2){
    
    igp <- data.frame(
      gear = "Demersales", 
      species = icom, 
      catchability = 1,
      sel_func = "double_sigmoid_length",
      l50 = ipars$l_mat,
      l25 = ipars$l_mat*0.8,
      l50_right = ipars$l_mat*1.05,
      l25_right = ipars$l_mat*1.2,
      yield_observed = sum(ilfd$catch)
    )
    
  } else if( ftype == 1){
    
    igp <- data.frame(
      gear = "Demersales", 
      species = icom, 
      catchability = 1,
      sel_func = "sigmoid_length",
      l50 = ipars$l_mat,
      l25 = ipars$l_mat*0.8,
      yield_observed = sum(ilfd$catch)
    )
  }

  gear_params(imodel) <- igp
  initial_effort(imodel) <- 1
  
  yield <- as.numeric(getYield(imodel))
  igp$catchability <- igp$yield_observed / yield
  gear_params(imodel) <- igp
  
  icatch <- ilfd |> mutate( dl = 1, species = icom, gear = "Demersales") |>
    select( length, catch, dl, species, gear)
  
  mu_mat_lim <- 4
  erepro <- 1
  
  while( erepro > 0.5){
    
    mu_mat_lim <- mu_mat_lim - 0.05
    imodel <- suppressWarnings( matchCatch( imodel, catch = icatch, mu_mat_lim = mu_mat_lim, map = NULL))
    erepro <- imodel@species_params$erepro
  
  }
  
  spp_mods[[i]] <- imodel
  catch_mods[[i]] <- icatch
  plots_mods[[i]] <- plotYieldVsSize(imodel, x_var = "Length", catch = icatch)
  
  print( data.frame(spp_mods[[i]]@gear_params)[,(-c(1,2))])
  
  cat('\n')
  
  print( data.frame( m = spp_mods[[i]]@species_params$m,
                     erepro = spp_mods[[i]]@species_params$erepro, 
                     mu_mat = spp_mods[[i]]@species_params$mu_mat))
  
  print(paste0('erepro for a reproduction_level=0.5: ',round( setBevertonHolt(imodel, reproduction_level = 0.5)@species_params$erepro,4)))
  cat('\n'); cat('\n')
  
  print(plotYieldVsSize( spp_mods[[i]], x_var = "Length", catch = catch_mods[[i]]))
  
}


for(i in spps) print(data.frame(SSB_obs = spp_mods[[i]]@species_params$biomass_observed/10^6, 
                                SSB_estwbio = getSSB(spp_mods[[i]])/10^6,
                                SSB_est = getBiomass(spp_mods[[i]], use_cutoff = T)/10^6))

print( data.frame( species = 'Hake', m = hake_model@species_params$m,
                   erepro = hake_model@species_params$erepro, 
                   mu_mat = hake_model@species_params$mu_mat))

plotYieldVsSizeByGear( hake_model, catch = hake_catch)


save( spp_mods, catch_mods, plots_mods, file = './output/other_spp.RData')


