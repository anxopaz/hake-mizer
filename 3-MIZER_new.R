
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~   Hake's MIZER model  ~~~~~~~~~~~~~~ #
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

source( './scripts/aux_functions.R')


# Biological parameters ----------------------------

load( './input/Bio_Pars.RData')     # './1-Bio_Pars.R' with Biological Parameters


# Fishing Mortality --------------------------

load( './input/Catch.RData')   # './2-Catch.R' with Catch and LFD data


# SSB / Bio ----------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

quantity <- 'biomass' # 'SSB'

ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)

obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm

species_params(bio_pars)$biomass_observed <- obs_q
species_params(bio_pars)$biomass_cutoff <- b_min

bio_pars <- setBevertonHolt( bio_pars,  # Rdd = Rdi * (Rmax/(Rdi+Rmax))
   reproduction_level = 0.001)          # rep_level = Rdd/Rmax (density dependance degree)



# MIZER model --------------------------------

library(mizerEcopath)


## Species params --------------------------

sp <- bio_pars@species_params |>
  select(species, w_mat, age_mat, w_max, a, b, n, pred_kernel_type, beta, sigma)

# Biomass
quantity <- 'biomass' # 'SSB'
ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)
obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm
sp$biomass_observed <- obs_q
sp$biomass_cutoff <- b_min

# Max size (from landings)
l_max <- 1.001 * max(corLFD$length + 1, na.rm = TRUE)
sp$w_max = lwf(l_max, sp$a, sp$b)

# Mortality exponent (from SS)
sp$d <- -0.1217


## Allometric params --------------------

# We now set up a Hake model with power law encounter rate and
# power law mortality rate.
hake_model <- newAllometricParams(sp)
# This is not yet calibrated


## Fishing (1 gear) -------------------------

# We work with only one "Total" gear
gp <- data.frame(
  gear = "Total", species = "Hake", catchability = 1,
  sel_func = "sigmoid_length",
  l50 = 30,
  l25 = 28,
  yield_observed = sum(corLFDs$catch)
)
gear_params(hake_model) <- gp
initial_effort(hake_model) <- 1

# Set initial catchability to get the observed yield
# This will be calibrated properly in `matchCatch()` below.
yield <- getYield(hake_model)
gp$catchability <- gp$yield_observed / yield
gear_params(hake_model) <- gp

# Transform landings data to required shape, adding landings from
# all gears to give total landings
catch <- corLFD |> group_by(length) |> summarise(catch = sum(number))
catch$dl = 1
catch$species = "Hake"
catch$gear = "Total"


### MatchCatch -------------

# Calibrate to total yield and landings size distribution
hake_model <- matchCatch(hake_model, catch = catch)

plotYieldVsSize(hake_model, x_var = "Length", catch = catch)


## Fishing (different gear) -------------------------
# Same for different gears
catch2 <- corLFD |>
  mutate( catch = number, gear = fleet, dl = 1, species = 'Hake') |>
  select( length, catch, dl, species, gear)

gear_names <- unique( catch2$gear)

gp <- data.frame(
  gear = gear_names, 
  species = "Hake",
  sel_func = ifelse( gear_names %in% c('palangre','vol'), 'sigmoid_length', 'double_sigmoid_length'),
  l50 =       c( 28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
  l25 =       c( 23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
  l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, NA,   58.0, 54.4, NA  ),
  l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, NA,   67.9, 60.8, NA  ),
  yield_observed = corLFDs$catch,
  catchability = corLFDs$catch/sum(corLFDs$catch))

gear_params(hake_model) <- gp

initial_effort(hake_model) <- 1

### MatchCatch -------------

source('./allometric/new_funs.R')
hake_model <- matchCatch(hake_model, catch = catch2)
plotYieldVsSizeByGear(hake_model, catch = catch2)


## Set metabolic loss rate --------------

# So far we ran the model without metabolic loss. If we now introduce this loss, we have to increase the encounter rate to make up for this.
species_params(hake_model)$ks <- bio_pars@species_params$ks
ext_encounter(hake_model) <- ext_encounter(hake_model) +
  metab(hake_model) / species_params(hake_model)$alpha


## Add density dependencies ------------------------

# We now have a model whose steady state matches the landings data, but we still need to calibrate its sensitivity to changes away from the steady state, for example its sensitivity to changes in fishing. We don't have a good way to choose the following three parameters yet,  so we'll just make up values for now.

# Set reproduction level
hake_model <- setBevertonHolt(hake_model, reproduction_level = 0.6)

# Set feeding level
hake_model <- setFeedingLevel(hake_model, 0.6)

# Set resource level
hake_model <- alignResource(hake_model) |>
  setResourceInteraction(resource_dynamics = "resource_semichemostat")
resource_level(hake_model) <- 1/2


# Cannibalism ----------------------

# I'll assume below that 14% of the total diet comes from cannibalism. You will get a warning if you try to increase this and we should discuss this.
diet_matrix <- matrix(c(0.14, 0.86), ncol = 2,
                      dimnames = list(predator = "Hake",
                                      prey = c("Hake", "other")))

cannibal_hake_model <- matchDiet(hake_model, diet_matrix)
plotDietX(cannibal_hake_model)


# Steady state ------------------------

sim <- project(hake_model, t_max = 8)
plotBiomass(sim)
sim <- project(cannibal_hake_model, t_max = 8)
plotBiomass(sim)

# Run simulations with increased fishing effort ---------------
sim12 <- project(hake_model, effort = 1.2, t_max = 12, t_save = 0.2)
sim_cannibal12 <- project(cannibal_hake_model, effort = 1.2, t_max = 12, t_save = 0.2)


# Compare biomasses --------------------------

bio12 <- melt(getBiomass(sim12))
bio12$Cannibalism <- "Off"
bio_cannibal12 <- melt(getBiomass(sim_cannibal12))
bio_cannibal12$Cannibalism <- "On"
bio <- rbind(bio12, bio_cannibal12)
ggplot(bio) + geom_line(aes(x = time, y = value, colour = Cannibalism)) 



# MSY curve --------------------

yieldf <- plotYieldVsF(hake_model, species = 'Hake', F_range = seq(0, 1.2, 0.02)); yieldf
max(yieldf$data$yield); yieldf$data$F[which.max(yieldf$data$yield)]

yieldf2 <- plotYieldVsF(cannibal_hake_model, species = 'Hake', F_range = seq(0, 1.2, 0.02)); yieldf2
max(yieldf2$data$yield); yieldf2$data$F[which.max(yieldf2$data$yield)]


# Growth vs NatMort (k/M) -------------------------

## Growth ---------------

mizergrowth <- plotGrowthCurves(hake_model)$data
gage <- mizergrowth$Age

ssgrowthf <- lwf( alf(gage, Linf_f, Kvb, al0_f), a, b)
ssgrowthm <- lwf( alf(gage, Linf_m, Kvb, al0_m), a, b)
ssgrowth <- lwf( alf(gage, Linf, Kvb, al0), a, b)

grdf <- data.frame( Age = gage, Mizer = mizergrowth$value, SS_f = ssgrowthf, ss_m = ssgrowthm, ss_mean = ssgrowth)

grdf %>% pivot_longer( cols = Mizer:ss_mean) %>% 
  ggplot( aes(x = Age, y = value, color = name)) + geom_line() + theme_bw()


## NatMort ------------------------

mizerextmort <- getExtMort(hake_model)
sizes <- as.numeric(colnames(mizerextmort))

testnm <- hake_model@species_params$mu_mat*sizes^hake_model@species_params$d
testnm2 <- hake_model@species_params$mu_mat*(sizes/w50)^hake_model@species_params$d

ssfitextmort <- mu0*sizes^(d)

M <- replist$Natural_Mortality_Bmark
M_female <- as.numeric(subset( M, Seas==1 & Settlement==1 & Sex==1)[-(1:4)])
M_male <- as.numeric(subset( M, Seas==1 & Settlement==1 & Sex==2)[-(1:4)])
NatM <- ( M_female + M_male)/2    # Mean of both sexes

weights <- lwf( c(4, alf( 0:15, Linf, Kvb, al0)[-1]),a,b)

mdf <- rbind(
  data.frame( Weight = weights, Sex = 'SS Female', Mortality = M_female),
  data.frame( Weight = weights, Sex = 'SS Male', Mortality = M_male),
  data.frame( Weight = weights, Sex = 'SS Average', Mortality = NatM),
  data.frame( Weight = sizes, Sex = 'SS refit', Mortality = ssfitextmort),
  data.frame( Weight = sizes, Sex = 'MIZER', Mortality = as.numeric(mizerextmort)),
  data.frame( Weight = sizes, Sex = 'MIZER test', Mortality = testnm),
  data.frame( Weight = sizes, Sex = 'MIZER test 2', Mortality = testnm2))

ggplot( mdf, aes( x = Weight, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point() + xlim(NA, max(weights)) + ylim(NA, 2) + 
  theme(legend.title=element_blank())


# Save ----------------------

hake_mizer <- hake_model

save.image( './output/hake_model.RData')

modelo <- hake_mizer
save( modelo, file = "./fit.RData")




# Time periods --------------------------------

# Loop for time periods

tpdir <- paste0( getwd(), '/output/time_periods/')
dir.create( path = tpdir, showWarnings = TRUE, recursive = TRUE)

tp_mods <- list()

years <- list( years[[2]], years[[3]], years[[4]])
for(i in 1:length(years)) names(years)[[i]] <- paste0( years[[i]][1],' - ', years[[i]][5]) 

# for( i in 1:length(years)){
#   
#   ny <- names(years)[i]
#   vy <- years[[i]]
#   ychar <- paste0( vy[1],'-',vy[length(vy)])
#   
#   ssbio_tp <- sum(assessment[,quantity][assessment$Year %in% vy]*1e6)/length(vy)
#   
#   itp_mod <- hake_mizer
#   
#   species_params(itp_mod)$biomass_observed <- ssbio_tp
#   
#   itp_mod <- itp_mod |>
#     calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
#     calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()
#   
#   gear_params( itp_mod)$yield_observed <- corLFD_sum[[i]]$catch
#   
#   tp_mods[[ychar]] <- MIZER( model =  itp_mod, catch = corLFD_list[[i]], compiler = compiler,
#                              plot = T, plot_dir = paste0(tpdir,ychar,'/'))
#   
#   plotSpectra( tp_mods[[ychar]]) + theme_bw()
#   ggsave(paste0(tpdir,ychar,'/spectra.jpg'), width = 9, height = 7) 
#   
#   plotSpectra( tp_mods[[ychar]], power = 2) + theme_bw()
#   ggsave(paste0(tpdir,ychar,'/spectra2.jpg'), width = 9, height = 7) 
#   
# }
# 
# 
# tpdf <- NULL
# for(i in names(tp_mods)){ tpdf <- rbind( tpdf, spf( tp_mods[[i]], name = i))}
# 
# ggplot( tpdf %>% filter(w>10), aes( x = w, y = value, color = Model)) + 
#   geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
#   labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Period')
# 
# ggsave( paste0( tpdir, 'spectra_comparison.jpg'), width = 12, height = 7)
# ggsave( paste0( tpdir, 'spectra_comparison_ppt.jpg'), width = 6, height = 4)
# 
# ggplot( tpdf %>% filter(w>10), aes( x = w, y = value2, color = Model)) + 
#   geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
#   labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Period')
# 
# ggsave( paste0( tpdir, 'spectra_comparison2.jpg'), width = 12, height = 7)
# ggsave( paste0( tpdir, 'spectra_comparison2_ppt.jpg'), width = 6, height = 4)
# 
# 
# tp_table <- matrix( NA, nrow = length(tp_mods), ncol = length(all_pars), 
#                     dimnames = list( names(tp_mods), all_pars))
# 
# for(i in names(tp_mods)){ tp_table[i,] <- parsf( tp_mods[[i]])}
# tp_table





# Multi-species model ------------------------

## Single spp fit for all spps -------------------

## Joint spps -------------------------------


# Save all ----------------------

save.image( './output/alldata.RData')




