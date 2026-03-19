
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

library(mizerEcopath)

source( './scripts/aux_functions.R')


# Species params -----------------

## Hake --------------------------------

load( './input/Bio_Pars.RData')     # './1a-BioPars_hake.R' with Biological Parameters
load( './input/Catch.RData')   # './2-Catch.R' with Catch and LFD data
load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

sp <- data.frame( bio_pars@species_params |>
  select(species, w_mat, age_mat, w_max, a, b, pred_kernel_type, beta, sigma))

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
hake_d <- -0.1217


## Other spp ---------------

load('./input/other_spp_rmd.RData')
load('./input/Catch_spp.RData')

LFD_data <- lapply(LFD_spp_tp, function(df) { split(df, df$spp)})
LFD_total <- lapply(LFD_spp_tp, function(df) { split(df, df$spp) |> lapply(function(x) sum(x$value, na.rm = TRUE))})
lfd <- LFD_data$aver_y; lfd
lfd_total <- LFD_total$aver_y; lfd_total
for(i in 1:length(lfd)) lfd[[i]] <- lfd[[i]] |> mutate(catch = value * 10^6)
df_max <- purrr::map_df(names(lfd), ~ { tibble( red = .x, max_length = max(lfd[[.x]]$length, na.rm = TRUE))})

ssb_long <- bind_rows(other_spp_ices, .id = "species") |> filter( Year >= 1992, Year <= 2023); ssb_long
ssb_spp <- ssb_long |> pivot_wider( names_from = species, values_from = SSB) |> arrange(Year); ssb_spp

ratio_sb <- mean( assessment$biomass/assessment$SSB)

bio <- ssb_spp |> filter(Year %in% years$aver_y) |> 
  summarise(across(-Year, ~ mean(.x, na.rm = TRUE))) |> 
  pivot_longer(cols = everything(), names_to = "species", values_to = "SSB") |>
  mutate( Bio = SSB * ratio_sb, red = species)  


pars <- other_spp$short_fb
tibble( pars)


sps <- data.frame(
  species = pars$common, 
  w_mat = lwf(pars$l_mat, pars$a, pars$b),
  age_mat = laf( pars$l_mat, pars$l_inf, pars$kvb, pars$al0),
  w_max = lwf( pars$l_max, pars$a, pars$b), 
  a = pars$a, b = pars$b, 
  pred_kernel_type = 'lognormal', beta = sp$beta, sigma = sp$sigma
)

for(i in 1:nrow(sps)) {
  sps[i,'w_max'] <- max(sps[i,'w_max'], 
    lwf( 1.001 * as.numeric(df_max[which(df_max$red==pars[i,'red']),'max_length']) + 5, 
         sps[i,'a'], sps[i,'b']))}

sps$biomass_observed <- bio$Bio[ match(pars$red, bio$red) ] * 1e6
sps$biomass_cutoff <- lwf( 4, pars$a, pars$b)

sps

sp <- rbind(sp, sps)

sp$d <- -0.25; sp[1,'d'] <- hake_d
sp$q <- sp$n <- 0.75

msm <- newAllometricParams(sp, no_w = 400)

# msm <- setBevertonHolt( msm, reproduction_level = 0.001)

sp <- msm@species_params

# celpars <- mizerEcopath::celtic_params@species_params
# 
# trcols <- c("pred_kernel_type", "kernel_exp", "kernel_l_l", "kernel_u_l", "kernel_l_r", "kernel_u_r")
# 
# for (col in trcols) {if (!col %in% names(sp)) { sp[[col]] <- NA}}
# 
# for (i in seq_len(nrow(sp))) {
#   sp_name <- sp$species[i]
#   if (sp_name %in% celpars$species) {
#     row_cel <- celpars[celpars$species == sp_name, trcols]
#     sp[i, trcols] <- row_cel
#   }
# }
# 
# sp[sp$species == "Four-spot megrim", trcols] <- sp[sp$species == "Megrim", trcols]
# 
# sp[which(sp$species=='Hake'), 'pred_kernel_type'] <- 'lognormal' 
# 
# species_params(msm) <- sp




# Gear params -------------------------

# hake_c <- corLFD |>
#   mutate( catch = number, gear = fleet, dl = 1, species = 'Hake') |>
#   select( length, catch, dl, species, gear)
# 
# gear_names <- unique( hake_c$gear)
# 
# gp <- data.frame(
#   gear = gear_names, 
#   species = "Hake",
#   sel_func = ifelse( gear_names %in% c('palangre','vol'), 'sigmoid_length', 'double_sigmoid_length'),
#   l50 =       c( 28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
#   l25 =       c( 23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
#   l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, NA,   58.0, 54.4, NA  ),
#   l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, NA,   67.9, 60.8, NA  ),
#   yield_observed = corLFDs$catch,
#   catchability = corLFDs$catch/sum(corLFDs$catch))

gp <- data.frame(
  gear = "demersales", species = "Hake",
  sel_func = "sigmoid_length",
  l50 = 30,
  l25 = 28,
  yield_observed = sum(corLFDs$catch), catchability = 1
)


spns <- sp$species[-1]

totalcatch <- purrr::map_df(names(lfd_total), ~ {
  tibble( species = .x, total_catch = lfd_total[[.x]])})

# gp <- rbind(gp, data.frame( gear = "demersales", species = spns, 
#   sel_func = "sigmoid_length", l50 = pars$l_mat, l25 = pars$l_mat*0.8, 
#   l50_right = NA, l25_right = NA,
#   yield_observed = totalcatch$total_catch*10^6, catchability = 1))

gp <- rbind(gp, data.frame( gear = "demersales", species = spns, 
                            sel_func = "sigmoid_length", l50 = pars$l_mat, l25 = pars$l_mat*0.8, 
                            yield_observed = totalcatch$total_catch*10^6, catchability = 1))

gear_params(msm) <- gp

names(lfd)  <- pars$common
for(i in 1:length(lfd)) lfd[[i]]$gear <- 'demersales'



hake_c <- corLFD |> group_by(length) |> summarise(catch = sum(number)) |>
  mutate( gear = "demersales")

lfds <- c( list(Hake = hake_c), lfd)

for(i in 1:length(lfds)){
  lfds[[i]]$dl <- 1
  lfds[[i]]$species <- names(lfds)[[i]]
  lfds[[i]] <- lfds[[i]] |> select(species, gear, length, catch, dl)}

catch <- data.frame( bind_rows(lfds))


msm@species_params
msm@gear_params
catch

initial_effort(msm) <- 1

# source('./allometric/new_funs.R')

fmsm <- matchCatch(msm, catch = catch)

plotYieldVsSizeByGear(fmsm, catch)

plotYieldVsSize(fmsm, x_var = "Length", catch = catch, species = 'Megrim')

load('./output/other_spp.RData')
plots_mods$Lep_whif


load( 'prefit.RData')


tibble(multi_sp@gear_params |> filter(species != 'Hake'))
tibble(fmsm@gear_params |> filter(species != 'Hake'))

tibble(multi_sp@gear_params |> filter(species == 'Hake'))
tibble(fmsm@gear_params |> filter(species == 'Hake'))

tibble(multi_sp@species_params)
tibble(fmsm@species_params)

sp1 <- tibble(multi_sp@species_params)
sp2 <- tibble(fmsm@species_params)
sp2_ordered <- sp2[rownames(sp1), colnames(sp1)]

sp1
sp2_ordered


sim <- project( fmsm, t_max=10)
plotBiomass( sim)


params <- msm
species = 'Anchovy'
lambda = 2.05
yield_lambda = 1
production_lambda = 1


hake_model_fleets <- matchCatch(hake_model, catch = hake_c)
plotYieldVsSizeByGear(hake_model_fleets, catch = hake_c)
