

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
library(TMB)

library(mizerEcopath)

source( './scripts/aux_functions.R')

source('./allometric/new_funs.R')

load( './output/alldata.RData')
load( './output/other_spp.RData')


multisp <- hake_model

params <- addSpecies( multisp, spp_mods[[1]])





