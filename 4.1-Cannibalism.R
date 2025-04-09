
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~   Hake's MIZER model with cannibalism  ~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rm(list=ls())

library(mizer)
library(mizerExperimental)
load( './output/hake_model.RData')
load( './data/Diet.RData')

plotSpectra( hake_model, power = 2)

cannibal <- hake_model2 <- hake_model

ext_mort( hake_model2)[1:60] <- 0

# Turn on cannibalism ----------------------------
interaction_matrix( cannibal)[] <- as.numeric(by_prey[which(by_prey$Prey=='M.merluccius'),2])
interaction_matrix( hake_model2)[] <- as.numeric(by_prey[which(by_prey$Prey=='M.merluccius'),2])

getPredMort( cannibal) - getPredMort( hake_model2)

# Reduce external mortality because Hake predation is now modeled explicitly
cannibal2 <- cannibal
ext_mort( cannibal) <- ext_mort( cannibal) - getPredMort( cannibal)

# Check that steady state has not changed
p_test <- steadySingleSpecies( cannibal)
all.equal( initialN(cannibal), initialN(p_test), tolerance = 1e-5)

plotSpectra( cannibal, power = 2, total = TRUE)
plotSpectra( cannibal2, power = 2, total = TRUE)

cannibal <- setBevertonHolt( cannibal, reproduction_level = 0.01)
cannibal2 <- setBevertonHolt( cannibal2, reproduction_level = 0.01)

plotYieldVsF( cannibal, species="Hake",F_max=15)
plotYieldVsF( cannibal2, species="Hake",F_max=15)


## Fit ------------

# tuneParams( cannibal, catch = catch_lengths)
# tuneParams( cannibal2, catch = catch_lengths)

cannibal <- readParams("./output/cannibal.rds")
cannibal2 <- readParams("./output/cannibal2.rds")

# Comparison -----------------

plotDiet( hake_model)
plotDiet( cannibal)
plotDiet( cannibal2)

plotDeath( hake_model)
plotDeath( cannibal)
plotDeath( cannibal2)

plotSpectra2(hake_model, name1 = "Base Model",
             cannibal2, name2 = "With cannibalism",
             power = 2, resource = FALSE, wlim = c(10, NA))

plotSpectra2(hake_model, name1 = "Base Model",
             cannibal, name2 = "With cannibalism - NM",
             power = 2, resource = FALSE, wlim = c(10, NA))

plotSpectra2(cannibal2, name1 = "With cannibalism",
             cannibal, name2 = "With cannibalism - NM",
             power = 2, resource = FALSE, wlim = c(10, NA))
