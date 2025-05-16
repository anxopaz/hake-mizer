
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~   Hake's MIZER model with cannibalism  ~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rm(list=ls())

library(mizer)
library(mizerExperimental)
load( './output/hake_model.RData')
load( './data/Diet.RData')

plotSpectra( hake_mizer, power = 2)

cannibal <- hake_mizer2 <- hake_mizer

ext_mort( hake_mizer2)[1:60] <- 0

# Turn on cannibalism ----------------------------

pcann <- as.numeric(by_prey[which(by_prey$Prey=='M.merluccius'),2])

interaction_matrix( cannibal)[] <- pcann
interaction_matrix( hake_mizer2)[] <- pcann

getPredMort( cannibal) - getPredMort( hake_mizer2)

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

pre_obj1 <- prefit( model = cannibal, catch = corLFD, dl = 1, yield_lambda = 1e7)
pre_obj2 <- prefit( model = cannibal2, catch = corLFD, dl = 1, yield_lambda = 1e7)

data_list1 <- pre_obj1$data_list
pars1 <- pre_obj1$pars

data_list2 <- pre_obj2$data_list
pars2 <- pre_obj2$pars

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
dyn.load(dynlib("./TMB/fit"))

obj1 <- MakeADFun(data = data_list1, parameters = pars1, DLL = "fit")
obj2 <- MakeADFun(data = data_list2, parameters = pars2, DLL = "fit")

optim_result1 <- nlminb(obj1$par, obj1$fn, obj1$gr, control = list(trace = 1, eval.max = 10000, iter.max = 10000))
optim_result2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control = list(trace = 1, eval.max = 10000, iter.max = 10000))

cannibal <- update_params( cannibal, optim_result1$par, data_list1$min_len, data_list1$max_len)
cannibal2 <- update_params( cannibal2, optim_result2$par, data_list2$min_len, data_list2$max_len)


# tuneParams( cannibal, catch = catch_lengths)
# tuneParams( cannibal2, catch = catch_lengths)

# cannibal <- readParams("./output/cannibal.rds")
# cannibal2 <- readParams("./output/cannibal2.rds")

# Comparison -----------------

plotDiet( hake_mizer)
plotDiet( cannibal)
plotDiet( cannibal2)

p1 <- plotDeath( hake_mizer) + theme_bw(); p1
p2 <- plotDeath( cannibal) + theme_bw(); p2
p3 <- plotDeath( cannibal2) + theme_bw(); p3

ggpubr::ggarrange(p1, p2, p3, nrow=1, common.legend = TRUE, legend="bottom")

plotSpectra2(hake_mizer, name1 = "Base Model",
             cannibal2, name2 = "With cannibalism",
             power = 2, resource = FALSE, wlim = c(10, NA)) + theme_bw()

plotSpectra2(hake_mizer, name1 = "Base Model",
             cannibal, name2 = "With cannibalism - NM",
             power = 2, resource = FALSE, wlim = c(10, NA)) + theme_bw()

plotSpectra2(cannibal2, name1 = "With cannibalism",
             cannibal, name2 = "With cannibalism - NM",
             power = 2, resource = FALSE, wlim = c(10, NA)) + theme_bw()
