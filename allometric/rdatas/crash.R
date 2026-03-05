load('./allometric/prefit.RData')

# data <- list( counts = data$counts, bin_index = data$bin_index, f_index=data$f_index, sel_type=gp$sel_type, 
#               coeff_fj=data$coeff_fj,coeff_fj1=data$coeff_fj1, dw=data$dw, w=data$w, l=data$l, 
#               yield=data$yield, biomass = data$biomass, biomass_cutoff_idx = data$biomass_cutoff_idx,
#               growth = data$growth, w_mat=data$w_mat, d=data$d, yield_lambda=1)

# initial_params$log_l50_right_offset[c(5,9)] <- c(1,1)
# initial_params$log_ratio_right[c(5,9)] <- c(-1,-1)


if ("objective_function4" %in% names(getLoadedDLLs())) dyn.unload( TMB::dynlib("./allometric/objective_function4"))
TMB::compile("./allometric/objective_function4.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
dyn.load( TMB::dynlib("./allometric/objective_function4"))

# data$counts <- round( data$counts)
obj <- TMB::MakeADFun( data = data, parameters = initial_params, DLL = "objective_function4", silent = FALSE, debug = TRUE)

optim_result <- nlminb( obj$par, obj$fn, obj$gr, lower = lower_bounds, upper = upper_bounds,
                        control = list(trace = 0, eval.max = 1000, iter.max = 1000))

optim_result$par[1:9]
initial_params$logit_l50

# rep <- obj$report()
# summary(rep$probs)
# sum(rep$probs)
# 
# grad <- obj$gr(obj$par)
# which(!is.finite(grad))
# grad[!is.finite(grad)]

# gdb R
# run
# source("D:/Usuarios/apaz/Documents/hake-mizer/allometric/crash.R")
# bt

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

newpars <- optim_result$par

gplist <- list()
gpnames <- c( 'logit_l50', 'log_ratio_left', 'log_l50_right_offset', 'log_ratio_right',
              'log_catchability')

for (i in gpnames) gplist[[i]] <- as.numeric(newpars[grep(i, names(newpars))])

lmin <- min(data$l)
lmax <- max(data$l)

l50 <- lmin + (lmax - lmin) * plogis(gplist$logit_l50)
l25 <- l50 * (1 - exp(gplist$log_ratio_left))
l50_right <- l50 + exp(gplist$log_l50_right_offset)
l25_right <- l50_right * (1 + exp(gplist$log_ratio_right))
catchability <- exp(gplist$log_catchability)

gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)

hake_model@gear_params[,'l50'] <- gp_res$l50
hake_model@gear_params[,'l25'] <- gp_res$l25
hake_model@gear_params[,'l50_right'] <- gp_res$l50_right
hake_model@gear_params[,'l25_right'] <- gp_res$l25_right
hake_model@gear_params[,'catchability'] <- gp_res$catchability

model@species_params <- sp

plot_lfd( hake_model, corLFD)
plot_lfd_gear( hake_model, corLFD)
if(plot == T){
  ggsave( paste0( plot_dir, 'LFD_gear.jpg'), width = 9, height = 7)
}