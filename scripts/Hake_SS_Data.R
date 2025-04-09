######## WGBIE24 southern hake data from ICES ##########

rm( list=ls())

assessment <- icesSAG::getSAG( 'hke.27.8c9a', 2024)

assessment <- assessment[,-c(14:ncol(assessment))]

for ( i in 1:nrow(assessment)) 
  assessment$catches[i] <- assessment$landings[i] + ifelse(is.na(assessment$discards[i]),0,assessment$discards[i])

save( assessment, file = './input/Hake_SS_Data.RData')

