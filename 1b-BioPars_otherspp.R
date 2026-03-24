# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~   Other spps biological parameters   ~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

rm(list=ls())

library(dplyr)
library(ggplot2)
library(plotly)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(rfishbase)
library(tidyr)

source( './scripts/aux_functions.R')


# Other spps ------------------------------

## Length - Weight relationship ------------------------

spp_pars <- data.frame( readxl::read_excel( "./data/Demersales/spp_pars.xlsx"))
spp_pars

FBdata <- list()

spp_pars$l_max <- spp_pars$w_max <- 
  spp_pars$age_mat <- spp_pars$age_mat2 <- 
  spp_pars$w_mat <- spp_pars$l_mat <- 
  spp_pars$l_inf <- spp_pars$kvb <- spp_pars$al0 <- 
  spp_pars$M <- NA

for(i in 1:nrow(spp_pars)){
  
  isp <- spp_pars$esp[i]
  
  ia <- spp_pars[i,'a']
  ib <- spp_pars[i,'b']
  
  i_gr <- popgrowth(isp)
  spp_pars[i,'l_inf'] <- mean( i_gr$Loo, na.rm = T)
  spp_pars[i,'kvb'] <- mean( i_gr$K, na.rm = T)
  spp_pars[i,'al0'] <- mean( i_gr$to, na.rm = T)
  
  i_max <- popchar(isp)
  spp_pars[i,'l_max'] <- mean( i_max$Lmax, na.rm = T)
  spp_pars[i,'w_max'] <- lwf( spp_pars[i,'l_max'], ia, ib)
  
  i_pred <- estimate(isp)
  
  i_mat <- maturity(isp)
  spp_pars[i,'l_mat'] <- mean( i_mat$Lm, na.rm = T)
  spp_pars[i,'w_mat'] <- lwf( spp_pars[i,'l_mat'], ia, ib)
  spp_pars[i,'age_mat'] <- laf( 
    spp_pars[i,'l_mat'], spp_pars[i,'l_inf'], spp_pars[i,'kvb'], spp_pars[i,'al0'])
  spp_pars[i,'age_mat2'] <- mean( i_mat$tm, na.rm = T)
  
  spp_pars[i,'M'] <- mean( i_gr$M, na.rm = T)
  
  FBdata[[isp]] <- list( growth = i_gr, max = i_max, pred = i_pred, maturity = i_mat)
  
}

save( spp_pars, FBdata, file = './input/other_spp.RData')
# load('./input/other_spp.RData')


spp_pars2 <- data.frame( readxl::read_excel( "./data/Demersales/spp_pars_wiki.xlsx"))
spp_pars2 <- spp_pars2 %>% mutate( l_max = l_inf)


load('./data/Demersales/not_included/other_spp.RData')


other_spp <- list( short = spp_pars, short_fb = spp_pars2, long = spps)

other_spp_plots <- list( length_weight = list(), von_bertalanffy = list())


for(i in 1:3){
  
  i_n <- names(other_spp)[[i]]
  
  LW_df <- other_spp[[i]] %>%
    mutate(length = purrr::map(l_max, ~ seq(1, .x, length.out = 200))) %>%
    unnest(length) %>%
    mutate(weight = lwf(length, a, b))
  
  
  VB_df <- other_spp[[i]] %>%
    mutate(age = purrr::map(l_max, ~ seq(0, 10, length.out = 200))) %>% 
    unnest(age) %>%
    mutate(length = alf(age, l_inf, kvb, al0))
  
  
  lwp <- LW_df %>%
    ggplot(aes(x = length, y = weight, color = common)) +
    geom_line(linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0("Length–Weight relationship (",i_n,")"),
      x = "Length (cm)",
      y = "Weight (g)",
      color = "Species"
    )
  
  lwp2 <- ggplotly(lwp)
  
  for(i in seq_along(lwp2$x$data)){ lwp2$x$data[[i]]$visible <- "legendonly"}
  
  other_spp_plots$length_weight[[i_n]] <- lwp2
  
  vbp <- VB_df %>%
    ggplot(aes(x = age, y = length, color = common)) +
    geom_line(linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0("Von Bertalanffy growth (",i_n,")"),
      x = "Age (years)",
      y = "Length (cm)",
      color = "Species"
    )
  
  vbp2<- ggplotly(vbp)
  
  for(i in seq_along(vbp2$x$data)){ vbp2$x$data[[i]]$visible <- "legendonly"}
  
  other_spp_plots$von_bertalanffy[[i_n]] <- vbp2
  
}



## ICES for biomass?? --------------

# load( './input/other_spp_rmd.RData')

library(icesSAG)

eng_enc_8 <- getSAG(stock="ane.27.8",year = 2025); tibble(eng_enc_8)
eng_enc_9a <- getSAG(stock="ane.27.9a",year = 2025); tibble(eng_enc_9a)
sar_pil_8c9a <- getSAG(stock="pil.27.8c9a",year = 2025); tibble(sar_pil_8c9a)
tra_tra_9a <- getSAG(stock="hom.27.9a",year = 2025); tibble(tra_tra_9a)
mic_pou_191214 <- getSAG(stock="whb.27.1-91214",year = 2025); tibble(mic_pou_191214)
sco_sco <- getSAG(stock="mac.27.nea",year = 2025); tibble(sco_sco)
meg_78abd <- getSAG("meg.27.7b-k8abd", year = 2025); tibble(meg_78abd)
meg_8c9a <- getSAG("meg.27.8c9a", year = 2025); tibble(meg_8c9a)
ldb_78abd <- getSAG("ldb.27.7b-k8abd", year = 2025); tibble(ldb_78abd)
ldb_8c9a <- getSAG("ldb.27.8c9a", year = 2025); tibble(ldb_8c9a)


engr_enc1 <- tibble(eng_enc_8)|> select(Year, SSB)
engr_enc2 <- tibble(eng_enc_9a)|> select(Year, SSB)

engr_enc <- bind_rows( engr_enc1 |> mutate(stock = "a"), engr_enc2 |> mutate(stock = "b")) |>
  group_by(Year) |> filter(n_distinct(stock) == 2) |>
  summarise(SSB = sum(SSB, na.rm = TRUE), .groups = "drop")


# lep_whi1 <- tibble(meg_78abd)|> select(Year, SSB)
# lep_whi2 <- tibble(meg_8c9a)|> select(Year, SSB)
# lep_bos1 <- tibble(ldb_78abd)|> select(Year, SSB)
# lep_bos2 <- tibble(ldb_8c9a)|> select(Year, SSB)
# 
# lepi_whi <- bind_rows(lep_whi1, lep_whi2) |> group_by(Year) |>
#   summarise(SSB = sum(SSB, na.rm = TRUE), .groups = "drop")
# 
# lepi_bos <- bind_rows(lep_bos1, lep_bos2) |> group_by(Year) |>
#   summarise(SSB = sum(SSB, na.rm = TRUE), .groups = "drop")

lepi_whi <- tibble(meg_8c9a)|> select(Year, SSB)

lepi_bos <- tibble(ldb_8c9a)|> select(Year, SSB)

sard_pil <- tibble(sar_pil_8c9a)|> select(Year, SSB)

trac_tra <- tibble(tra_tra_9a)|> select(Year, SSB)|> mutate(SSB=SSB*2.2)

micr_pou <- tibble(mic_pou_191214)|> select(Year, SSB)|> mutate(SSB=SSB*0.04)

scom_sco <- tibble(sco_sco)|> select(Year, SSB)|> mutate(SSB=SSB*0.04)



other_spp_ices <- list( Eng_encr = engr_enc, Sar_pilc = sard_pil, Tra_trac = trac_tra,
  Mic_pout = micr_pou, Sco_scom = scom_sco, Lep_bosc = lepi_bos, Lep_whif = lepi_whi)

other_spp_ices2 <- list( eng_enc_8 = eng_enc_8, eng_enc_9a = eng_enc_9a, sar_pil_8c9a = sar_pil_8c9a, 
                         tra_tra_9a = tra_tra_9a, mic_pou_91214 = mic_pou_91214, sco_sco = sco_sco)

save( other_spp, other_spp_ices, other_spp_plots, file = './input/other_spp_rmd.RData')

