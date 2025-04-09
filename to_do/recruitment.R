# This code is intended to compare Mizer and SS different metrics of hake abundance
# by calibrating the Mizer model from biomass and yield values of a base period of years and projecting
# the model into a time frame that is to be compared with SS biomass, recruitment and SSB values
# As Mizer does not have an intrinsic definition of recruitment, abundance levels (number density) for
# different lengths were used as possible proxies.

library(mizer)
library(mizerExperimental)
library(dplyr)
library(reshape)
library(ggplot2)
p <- readRDS("tuned_params.rds")
p <- scaleDownBackground(p, 1/600000)

load("params_data.RData")
load("M_ext.RData")
load("ss_table.RData")

#- Create model taking 2003-2007 avgs. as base for 2008-2023 projection
#- and compare it with the SS model outputs

params<-newSingleSpeciesParams(species_name="Hake",no_w=129,
                               w_max=18310,w_mat=w_mat,lambda=2,h=h,
                               beta=2^b)
species_params(params)$w_mat25<-a*l_mat25^b
species_params(params)$a<-a
species_params(params)$b<-b
species_params(params)$age_mat<-age_mat
params_bygear<-params
gear_params(params_bygear)<-data.frame(
  gear=gear_names,
  species="Hake",catchability=1,
  sel_func="double_sigmoid_length",
  l50=c(26.4,27.5,27.4,14.1,27.0,28.2,47.1,53.4,14.3),
  l50_right=c(32.4,33.2,38.3,20.6,24.6,34.3,57.9,54.4,26.7),
  l25=c(24.1,26.2,25.5,12.6,26.3,26.7,42.5,48.5,12.1),
  l25_right=c(38,43.5,44.1,27,30.4,41.4,66.2,60.8,28.4)
)
initial_effort(params_bygear)<-mean(ss_table$F_val[ss_table$years%in%2003:2007])
nls_mort<-.7629334*w(params)^(-0.1094125) #NLS power law fit
ext_mort(params_bygear)<-array(nls_mort,dim=c(1,129))
ss_biomass<-ss_table$Bio[ss_table$years%in%2003:2007]*1e6 #SS biomass is in tonnes, 1e6 converts to grams
species_params(params_bygear)$biomass_observed<-sum(ss_biomass)/5 #average over 5 years
species_params(params_bygear)$biomass_cutoff<-a*4^b #SS model smallest size is 4 cm

hake_model<- params_bygear |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 

yield_gear<- catch3 %>% filter(year%in%2003:2007) %>% group_by(fleet) %>%
  summarise(avg_yield=sum(amount)/5e-3) %>% slice(c(1,7,2,3,8,5,6,9,4))
gear_params(hake_model)$yield_observed<-yield_gear$avg_yield

#hake_model<-matchYield(hake_model)
hake_model<-setBevertonHolt(hake_model,reproduction_level = 0.001)
hake_model<-steady(hake_model)

effort_time<-ss_table$F_val[ss_table$years%in%2008:2023]
effort_array<-array(replicate(length(gear_names),effort_time),
                    dim=c(16,length(gear_names)),dimnames=list(time=2008:2023,gear=gear_names))

sim<-project(hake_model,t_max=16,effort = effort_array)

# recruits are between ages 0 and 1
# age=1 --> length=13.3652cm

w10<-a*10^b;iw10<-last(which(w(hake_model)<w10))
w11<-a*11^b;iw11<-last(which(w(hake_model)<w11))
w12<-a*12^b;iw12<-last(which(w(hake_model)<w12))
w13<-a*13^b;iw13<-last(which(w(hake_model)<w13))
w20<-a*20^b;iw20<-last(which(w(hake_model)<w20))

rec10<-sim@n[,,iw10]
rec11<-sim@n[,,iw11]
rec12<-sim@n[,,iw12]
rec13<-sim@n[,,iw13]
rec20<-sim@n[,,iw20]

rec_df<-as.data.frame(cbind(rec10,rec11,rec12,rec13,rec20,
                            ss_table$rec_value[ss_table$years%in%2008:2023]))
colnames(rec_df)[6]<-"SS"
rec_df["year"]<-rownames(rec_df)
rec_df_long<-melt(rec_df,id.vars="year")

rec_df_long %>% ggplot(aes(x=year)) +
  geom_line(aes(y=value,group=variable,col=variable),lwd=0.7) +
  scale_y_log10()

#- biomass and SSB mizer vs SS comparison
bio_array<-getBiomass(sim,min_l=4)
bio_ss<-ss_table$Bio[ss_table$years%in%2008:2023]*1e6
bio_df<-data.frame(year=2008:2023,mizer=as.vector(bio_array),ss=bio_ss)

bio_df %>% ggplot(aes(x=year)) +
  geom_line(aes(y=mizer,col="mizer"),lwd=0.7)+
  geom_line(aes(y=ss,col="SS"),lwd=0.7)+
  scale_color_manual(values=c("mizer"="blue","SS"="red"))+
  scale_y_continuous(name="Biomass (g)")+scale_x_continuous(breaks=seq(2008,2023,by=2))+
  theme_bw()
  
ssb_array<-getSSB(sim)
ssb_ss<-ss_table$ssb_val[ss_table$years%in%2008:2023]*1e6
ssb_df<-data.frame(year=2008:2023,mizer=as.vector(ssb_array),ss=ssb_ss)

ssb_df %>% ggplot(aes(x=year)) +
  geom_line(aes(y=mizer,col="mizer"),lwd=0.7)+
  geom_line(aes(y=ss,col="SS"),lwd=0.7)+
  scale_color_manual(values=c("mizer"="blue","SS"="red"))+
  scale_y_continuous(name="SSB (g)")+scale_x_continuous(breaks=seq(2008,2023,by=2))+
  theme_bw()


#- Same for '94-'03 as base period -----
params<-newSingleSpeciesParams(species_name="Hake",no_w=129,
                               w_max=18310,w_mat=w_mat,lambda=2,h=h,
                               beta=2^b)
species_params(params)$w_mat25<-a*l_mat25^b
species_params(params)$a<-a
species_params(params)$b<-b
species_params(params)$age_mat<-age_mat
params_bygear2<-params
gear_params(params_bygear2)<-data.frame(
  gear=gear_names,
  species="Hake",catchability=1,
  sel_func="double_sigmoid_length",
  l50=c(26.4,27.5,27.4,14.1,27.0,28.2,47.1,53.4,14.3),
  l50_right=c(32.4,33.2,38.3,20.6,24.6,34.3,57.9,54.4,26.7),
  l25=c(24.1,26.2,25.5,12.6,26.3,26.7,42.5,48.5,12.1),
  l25_right=c(38,43.5,44.1,27,30.4,41.4,66.2,60.8,28.4)
)
initial_effort(params_bygear2)<-mean(ss_table$F_val[ss_table$years%in%1994:2003])
nls_mort<-.7629334*w(params)^(-0.1094125) #NLS power law fit
ext_mort(params_bygear2)<-array(nls_mort,dim=c(1,129))
ss_biomass<-ss_table$Bio[ss_table$years%in%1994:2003]*1e6 #SS biomass is in tonnes, 1e6 converts to grams
species_params(params_bygear2)$biomass_observed<-sum(ss_biomass)/10 #average over 5 years
species_params(params_bygear2)$biomass_cutoff<-a*4^b #SS model smallest size is 4 cm

hake_model<- params_bygear2 |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 

yield_gear<- catch3 %>% filter(year%in%1994:2003) %>% group_by(fleet) %>%
  summarise(avg_yield=sum(amount)/10e-3) %>% slice(c(1,7,2,3,8,5,6,9,4))
gear_params(hake_model)$yield_observed<-yield_gear$avg_yield

#hake_model<-matchYield(hake_model)
hake_model<-setBevertonHolt(hake_model,reproduction_level = 0.001)
hake_model<-steady(hake_model)

effort_time<-ss_table$F_val[ss_table$years%in%2004:2023]
effort_array<-array(replicate(length(gear_names),effort_time),
                    dim=c(20,length(gear_names)),dimnames=list(time=2004:2023,gear=gear_names))

sim2<-project(hake_model,t_max=20,effort = effort_array)

# recruits are between ages 0 and 1
# age=1 --> length=13.3652cm

w10<-a*10^b;iw10<-last(which(w(hake_model)<w10))
w11<-a*11^b;iw11<-last(which(w(hake_model)<w11))
w12<-a*12^b;iw12<-last(which(w(hake_model)<w12))
w13<-a*13^b;iw13<-last(which(w(hake_model)<w13))

rec10<-sim2@n[,,iw10]
rec11<-sim2@n[,,iw11]
rec12<-sim2@n[,,iw12]
rec13<-sim2@n[,,iw13]

rec_df2<-as.data.frame(cbind(rec10,rec11,rec12,rec13,
                            ss_table$rec_value[ss_table$years%in%2004:2023]))
colnames(rec_df2)[5]<-"SS"
rec_df2["year"]<-rownames(rec_df2)
rec_df_long2<-melt(rec_df2,id.vars="year")

rec_df_long2 %>% ggplot(aes(x=year)) +
  geom_line(aes(y=value,group=variable,col=variable),lwd=0.7) +
  scale_y_log10()

#- biomass and SSB mizer vs SS comparison
bio_array2<-getBiomass(sim2,min_l=4)
bio_ss2<-ss_table$Bio[ss_table$years%in%2004:2023]*1e6
bio_df2<-data.frame(year=2004:2023,mizer=as.vector(bio_array2),ss=bio_ss2)

bio_df2 %>% ggplot(aes(x=year)) +
  geom_line(aes(y=mizer,col="mizer"),lwd=0.7)+
  geom_line(aes(y=ss,col="SS"),lwd=0.7)+
  scale_color_manual(values=c("mizer"="blue","SS"="red"))+
  scale_y_continuous(name="Biomass (g)")+scale_x_continuous(breaks=seq(2004,2023,by=2))

ssb_array2<-getSSB(sim2)
ssb_ss2<-ss_table$ssb_val[ss_table$years%in%2004:2023]*1e6
ssb_df2<-data.frame(year=2004:2023,mizer=as.vector(ssb_array2),ss=ssb_ss2)

ssb_df2 %>% ggplot(aes(x=year)) +
  geom_line(aes(y=mizer,col="mizer"),lwd=0.7)+
  geom_line(aes(y=ss,col="SS"),lwd=0.7)+
  scale_color_manual(values=c("mizer"="blue","SS"="red"))+
  scale_y_continuous(name="SSB (g)")+scale_x_continuous(breaks=seq(2006,2023,by=2))

#- Same for '12-'16 as base period -----
params<-newSingleSpeciesParams(species_name="Hake",no_w=129,
                               w_max=18310,w_mat=w_mat,lambda=2,h=h,
                               beta=2^b)
species_params(params)$w_mat25<-a*l_mat25^b
species_params(params)$a<-a
species_params(params)$b<-b
species_params(params)$age_mat<-age_mat
params_bygear3<-params
gear_params(params_bygear3)<-data.frame(
  gear=gear_names,
  species="Hake",catchability=1,
  sel_func="double_sigmoid_length",
  l50=c(26.4,27.5,27.4,14.1,27.0,28.2,47.1,53.4,14.3),
  l50_right=c(32.4,33.2,38.3,20.6,24.6,34.3,57.9,54.4,26.7),
  l25=c(24.1,26.2,25.5,12.6,26.3,26.7,42.5,48.5,12.1),
  l25_right=c(38,43.5,44.1,27,30.4,41.4,66.2,60.8,28.4)
)
initial_effort(params_bygear3)<-mean(ss_table$F_val[ss_table$years%in%2012:2016])
nls_mort<-.7629334*w(params)^(-0.1094125) #NLS power law fit
ext_mort(params_bygear3)<-array(nls_mort,dim=c(1,129))
ss_biomass<-ss_table$Bio[ss_table$years%in%2012:2016]*1e6 #SS biomass is in tonnes, 1e6 converts to grams
species_params(params_bygear3)$biomass_observed<-sum(ss_biomass)/5 #average over 5 years
species_params(params_bygear3)$biomass_cutoff<-a*4^b #SS model smallest size is 4 cm

hake_model<- params_bygear3 |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() 

yield_gear<- catch3 %>% filter(year%in%2012:2016) %>% group_by(fleet) %>%
  summarise(avg_yield=sum(amount)/5e-3) %>% slice(c(1,7,2,3,8,5,6,9,4))
gear_params(hake_model)$yield_observed<-yield_gear$avg_yield

hake_model<-matchYield(hake_model)
hake_model<-setBevertonHolt(hake_model,reproduction_level = 0.001)
hake_model<-steady(hake_model)

effort_time<-ss_table$F_val[ss_table$years%in%2017:2023]
effort_array<-array(replicate(length(gear_names),effort_time),
                    dim=c(7,length(gear_names)),dimnames=list(time=2017:2023,gear=gear_names))

sim3<-project(hake_model,t_max=7,effort = effort_array)

# recruits are between ages 0 and 1
# age=1 --> length=13.3652cm

w10<-a*10^b;iw10<-last(which(w(hake_model)<w10))
w11<-a*11^b;iw11<-last(which(w(hake_model)<w11))
w12<-a*12^b;iw12<-last(which(w(hake_model)<w12))
w13<-a*13^b;iw13<-last(which(w(hake_model)<w13))

rec10<-sim3@n[,,iw10]
rec11<-sim3@n[,,iw11]
rec12<-sim3@n[,,iw12]
rec13<-sim3@n[,,iw13]

rec_df3<-as.data.frame(cbind(rec10,rec11,rec12,rec13,
                            ss_table$rec_value[ss_table$years%in%2017:2023]))
colnames(rec_df3)[5]<-"SS"
rec_df3["year"]<-rownames(rec_df3)
rec_df_long3<-melt(rec_df3,id.vars="year")

rec_df_long3 %>% ggplot(aes(x=year)) +
  geom_line(aes(y=value,group=variable,col=variable),lwd=0.7) +
  scale_y_log10()

#- biomass and SSB mizer vs SS comparison
bio_array3<-getBiomass(sim3,min_l=4)
bio_ss3<-ss_table$Bio[ss_table$years%in%2017:2023]*1e6
bio_df3<-data.frame(year=2017:2023,mizer=as.vector(bio_array3),ss=bio_ss3)

bio_df3 %>% ggplot(aes(x=year)) +
  geom_line(aes(y=mizer,col="mizer"),lwd=0.7)+
  geom_line(aes(y=ss,col="SS"),lwd=0.7)+
  scale_color_manual(values=c("mizer"="blue","SS"="red"))+
  scale_y_log10(name="Biomass (g)")+scale_x_continuous(breaks=seq(2017,2023))

ssb_array3<-getSSB(sim3)
ssb_ss3<-ss_table$ssb_val[ss_table$years%in%2017:2023]*1e6
ssb_df3<-data.frame(year=2017:2023,mizer=as.vector(ssb_array3),ss=ssb_ss3)

ssb_df3 %>% ggplot(aes(x=year)) +
  geom_line(aes(y=mizer,col="mizer"),lwd=0.7)+
  geom_line(aes(y=ss,col="SS"),lwd=0.7)+
  scale_color_manual(values=c("mizer"="blue","SS"="red"))+
  scale_y_log10(name="SSB (g)")+scale_x_continuous(breaks=seq(2017,2023))

