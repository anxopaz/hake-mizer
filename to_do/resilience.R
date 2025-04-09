# This script was aimed as a study of the relationship between fishing effort (F) and different abundance
# metrics such as biomass, recruitment and SSB in the SS data.
# This relationship is crucial when trying to adjust the reproduction parameters in mizer

library(mizer)
library(mizerExperimental)
library(dplyr)
library(reshape)
library(ggplot2)
p <- readRDS("tuned_params.rds")
p <- scaleDownBackground(p, 1/600000)

#######
# Changing the reproduction_level and interaction_resource parameters in the mizer model
# to check for their impact in the size spectrum

p<-setBevertonHolt(p,reproduction_level=0.9)
species_params(p)$interaction_resource<-1
sim_p<-project(p,t_max=25)
animateSpectra(sim_p,power=2)

p1<-setBevertonHolt(p,reproduction_level = 0.9)
species_params(p1)$interaction_resource<-0.8
sim_p1<-project(p1,t_max=25)
animateSpectra(sim_p1,power=2)

p0<-setBevertonHolt(p,reproduction_level = 0.01)
species_params(p0)$interaction_resource<-0.8
sim_p0<-project(p0,t_max=25)
animateSpectra(sim_p1,power=2)

plotSpectra2(sim_p1,sim_p0,power=2)

#########
# Growth curve plot

plotGrowthCurves(p)
growth_curve<-getGrowthCurves(p,max_age=14)
growth_curve_conv<-(growth_curve/species_params(p)$a)^(1/species_params(p)$b) #conversion from weight to length
growth_curve_ages<-as.numeric(names(as.data.frame(growth_curve_conv)))
growth_curve_lengths<-c(growth_curve_conv)
growth_curve_dataframe<-data.frame(age=growth_curve_ages,length=growth_curve_lengths,weight=c(growth_curve))
plot(growth_curve_ages,growth_curve_lengths)

growth_curve_dataframe %>% ggplot(aes(x=age)) +
  geom_line(aes(y=length),lwd=0.7)+theme_bw()+
  scale_y_continuous(limits=c(0,120),breaks=c(0,20,40,60,80,100,120))

#########
# SS data study of cross correlation between F and biomass/recruitment/SSB time series

load("ss_table.RData")

scaleFactor <- max(ss_table[30:64,]$F_val,na.rm=T) / max(ss_table[30:64,]$Bio,na.rm=T)
scaleFactor2 <- max(ss_table[30:64,]$F_val,na.rm=T) / max(ss_table[30:64,]$rec_value,na.rm=T)
scaleFactor3 <- max(ss_table[30:64,]$F_val,na.rm=T) / max(ss_table[30:64,]$ssb_val,na.rm=T)

ss_table[30:64,] %>% ggplot(aes(x=years)) +
  geom_line(aes(y=F_val),col="blue",lwd=0.7) +
  geom_line(aes(y=Bio*scaleFactor),col="red",lwd=0.7) +
  scale_y_continuous(name="F",sec.axis = sec_axis(~./scaleFactor,name="Biomass"))+
  #geom_line(aes(y=rec_value*scaleFactor2),col="red",lwd=0.7)+
  #scale_y_continuous(name="F",sec.axis = sec_axis(~./scaleFactor2,name="Recruitment"))+
  #geom_line(aes(y=ssb_val*scaleFactor3),col="red",lwd=0.7)+
  #scale_y_continuous(name="F",sec.axis = sec_axis(~./scaleFactor3,name="Spawning stock biomass"))+
  theme(
    axis.title.y.left=element_text(color="blue"),
    axis.text.y.left=element_text(color="blue"),
    axis.title.y.right=element_text(color="red"),
    axis.text.y.right=element_text(color="red")
  )

# F - biomass cross-correlation
ccf(ss_table[34:64,]$F_val,ss_table[34:64,]$Bio,type="correlation")
# Lag of 3 years from F to biomass

model<-lm(ss_table[37:64,]$Bio~ss_table[34:61,]$F_val)
summary(model)

plot(ss_table[37:64,]$F_val,ss_table[34:61,]$Bio,type="p",
     ylab="Biomass",xlab="F",pch=16,cex=1)
abline(model,lwd=2,col="red")

# F - Spawning stock biomass correlation: also 3 year lag
ccf(ss_table[34:64,]$F_val,ss_table[34:64,]$ssb_val,type="correlation")
model<-lm(ss_table[37:64,]$ssb_val~ss_table[34:61,]$F_val)
summary(model)
plot(ss_table[37:64,]$F_val,ss_table[34:61,]$ssb_val,type="p",
     ylab="Spawning stock biomass",xlab="F",pch=16,cex=1)
abline(model,lwd=2,col="red")

# F - catch: 5 year lag
ccf(ss_table[34:64,]$F_val,ss_table[34:64,]$catch,type="correlation",na.action = na.pass)
model<-lm(ss_table[39:64,]$catch~ss_table[34:59,]$F_val)
summary(model)
  plot(ss_table[39:64,]$F_val,ss_table[34:59,]$catch,type="p",
     ylab="Catch",xlab="F",pch=16,cex=1)
abline(model,lwd=2,col="red")

# CPUE as catch / F (?)
ss_table$cpue<-ss_table$catch/ss_table$F_val

ss_table[30:64,] %>% ggplot(aes(x=years)) +
  geom_line(aes(y=cpue,color="CPUE"),lwd=1) +
  geom_line(aes(y=Bio,col="Biomass",),lwd=1) +
  scale_y_continuous(limits=c(0,60000))+labs(y="",x="Year",colour="Legend")

ccf(ss_table[34:64,]$catch,ss_table[34:64,]$cpue,type="correlation",na.action = na.pass)
model<-lm(ss_table[34:63,]$cpue~ss_table[35:64,]$catch)
summary(model)
plot(ss_table[35:64,]$catch,ss_table[34:63,]$cpue,type="p",
     ylab="CPUE",xlab="Catch",pch=16,cex=1)
abline(model,lwd=2,col="red")
