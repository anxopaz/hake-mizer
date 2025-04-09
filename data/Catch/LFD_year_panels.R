# Construction of the length distribution through the years and for each gear panels contained in
# "LFD_freq_panels.pdf"

# LFD3<-read.csv("D:/Usuarios/dbamio/Hake/hakeData/LFDs 1994-2023.csv",header=T)
# names(LFD3)[6]<-"length"
# LFD3<-LFD3[,-(4:5)]
# LFD3$length<-as.numeric(gsub("len","",LFD3$length))
# LFD3$fleet<-gsub("_fem","",LFD3$fleet)
# LFD3$fleet<-gsub("_mal","",LFD3$fleet)
# LFD3$fleet<-gsub("_ind","",LFD3$fleet)
# LFD3<-LFD3 %>% mutate(weight=a*length^b) %>% mutate(nweight=weight*number)
load("LFD3.RData")

library(dplyr)
library(ggplot2)

#- Preparando los datos
a=0.0037;b=3.168;
by_year<-LFD3 %>% filter(year>1993) %>% filter(!(fleet=="PtSurv"|fleet=="SpSurv"|fleet=="cdSurv")) %>%
  group_by(length,year,fleet) %>% 
  summarise(total=sum(number),.groups="keep") %>% mutate(weight=a*length^b)
gear_list<-split(by_year,f=by_year$fleet)
gear_list<-gear_list[c(1,7,2,3,8,5,6,9,4)] #reordenando los elementos
discNA<-data.frame(length=1,year=c(1995,1996,1998,2001,2002),fleet="disc",total=NA,weight=0.0037)
gear_list[[9]]<-rbind(gear_list[[9]],discNA)
#- Paleta de colores
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols<-gg_color_hue(9)

gear_list2<-list()
for (i in 1:length(gear_list)) {
gear_list2[[i]] <- inner_join(gear_list[[i]] %>% 
                        group_by(year) %>% 
                        summarise(sumtotal = sum(total)),
                        gear_list[[i]] %>% 
                        group_by(length,year),
                      by = "year") %>%
  mutate(prop = total/sumtotal)
names(gear_list2)[[i]]<-gear_list2[[i]]$fleet[1]
}

#### Todos los paneles en un mismo archivo -----
#- Frecuencia relativa
pdf("D:/Usuarios/dbamio/Hake/mizer-hake/imgs/LFD_freq_panels.pdf", onefile = TRUE)
for(i in 1:length(gear_list2)){
  dat <- gear_list2[[i]]
  var <- names(gear_list2)[i]
  l.plot <- ggplot(dat,aes(x=length)) +
    geom_area(aes(y=prop),fill=cols[i],alpha=.6)+
    geom_line(aes(y=prop),lwd=.2)+
    scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7),legend.position = "none")+
    facet_wrap(~year,nrow=6,ncol=5,drop=F,dir="v")+labs(title=var)+xlab("length (cm)")+
    ylab("proportion")
  print(l.plot)
  w.plot <- ggplot(dat,aes(x=weight)) +
    geom_area(aes(y=prop),fill=cols[i],alpha=.6)+
    geom_line(aes(y=prop),lwd=.2)+
    scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),legend.position = "none")+
    facet_wrap(~year,nrow=6,ncol=5,drop=F,dir="v")+labs(title=var)+xlab("weight (g)")+
    ylab("proportion")
  print(w.plot)
}
dev.off()

#- Densidad
pdf("D:/Usuarios/dbamio/Hake/mizer-hake/imgs/LFD_density_panels.pdf", onefile = TRUE)
for(i in 1:length(gear_list)){
  dat <- gear_list[[i]]
  var <- names(gear_list)[i]
  l.plot <- ggplot(dat, aes(x=length)) +
    geom_density(aes(weight=total),fill=cols[i],alpha=.6,adjust=0.2,position="identity")+
    scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
    facet_wrap(~year,nrow=6,ncol=5,drop=F,dir="v")+labs(title=var)+xlab("length (cm)")
  print(l.plot)
  w.plot <- ggplot(dat,aes(x=weight)) +
    geom_density(aes(weight=total),fill=cols[i],alpha=.6,adjust=0.2,position="identity")+
    scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
    facet_wrap(~year,nrow=6,ncol=5,drop=F,dir="v")+labs(title=var)+xlab("weight (g)")
  print(w.plot)
}
dev.off()



#### Paneles de gráficas de frecuencia por año -----
# Por longitud
gear_list2[[9]] %>% ggplot(aes(x=length)) +
  geom_area(aes(y=prop),fill=cols[9],alpha=.6)+
  geom_line(aes(y=prop),lwd=.2)+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7),legend.position = "none")+
  facet_wrap(~year,nrow=6,ncol=5,drop=F,dir="v")+labs(title=names(gear_list)[9])+xlab("length (cm)")+
  ylab("proportion")
# Por peso
gear_list2[[1]] %>% ggplot(aes(x=weight)) +
  geom_area(aes(y=prop,fill=cols[1],alpha=.6))+
  geom_line(aes(y=prop),lwd=.2)+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),legend.position = "none")+
  facet_wrap(~year,nrow=6,ncol=5,drop=F,dir="v")+labs(title=names(gear_list)[1])+xlab("weight (g)")+
  ylab("proportion")


#### Paneles de gráficas de densidad de talla por año -------
#### 1. Gear: Art
#- Por longitud
gear_list[[1]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[1],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: Art")+xlab("length (cm)")
#- Por peso
gear_list[[1]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[1],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: Art")+xlab("weight (g)")

#### 2. Gear: ptArt
#- Por longitud
gear_list[[2]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[2],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: ptArt")+xlab("length (cm)")
#- Por peso
gear_list[[2]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[2],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: ptArt")+xlab("weight (g)")

#### 3. Gear: bakka
#- Por longitud
gear_list[[3]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[3],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: bakka")+xlab("length (cm)")
#- Por peso
gear_list[[3]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[3],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: bakka")+xlab("weight (g)")

#### 4. Gear: cdTrw
#- Por longitud
gear_list[[4]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[4],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: cdTrw")+xlab("length (cm)")
#- Por peso
gear_list[[4]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[4],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: cdTrw")+xlab("weight (g)")

#### 5. Gear: ptTrw
#- Por longitud
gear_list[[5]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[5],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: ptTrw")+xlab("length (cm)")
#- Por peso
gear_list[[5]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[5],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: ptTrw")+xlab("weight (g)")

#### 6. Gear: pair
#- Por longitud
gear_list[[6]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[6],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: pair")+xlab("length (cm)")
#- Por peso
gear_list[[6]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[6],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: pair")+xlab("weight (g)")

#### 7. Gear: palangre
#- Por longitud
gear_list[[7]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[7],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: palangre")+xlab("length (cm)")
#- Por peso
gear_list[[7]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[7],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: palangre")+xlab("weight (g)")

#### 8. Gear: vol
#- Por longitud
gear_list[[8]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[8],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: vol")+xlab("length (cm)")
#- Por peso
gear_list[[8]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[8],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: vol")+xlab("weight (g)")


#### 9. Gear: disc 
discNA<-data.frame(length=1,year=c(1995,1996,1998,2001,2002),fleet="disc",total=NA,weight=0.0037)
gear_list[[9]]<-rbind(gear_list[[9]],discNA)
#- Por longitud
gear_list[[9]] %>% ggplot(aes(x=length)) +
  geom_density(aes(weight=total),fill=cols[9],alpha=.6,adjust=0.2,position="identity")+
  scale_x_continuous(breaks=c(0,20,40,60,80,100,120))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: disc")+xlab("length (cm)")
#- Por peso
gear_list[[9]] %>% ggplot(aes(x=weight)) +
  geom_density(aes(weight=total),fill=cols[9],alpha=.6,adjust=0.2,position="identity")+
  scale_x_log10(breaks=c(0,1,10,100,1000,10000),labels=scales::scientific)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~year,nrow=6,ncol=5,drop=F)+labs(title="Gear: disc")+xlab("weight (g)")


