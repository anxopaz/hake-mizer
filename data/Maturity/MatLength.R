
rm(list=ls())

fdata <- read.csv("./data/Maturity/mat_data.csv", header = T, check.names = FALSE, sep = ";" ,dec = ".", stringsAsFactors = F)
fdata$year_mat <- as.factor(fdata$year_mat)
fdata$sex <- as.factor(fdata$sex)
fdata$lab <- as.factor(fdata$lab)
fdata$month <- as.factor(fdata$month)

last <- unique(fdata$year)

MatSize <- matrix( NA, 2, 2, dimnames = list( c('L50','L25'), c('Males','Females')))

for( i in c(1,2)){
  
  data <- subset( fdata, fdata$sex == i)
  head(data)
  
  ind <- is.na(data$mat)
  ind <- which(ind==TRUE)
  data <- data[-ind,]
  
  cutoff_lengths <- c(seq(min(data$lt),19,by=1),seq(from=20, to=40, by=1), seq(from=42, to=70, by=2),seq(71,max(data$lt),by=4))
  data$bin <- cut(data$lt, cutoff_lengths, labels = cutoff_lengths[-1])
  
  aux <- subset(data,data$lt<21)[,c(3,5)]
  ind <- which(aux$mat==1)
  
  data$mat[data$lt < 21 ] <- 0
  
  NLbins<-c(seq(from=20, to=40, by=1),seq(from=42, to=70, by=2)) # Desired bins (SS model) 67
  l_b <- length(NLbins)
  
  len <- data$lt
  l_len <- length(len); aux <- rep(0,l_len)
  
  years <- (min(as.numeric(as.character(data$year_mat))):max(as.numeric(as.character(data$year_mat))))
  
  data_ieo <- subset(data,data$lab=="ieo")
  data_ipma <- subset(data,data$lab=="ipma")
  data <- rbind(data_ieo,data_ipma)
  
  ind_ieo <- which(data$lab=="ieo"); ind_ipma <- which(data$lab=="ipma")
  len <- length(data$lab)
  
  len_ieo <- length(ind_ieo)
  len_ipma <- length(ind_ipma)
  
  YCombined <- matrix(NA, nrow = len, ncol = 2)
  YCombined[1:len_ieo, 1] <- (data$mat[ind_ieo])
  YCombined[(len_ieo+1):(len_ipma+len_ieo), 2] <- (data$mat[ind_ipma])
  
  data$Gyear_mat <- as.character(data$year_mat)
  ind <- which(as.numeric(as.character(data$year_mat))<2001)
  data$Gyear_mat[ind] <- "1980_2000"
  
  ind <- which(as.numeric(as.character(data$year_mat))>2016)
  data$Gyear_mat[ind] <- "2017-2019"
  
  data$Gyear_mat <- as.factor(data$Gyear_mat)
  
  f3 <-  YCombined ~ 1 + lt + f( Gyear_mat, model = "iid")
  
  I3 <- inla( f3, control.compute = list(config=TRUE, dic = TRUE, cpo=TRUE),
              family = c("binomial","binomial"), data = data, 
              control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads = 1)
  
  intercept <- I3$summary.fixed[1, 1]; coef_lt <- I3$summary.fixed[2, 1]
 
  MatSize['L50',i] <- - intercept/coef_lt
  MatSize['L25',i] <- (log(0.25 / (1 - 0.25)) - intercept) / coef_lt
  
}
  
  
save( MatSize, file = './input/MatLength.RData')
  
