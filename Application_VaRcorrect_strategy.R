library(PerformanceAnalytics)
library(xts)
library(matrixStats)
library(xtable)


portp <- function(portre){
  N = length(portre)
  
  # 208 day per year, given 4 days a year
  annulre = (mean(portre/100)+1)^208-1
  sr = (mean(portre)/100)/sd(portre) * sqrt(208)

  cumulative_returns <- cumprod(1 + portre/100) 
  running_max <- cummax(cumulative_returns)
  drawdowns <- (cumulative_returns - running_max) / running_max
  max_drawdown <- min(drawdowns) 

  dd = max_drawdown
  return(list(annulre, sr, dd))
}



setwd('your working directory')
euspot = read.csv('euspot_24.csv',sep=",",header=FALSE)
ukspot = read.csv('ukspot_24.csv',sep=",",header=FALSE)
jpspot = read.csv('jpspot_24.csv',sep=",",header=FALSE)

intra.return <- function(yd){
  grid_point=nrow(yd)
  N=ncol(yd)
  idr_re=matrix(0,grid_point,N)
  cidr_re=matrix(0,grid_point,N)
  ocidr_re=matrix(0,grid_point,N-1)
  for(j in c(2:grid_point)){
    for(i in c(1:N)){
      idr_re[j,i]=log(yd[j,i])-log(yd[j-1,i])
    }
  }
  for(j in c(1:grid_point)){
    for(i in c(1:N)){
      cidr_re[j,i]=log(yd[j,i])-log(yd[1,i])
    }
  }
  for(j in c(1:grid_point)){
    for(i in c(2:N)){
      ocidr_re[j,i-1]=log(yd[j,i])-log(yd[grid_point,i-1])
    }
  }
  return(list(idr = idr_re*100, cidr = cidr_re, ocidr = ocidr_re*100 ))
}

eu_re = intra.return(euspot)$ocidr
uk_re = intra.return(ukspot)$ocidr
jp_re = intra.return(jpspot)$ocidr

win=102
eu_true = eu_re[,(1274-674-win+1):1274]
uk_true = uk_re[,(1274-674-win+1):1274]
jp_true = jp_re[,(1274-674-win+1):1274]





# Application intraday trading strategy
#################### long
##############
### eu
##############
#for(country in c("eu","jp","uk")){ #Loop through different currencies
country = "eu"
true_re = eu_true
setwd('your working directory')

load(file.path(paste0(country,"_forecast_intraday_g.RData"))) # g or x
profit_bench<-profit_tfpca<-profit_dfpca<-profit_lfpca<-profit_mfpca<-matrix(nrow=674, ncol=1, byrow=TRUE) #Set up variables to hold daily trading strategy profits


for(day in 1:674){
  # benchmark - functional mean
  roll_fmean = rowMeans(true_re[,day:(day+win-1)]) # historical functional mean
  
  iBuy <- which.min(roll_fmean[2:287])

  iSell <- iBuy-1+which.max(roll_fmean[iBuy:length(roll_fmean)])

  if((roll_fmean[iSell]-roll_fmean[iBuy])>=mean(roll_fmean) & iBuy!=288){# 287
    Ret0 <- true_re[iSell,day+win] - true_re[iBuy,day+win]
    fee <- 0.0003/100
    Ret <- Ret0 - 2*fee - fee*Ret0
  }else{
    Ret <-0
  }

  # MFPCA - VaR
  if(any(true_re[2:287,day+win]<=g_var_mfpca_eu01[2:287,day])==TRUE){
    iSell_mfpca = which(true_re[(iBuy):287,day+win]<g_var_mfpca_eu01[(iBuy):287,day])[1]
    if(is.na(iSell_mfpca)==TRUE){iSell_mfpca = iSell}
    Ret0_mfpca <- true_re[iSell_mfpca,day+win] - true_re[iBuy,day+win]
    Ret_mfpca <- Ret0_mfpca - 2*fee - fee*Ret0_mfpca
    if(Ret_mfpca<0){Ret_mfpca = 0}
  }else{
    Ret_mfpca <-Ret
  }

  profit_bench[day]<-Ret
  Ret_mfpca[is.na(Ret_mfpca)] <- 0
  profit_mfpca[day]<-Ret_mfpca
}

tablep1 = matrix(NA,2,3)
tablep1[1,1] = portp(profit_bench)[[1]]
tablep1[2,1] = portp(profit_mfpca)[[1]]

tablep1[1,2] = portp(profit_bench)[[2]]
tablep1[2,2] = portp(profit_mfpca)[[2]]

tablep1[1,3] = portp(profit_bench)[[3]]
tablep1[2,3] = portp(profit_mfpca)[[3]]

profit_fgmfpca=profit_mfpca

#### fgarch-x

load(file.path(paste0(country,"_forecast_intraday_x.RData"))) # g or x
profit_bench<-profit_tfpca<-profit_dfpca<-profit_lfpca<-profit_mfpca<-matrix(nrow=674, ncol=1, byrow=TRUE) #Set up variables to hold daily trading strategy profits

for(day in 1:674){
  # benchmark - functional mean
  roll_fmean = rowMeans(true_re[,day:(day+win-1)]) # historical functional mean
  iBuy <- which.min(roll_fmean[2:287])
  iSell <- iBuy-1+which.max(roll_fmean[iBuy:length(roll_fmean)])

  if((roll_fmean[iSell]-roll_fmean[iBuy])>=mean(roll_fmean) & iBuy!=288){# 287
    Ret0 <- true_re[iSell,day+win] - true_re[iBuy,day+win]
    fee <- 0.0003/100
    Ret <- Ret0 - 2*fee - fee*Ret0
  }else{
    Ret <-0
  }

  # MFPCA - VaR
  if(any(true_re[2:287,day+win]<=x_var_mfpca_eu01[2:287,day])==TRUE){
    iSell_mfpca = which(true_re[(iBuy+1):287,day+win]<x_var_mfpca_eu01[(iBuy+1):287,day])[1]
    if(is.na(iSell_mfpca)==TRUE){iSell_mfpca = iSell}
    Ret0_mfpca <- true_re[iSell_mfpca,day+win] - true_re[iBuy,day+win]
    Ret_mfpca <- Ret0_mfpca - 2*fee - fee*Ret0_mfpca
    if(Ret_mfpca<0){Ret_mfpca = 0}
  }else{
    Ret_mfpca <-Ret
  } 

  profit_bench[day]<-Ret
  Ret_mfpca[is.na(Ret_mfpca)] <- 0
  profit_mfpca[day]<-Ret_mfpca
}

tablep = matrix(NA,2,3)
tablep[1,1] = portp(profit_bench)[[1]]
tablep[2,1] = portp(profit_mfpca)[[1]]

tablep[1,2] = portp(profit_bench)[[2]]
tablep[2,2] = portp(profit_mfpca)[[2]]

tablep[1,3] = portp(profit_bench)[[3]]
tablep[2,3] = portp(profit_mfpca)[[3]]

xtable(tablep1,digits=4)

xtable(tablep,digits=4)


### plot cumsum
plot(cumsum(profit_bench/100), col="black", type="l", pch =0, ylim=c(0, 0.28), 
     lwd=3, xlab=" ", ylab="Value", main="Cumulative Returns",xaxt = "n")
lines(cumsum(profit_fgmfpca/100), col="red", pch =1,lwd=3)
lines(cumsum(profit_mfpca/100), col="blue", pch =2,lwd=3)
legend("topleft", legend=c("Benchmark", "FG-MFPCA", "X-MFPCA"),
       col=c("black", "red", "blue"), lwd=2, cex=1.8, bty="n")
manual_labels <- c("May-2017", "Jan-2018", "Sep-2018", "May-2019", "Jan-2020", "Sep-2020")  # Your custom labels
axis(1, at = round(seq(1,674,134.5)), labels = manual_labels, las = 1, cex.axis = 1.5)


###################################
## uk
###################################
country = "uk"
true_re = uk_true
setwd('your working directory')

load(file.path(paste0(country,"_forecast_intraday_g.RData"))) # g or x
profit_bench<-profit_tfpca<-profit_dfpca<-profit_lfpca<-profit_mfpca<-matrix(nrow=674, ncol=1, byrow=TRUE) #Set up variables to hold daily trading strategy profits


# hit=c(rep(NA,721))
for(day in 1:674){
  # benchmark - functional mean
  roll_fmean = rowMeans(true_re[,day:(day+win-1)]) # historical functional mean
  
  iBuy <- which.min(roll_fmean[2:287])

  iSell <- iBuy-1+which.max(roll_fmean[iBuy:length(roll_fmean)])

  if((roll_fmean[iSell]-roll_fmean[iBuy])>=mean(roll_fmean) & iBuy!=288){# 287
    Ret0 <- true_re[iSell,day+win] - true_re[iBuy,day+win]
    fee <- 0.0003/100
    Ret <- Ret0 - 2*fee - fee*Ret0
  }else{
    Ret <-0
  }

  # MFPCA - VaR
  if(any(true_re[2:287,day+win]<=g_var_mfpca_uk01[2:287,day])==TRUE){
    iSell_mfpca = which(true_re[(iBuy):287,day+win]<g_var_mfpca_uk01[(iBuy):287,day])[1]
    if(is.na(iSell_mfpca)==TRUE){iSell_mfpca = iSell}
    Ret0_mfpca <- true_re[iSell_mfpca,day+win] - true_re[iBuy,day+win]
    Ret_mfpca <- Ret0_mfpca - 2*fee - fee*Ret0_mfpca
    if(Ret_mfpca<0){Ret_mfpca = 0}
  }else{
    Ret_mfpca <-Ret
  } 

  profit_bench[day]<-Ret
  Ret_mfpca[is.na(Ret_mfpca)] <- 0
  profit_mfpca[day]<-Ret_mfpca
}

tablep1 = matrix(NA,2,3)
tablep1[1,1] = portp(profit_bench)[[1]]
tablep1[2,1] = portp(profit_mfpca)[[1]]

tablep1[1,2] = portp(profit_bench)[[2]]
tablep1[2,2] = portp(profit_mfpca)[[2]]

tablep1[1,3] = portp(profit_bench)[[3]]
tablep1[2,3] = portp(profit_mfpca)[[3]]

profit_fgmfpca=profit_mfpca

#### fgarch-x

load(file.path(paste0(country,"_forecast_intraday_x.RData"))) # g or x
profit_bench<-profit_tfpca<-profit_dfpca<-profit_lfpca<-profit_mfpca<-matrix(nrow=674, ncol=1, byrow=TRUE) #Set up variables to hold daily trading strategy profits

for(day in 1:674){
  # benchmark - functional mean
  roll_fmean = rowMeans(true_re[,day:(day+win-1)]) # historical functional mean
  iBuy <- which.min(roll_fmean[2:287])
  iSell <- iBuy-1+which.max(roll_fmean[iBuy:length(roll_fmean)])

  if((roll_fmean[iSell]-roll_fmean[iBuy])>=mean(roll_fmean) & iBuy!=288){# 287
    Ret0 <- true_re[iSell,day+win] - true_re[iBuy,day+win]
    fee <- 0.0003/100
    Ret <- Ret0 - 2*fee - fee*Ret0
  }else{
    Ret <-0
  }

  # MFPCA - VaR
  if(any(true_re[2:287,day+win]<=x_var_mfpca_uk01[2:287,day])==TRUE){
    iSell_mfpca = which(true_re[(iBuy+1):287,day+win]<x_var_mfpca_uk01[(iBuy+1):287,day])[1]
    if(is.na(iSell_mfpca)==TRUE){iSell_mfpca = iSell}
    Ret0_mfpca <- true_re[iSell_mfpca,day+win] - true_re[iBuy,day+win]
    Ret_mfpca <- Ret0_mfpca - 2*fee - fee*Ret0_mfpca
    if(Ret_mfpca<0){Ret_mfpca = 0}
  }else{
    Ret_mfpca <-Ret
  }
   

  profit_bench[day]<-Ret
  Ret_mfpca[is.na(Ret_mfpca)] <- 0
  profit_mfpca[day]<-Ret_mfpca
}

tablep = matrix(NA,2,3)
tablep[1,1] = portp(profit_bench)[[1]]
tablep[2,1] = portp(profit_mfpca)[[1]]

tablep[1,2] = portp(profit_bench)[[2]]
tablep[2,2] = portp(profit_mfpca)[[2]]

tablep[1,3] = portp(profit_bench)[[3]]
tablep[2,3] = portp(profit_mfpca)[[3]]

xtable(tablep1,digits=4)

xtable(tablep,digits=4)


### plot cumsum
plot(cumsum(profit_bench/100), col="black", type="l", pch =0, ylim=c(-0.15, 0.35), 
     lwd=3, xlab=" ", ylab="Value", main="Cumulative Returns",xaxt = "n")
lines(cumsum(profit_fgmfpca/100), col="red", pch =1,lwd=3)
lines(cumsum(profit_mfpca/100), col="blue", pch =2,lwd=3)
legend("topleft", legend=c("Benchmark", "FG-MFPCA", "X-MFPCA"),
       col=c("black", "red", "blue"), lwd=2, cex=1.8, bty="n")
manual_labels <- c("May-2017", "Jan-2018", "Sep-2018", "May-2019", "Jan-2020", "Sep-2020")  # Your custom labels
axis(1, at = round(seq(1,674,134.5)), labels = manual_labels, las = 1, cex.axis = 1.5)
abline(h = 0, col = "forestgreen", lty = 2, lwd = 2)  

########################
#### jp
########################
country = "jp"
true_re = jp_true
setwd('your working directory')

load(file.path(paste0(country,"_forecast_intraday_g.RData"))) # g or x
profit_bench<-profit_tfpca<-profit_dfpca<-profit_lfpca<-profit_mfpca<-matrix(nrow=674, ncol=1, byrow=TRUE) #Set up variables to hold daily trading strategy profits


# hit=c(rep(NA,721))
for(day in 1:674){
  # benchmark - functional mean
  roll_fmean = rowMeans(true_re[,day:(day+win-1)]) # historical functional mean
  
  iBuy <- which.min(roll_fmean[2:287])

  iSell <- iBuy-1+which.max(roll_fmean[iBuy:length(roll_fmean)])

  if((roll_fmean[iSell]-roll_fmean[iBuy])>=mean(roll_fmean) & iBuy!=288){# 287
    Ret0 <- true_re[iSell,day+win] - true_re[iBuy,day+win]
    fee <- 0.0003/100
    Ret <- Ret0 - 2*fee - fee*Ret0
  }else{
    Ret <-0
  }

  # MFPCA - VaR
  if(any(true_re[2:287,day+win]<=g_var_mfpca_jp01[2:287,day])==TRUE){
    iSell_mfpca = which(true_re[(iBuy):287,day+win]<g_var_mfpca_jp01[(iBuy):287,day])[1]
    if(is.na(iSell_mfpca)==TRUE){iSell_mfpca = iSell}
    Ret0_mfpca <- true_re[iSell_mfpca,day+win] - true_re[iBuy,day+win]
    Ret_mfpca <- Ret0_mfpca - 2*fee - fee*Ret0_mfpca
    if(Ret_mfpca<0){Ret_mfpca = 0}
  }else{
    Ret_mfpca <-Ret
  }
  
  # hit[day]=iSell_mfpca

  profit_bench[day]<-Ret
  Ret_mfpca[is.na(Ret_mfpca)] <- 0
  profit_mfpca[day]<-Ret_mfpca
}

tablep1 = matrix(NA,2,3)
tablep1[1,1] = portp(profit_bench)[[1]]
tablep1[2,1] = portp(profit_mfpca)[[1]]

tablep1[1,2] = portp(profit_bench)[[2]]
tablep1[2,2] = portp(profit_mfpca)[[2]]

tablep1[1,3] = portp(profit_bench)[[3]]
tablep1[2,3] = portp(profit_mfpca)[[3]]

profit_fgmfpca=profit_mfpca

#### fgarch-x

load(file.path(paste0(country,"_forecast_intraday_x.RData"))) # g or x
profit_bench<-profit_tfpca<-profit_dfpca<-profit_lfpca<-profit_mfpca<-matrix(nrow=674, ncol=1, byrow=TRUE) #Set up variables to hold daily trading strategy profits

for(day in 1:674){
  # benchmark - functional mean
  roll_fmean = rowMeans(true_re[,day:(day+win-1)]) # historical functional mean
  iBuy <- which.min(roll_fmean[2:287])
  iSell <- iBuy-1+which.max(roll_fmean[iBuy:length(roll_fmean)])

  if((roll_fmean[iSell]-roll_fmean[iBuy])>=mean(roll_fmean) & iBuy!=288){# 287
    Ret0 <- true_re[iSell,day+win] - true_re[iBuy,day+win]
    fee <- 0.0003/100
    Ret <- Ret0 - 2*fee - fee*Ret0
  }else{
    Ret <-0
  }

  # MFPCA - VaR
  if(any(true_re[2:287,day+win]<=x_var_mfpca_jp01[2:287,day])==TRUE){
    iSell_mfpca = which(true_re[(iBuy+1):287,day+win]<x_var_mfpca_jp01[(iBuy+1):287,day])[1]
    if(is.na(iSell_mfpca)==TRUE){iSell_mfpca = iSell}
    Ret0_mfpca <- true_re[iSell_mfpca,day+win] - true_re[iBuy,day+win]
    Ret_mfpca <- Ret0_mfpca - 2*fee - fee*Ret0_mfpca
    if(Ret_mfpca<0){Ret_mfpca = 0}
  }else{
    Ret_mfpca <-Ret
  }
  
  # hit[day]=iSell_mfpca

  profit_bench[day]<-Ret
  Ret_mfpca[is.na(Ret_mfpca)] <- 0
  profit_mfpca[day]<-Ret_mfpca
}

tablep = matrix(NA,2,3)
tablep[1,1] = portp(profit_bench)[[1]]
tablep[2,1] = portp(profit_mfpca)[[1]]

tablep[1,2] = portp(profit_bench)[[2]]
tablep[2,2] = portp(profit_mfpca)[[2]]

tablep[1,3] = portp(profit_bench)[[3]]
tablep[2,3] = portp(profit_mfpca)[[3]]

# get the result for Table 6.2
xtable(tablep1,digits=4)

xtable(tablep,digits=4)

