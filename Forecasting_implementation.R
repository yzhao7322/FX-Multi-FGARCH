library('sde')
library('xtable')
library('mvtnorm')
library('tseries')
library('optimx')
library('compiler')
library('expm')
library('quantreg')
library('fGarch')
library('rugarch')
library('DescTools')
library('Jmisc')
library('matrixStats')
library('POT')
library('xts')
library('sandwich')
library('reshape2')
library('zoo')
library("nsprcomp")
library("QZ")
library('ftsa')
library('LongMemoryTS')

#############
### load data
#############
setwd('/your working directory/')

euspot = read.csv('euspot_24.csv',sep=",",header=FALSE)
ukspot = read.csv('ukspot_24.csv',sep=",",header=FALSE)
jpspot = read.csv('jpspot_24.csv',sep=",",header=FALSE)
eubid = read.csv('eubid_24.csv',sep=",",header=FALSE)
euask = read.csv('euask_24.csv',sep=",",header=FALSE)
ukbid = read.csv('ukbid_24.csv',sep=",",header=FALSE)
ukask = read.csv('ukask_24.csv',sep=",",header=FALSE)
jpbid = read.csv('jpbid_24.csv',sep=",",header=FALSE)
jpask = read.csv('jpask_24.csv',sep=",",header=FALSE)


###############
### construct return and bid-ask spread
###############
intra.return <- function(yd){
  grid_point=nrow(yd)
  N=ncol(yd)
  idr_re=matrix(0,grid_point-1,N)
  cidr_re=matrix(0,grid_point,N)
  ocidr_re=matrix(0,grid_point,N-1)
  for(j in c(2:grid_point)){
    for(i in c(1:N)){
      idr_re[j-1,i]=log(yd[j,i])-log(yd[j-1,i])
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
  return(list(idr = idr_re, cidr = cidr_re, ocidr = (ocidr_re * 100) ))
}

bid.ask.spread <- function(bid, ask){
   grid_point=nrow(bid) 
   N = ncol(bid)
   spread = spr_re = matrix(NA,nrow(bid),ncol(bid))
   
   for (i in c(1:N)){
     ask_d = ask[,i]
     bid_d = bid[,i]
     spread[,i] = (ask_d-bid_d) #/mid_d
   }
   
    for(i in c(2:N)){
      ae = (spread[,i] - spread[grid_point,i-1])
      # ae = Winsorize(ae, probs = c(0.05, 0.95) )
      spr_re[,i]= ae
    }
    spr_re = spr_re[,2:N]
   return(spr_re)
}


eu_re = intra.return(euspot)$ocidr
uk_re = intra.return(ukspot)$ocidr
jp_re = intra.return(jpspot)$ocidr

eu_spr =bid.ask.spread(eubid,euask)
uk_spr =bid.ask.spread(ukbid,ukask)
jp_spr =bid.ask.spread(jpbid,jpask)

# triming data
for(i in 1:ncol(eu_re) ){
   eu_re[,i] = Winsorize(eu_re[,i], probs = c(0.02, 0.98) )   
}
for(i in 1:ncol(uk_re) ){
   uk_re[,i] = Winsorize(uk_re[,i], probs = c(0.02, 0.98) )   
}
for(i in 1:ncol(jp_re) ){
   jp_re[,i] = Winsorize(jp_re[,i], probs = c(0.02, 0.98) )    
}

for(i in 1:ncol(eu_spr) ){
   eu_spr[,i] = Winsorize(eu_spr[,i], probs = c(0.02, 0.98) )   
}
for(i in 1:ncol(uk_spr) ){
   uk_spr[,i] = Winsorize(uk_spr[,i], probs = c(0.02, 0.98) )   
}
for(i in 1:ncol(jp_spr) ){
   jp_spr[,i] = Winsorize(jp_spr[,i], probs = c(0.02, 0.98) )    
}

# demean data
eu_re=eu_re-rowMeans(eu_re)
uk_re=uk_re-rowMeans(uk_re)
jp_re=jp_re-rowMeans(jp_re)
eu_spr=eu_spr-rowMeans(eu_spr)
uk_spr=uk_spr-rowMeans(uk_spr)
jp_spr=jp_spr-rowMeans(jp_spr)

eu_x = abs(eu_spr)
uk_x = abs(uk_spr)
jp_x = abs(jp_spr)

# here we take USD/EUR as an example.
# USD/GBP and USD/JPY can be performed in a similar manner.
# However, for USD/GBP and USD/JPY, the function forecast_multi_X<-function(cidr_data, cidr_data2, cidr_data3, x_data)
# needs to be carefully checked, as the objective bases should be speficied for the corresponding currency pairs.
# especially for the FGARCH-X model, the covariate X should be replaced correspondingly.

# FGARCH - TFPCA 
pred_eu = forecast_FGARCH(eu_re, "tfpca")
save(pred_eu, file = "fgarch_tfpca_eu.RData")

# FGARCH - DFPCA
pred_uk = forecast_FGARCH(uk_re, "dfpca")
save(pred_uk, file = "fgarch_dpca_uk.RData")

# FGARCH - LFPCA
pred_eu = forecast_FGARCH(eu_re, "lfpca")
save(pred_eu, file = "fgarch_lpca_eu.RData")

# FGARCH - MFPCA
pred_eu = forecast_multi_FGARCH(eu_re,uk_re,jp_re)
save(pred_eu, file = "fgarch_mfpca_eu.RData")



# FGARCHX - TFPCA 
pred_eu = forecast_X(eu_re,eu_x, "tfpca")
save(pred_eu, file = "x_tfpca_eu.RData")

# FGARCHX - DFPCA
pred_eu = forecast_X(eu_re,eu_x, "dfpca")
save(pred_eu, file = "x_dpca_eu.RData")

# FGARCHX - LFPCA
pred_eu = forecast_X(eu_re,eu_x, "lfpca")
save(pred_eu, file = "x_lpca_eu.RData")

# FGARCHX - MFPCA
pred_eu = forecast_multi_X(eu_re, uk_re, jp_re, eu_x)
save(pred_eu, file = "x_mfpca_eu.RData")


###########################################################################
### evaluation

eu_res = intra.return(euspot)
uk_res = intra.return(ukspot)
jp_res = intra.return(jpspot)

eu_ocidr = eu_res$ocidr
uk_ocidr = uk_res$ocidr
jp_ocidr = jp_res$ocidr

full_N = ncol(eu_ocidr)
grid_point = nrow(eu_ocidr)

####################################
#### intraday volatility evaluation
####################################

eu_rv = (eu_ocidr[,601:1274]^2)
uk_rv = (uk_ocidr[,601:1274]^2)
jp_rv = (jp_ocidr[,601:1274]^2)

#### conditional volatility forecast

load("fgarch_tfpca_eu.RData")
g_tfpca_eu = pred_eu[[2]][,1:674]
load("fgarch_tfpca_uk.RData")
g_tfpca_uk = pred_uk[[2]][,1:674]
load("fgarch_tfpca_jp.RData")
g_tfpca_jp = pred_jp[[2]][,1:674]

load("fgarch_dpca_eu.RData")
g_dpca_eu = pred_eu[[2]][,1:674]
load("fgarch_dpca_uk.RData")
g_dpca_uk = pred_uk[[2]][,1:674]
load("fgarch_dpca_jp.RData")
g_dpca_jp = pred_jp[[2]][,1:674]

load("fgarch_lpca_eu.RData")
g_lpca_eu = pred_eu[[2]][,1:674]
load("fgarch_lpca_uk.RData")
g_lpca_uk = pred_uk[[2]][,1:674]
load("fgarch_lpca_jp.RData")
g_lpca_jp = pred_jp[[2]][,1:674]

load("fgarch_mfpca_eu.RData")
g_mfpca_eu = pred_eu[[2]][,1:674]
load("fgarch_mfpca_uk.RData")
g_mfpca_uk = pred_uk[[2]][,1:674]
load("fgarch_mfpca_jp.RData")
g_mfpca_jp = pred_jp[[2]][,1:674]


load("x_tfpca_eu.RData")
x_tfpca_eu = pred_eu[[2]][,1:674]
load("x_tfpca_uk.RData")
x_tfpca_uk = pred_uk[[2]][,1:674]
load("x_tfpca_jp.RData")
x_tfpca_jp = pred_jp[[2]][,1:674]

load("x_dpca_eu.RData")
x_dpca_eu = pred_eu[[2]][,1:674]
load("x_dpca_uk.RData")
x_dpca_uk = pred_uk[[2]][,1:674]
load("x_dpca_jp.RData")
x_dpca_jp = pred_jp[[2]][,1:674]

load("x_lpca_eu.RData")
x_lpca_eu = pred_eu[[2]][,1:674]
load("x_lpca_uk.RData")
x_lpca_uk = pred_uk[[2]][,1:674]
load("x_lpca_jp.RData")
x_lpca_jp = pred_jp[[2]][,1:674]

load("x_mfpca_eu.RData")
x_mfpca_eu = pred_eu[[2]][,1:674]
load("x_mfpca_uk.RData")
x_mfpca_uk = pred_uk[[2]][,1:674]
load("x_mfpca_jp.RData")
x_mfpca_jp = pred_jp[[2]][,1:674]


### forecasting errors
### mse/qlike
forecast_error<-function(true, pred){
  out_sampleT = ncol(true)
  grid_point = nrow(true)
  maein = colSums(abs(true-pred))/grid_point
  msein = colSums((true-pred)^2)/grid_point
  qlikein = colSums( log(pred) + true/pred )/grid_point
  mae = sum( maein )/out_sampleT
  mse = sum( msein )/out_sampleT
  qlike = mean( qlikein, na.rm = TRUE)
  return(list(mae=mae, mse=mse, qlike=qlike))
}

eu_rv1=eu_rv
uk_rv1=uk_rv
jp_rv1=jp_rv

fore_err_table = matrix(NA,8,6)
fore_err = forecast_error(eu_rv1,g_tfpca_eu)
fore_err_table[1,1] = fore_err$mse
fore_err_table[1,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_tfpca_uk)
fore_err_table[1,3] = fore_err$mse
fore_err_table[1,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_tfpca_jp)
fore_err_table[1,5] = fore_err$mse
fore_err_table[1,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,g_dpca_eu)
fore_err_table[2,1] = fore_err$mse
fore_err_table[2,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_dpca_uk)
fore_err_table[2,3] = fore_err$mse
fore_err_table[2,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_dpca_jp)
fore_err_table[2,5] = fore_err$mse
fore_err_table[2,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,g_lpca_eu)
fore_err_table[3,1] = fore_err$mse
fore_err_table[3,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_lpca_uk)
fore_err_table[3,3] = fore_err$mse
fore_err_table[3,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_lpca_jp)
fore_err_table[3,5] = fore_err$mse
fore_err_table[3,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,g_mfpca_eu)
fore_err_table[4,1] = fore_err$mse
fore_err_table[4,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_mfpca_uk)
fore_err_table[4,3] = fore_err$mse
fore_err_table[4,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_mfpca_jp)
fore_err_table[4,5] = fore_err$mse
fore_err_table[4,6] = fore_err$qlike



fore_err = forecast_error(eu_rv1,x_tfpca_eu)
fore_err_table[5,1] = fore_err$mse
fore_err_table[5,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_tfpca_uk)
fore_err_table[5,3] = fore_err$mse
fore_err_table[5,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_tfpca_jp)
fore_err_table[5,5] = fore_err$mse
fore_err_table[5,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,x_dpca_eu)
fore_err_table[6,1] = fore_err$mse
fore_err_table[6,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_dpca_uk)
fore_err_table[6,3] = fore_err$mse
fore_err_table[6,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_dpca_jp)
fore_err_table[6,5] = fore_err$mse
fore_err_table[6,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,x_lpca_eu)
fore_err_table[7,1] = fore_err$mse
fore_err_table[7,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_lpca_uk)
fore_err_table[7,3] = fore_err$mse
fore_err_table[7,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_lpca_jp)
fore_err_table[7,5] = fore_err$mse
fore_err_table[7,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,x_mfpca_eu)
fore_err_table[8,1] = fore_err$mse
fore_err_table[8,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_mfpca_uk)
fore_err_table[8,3] = fore_err$mse
fore_err_table[8,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_mfpca_jp)
fore_err_table[8,5] = fore_err$mse
fore_err_table[8,6] = fore_err$qlike

xtable(fore_err_table,digits=4)



###################
#### comparasion methods
# again, take USD/EUR as an example, USD/GBP, USD/JPY are the same
###################
eu_rv1=eu_rv
uk_rv1=uk_rv
jp_rv1=jp_rv

loss_f <- function(true, pred){
  grid_point = nrow(true)
  mse = colSums((pred-true)^2)/grid_point
  qlike = colSums( log(pred) + true/pred  )/grid_point

  return(list(mse=mse, qlike=qlike))
}

eu_loss = loss_f(eu_rv1,g_tfpca_eu)
g_tfpca_eu_mse = eu_loss$mse
g_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_dpca_eu)
g_dpca_eu_mse = eu_loss$mse
g_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_lpca_eu)
g_lpca_eu_mse = eu_loss$mse
g_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_mfpca_eu)
g_mfpca_eu_mse = eu_loss$mse
g_mfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_tfpca_eu)
x_tfpca_eu_mse = eu_loss$mse
x_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_dpca_eu)
x_dpca_eu_mse = eu_loss$mse
x_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_lpca_eu)
x_lpca_eu_mse = eu_loss$mse
x_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_mfpca_eu)
x_mfpca_eu_mse = eu_loss$mse
x_mfpca_eu_qlike = eu_loss$qlike

loss_eu_mse = cbind(g_tfpca_eu_mse,g_dpca_eu_mse,g_lpca_eu_mse,g_mfpca_eu_mse,x_tfpca_eu_mse,x_dpca_eu_mse,x_lpca_eu_mse,x_mfpca_eu_mse)
loss_eu_qlike = cbind(g_tfpca_eu_qlike,g_dpca_eu_qlike,g_lpca_eu_qlike,g_mfpca_eu_qlike,x_tfpca_eu_qlike,x_dpca_eu_qlike,x_lpca_eu_qlike,x_mfpca_eu_qlike)

MCSprocedure(loss_eu_mse, alpha = 0.95, B = 5000, cl = NULL,
ram.allocation = TRUE, statistic = "Tmax", k = NULL, min.k = 3,
verbose = TRUE)
MCSprocedure(as.matrix(na.omit(loss_eu_qlike)), alpha = 0.95, B = 5000, cl = NULL,
ram.allocation = TRUE, statistic = "Tmax", k = NULL, min.k = 3,
verbose = TRUE)



loss_f <- function(true, pred){
  # true = Winsorize(true, quantile(true, probs = c(0.01, 0.99), na.rm = FALSE))
  grid_point = nrow(true)
  mse = colSums((pred-true)^2)/grid_point
  qlike = colSums( log(pred) + true/pred  )/grid_point
  return(list(mse=mse, qlike=qlike))
}

eu_loss = loss_f(eu_rv1,g_tfpca_eu)
g_tfpca_eu_mse = eu_loss$mse
g_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_dpca_eu)
g_dpca_eu_mse = eu_loss$mse
g_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_lpca_eu)
g_lpca_eu_mse = eu_loss$mse
g_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_mfpca_eu)
g_mfpca_eu_mse = eu_loss$mse
g_mfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_tfpca_eu)
x_tfpca_eu_mse = eu_loss$mse
x_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_dpca_eu)
x_dpca_eu_mse = eu_loss$mse
x_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_lpca_eu)
x_lpca_eu_mse = eu_loss$mse
x_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_mfpca_eu)
x_mfpca_eu_mse = eu_loss$mse
x_mfpca_eu_qlike = eu_loss$qlike

dm_tab=matrix(NA,16,3)
hd = 1
dm_tab[1,1]=dm.test(g_tfpca_eu_mse,g_dpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[2,1]=dm.test(g_tfpca_eu_mse,g_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[3,1]=dm.test(g_tfpca_eu_mse,g_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[4,1]=dm.test(g_dpca_eu_mse,g_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[5,1]=dm.test(g_dpca_eu_mse,g_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[6,1]=dm.test(g_lpca_eu_mse,g_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[7,1]=dm.test(x_tfpca_eu_mse,x_dpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[8,1]=dm.test(x_tfpca_eu_mse,x_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[9,1]=dm.test(x_tfpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[10,1]=dm.test(x_dpca_eu_mse,x_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[11,1]=dm.test(x_dpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[12,1]=dm.test(x_lpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[13,1]=dm.test(g_tfpca_eu_mse,x_tfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[14,1]=dm.test(g_dpca_eu_mse,x_dpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[15,1]=dm.test(g_lpca_eu_mse,x_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]
dm_tab[16,1]=dm.test(g_mfpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = hd, power = 2)[[1]]



load("fgarch_tfpca_eu.RData")
g_tfpca_eu = pred_eu[[2]][,1:674]
load("fgarch_lpca_eu.RData")
g_lpca_eu = pred_eu[[2]][,1:674]
load("fgarch_dpca_eu.RData")
g_dpca_eu = pred_eu[[2]][,1:674]
load("fgarch_mfpca_eu.RData")
g_mfpca_eu = pred_eu[[2]][,1:674]
# load("fgarch_lfpca_eu.RData")
# g_lfpca_eu = pred_eu[[2]][,1:674]
load("x_mfpca_eu.RData")
x_mfpca_eu = pred_eu[[2]][,1:674]

##################################################### Figure 5.1
library(rgl)

plot_tar = x_mfpca_eu # g_mfpca_eu x_mfpca_eu
x <- seq(1, nrow(plot_tar), length.out = nrow(plot_tar))
y <- seq(1, ncol(plot_tar), length.out = ncol(plot_tar))
z <- plot_tar
# Generate color gradient using terrain.colors, heat.colors, or rainbow
# col_matrix <- terrain.colors(100)[cut(z, breaks = 100)] 
col_matrix <- rainbow(100)[cut(z, breaks = 100, labels = FALSE)]  # Corrected color mapping
open3d()
persp3d(x, y, z, col = col_matrix, alpha = 0.9, 
        xlab = "Intraday Interval", ylab = "", zlab = "Value",
        axes = FALSE, main = " ")  # Keep default axes (except y)
manual_labels <- c("May-2017", "Jan-2018", "Sep-2018", "May-2019", "Jan-2020", "Sep-2020")
y_positions <- round(seq(1, ncol(plot_tar), length.out = length(manual_labels)))  # Positioning labels
axis3d("x", col = "black", lwd = 3)  
axis3d("z", col = "deeppink1", lwd = 3, cex = 2)  
axis3d("y", at = y_positions, labels = manual_labels, col = "black", cex.axis = 1.5, lwd = 3)



##################################################### Figure 5.2

load("fgarch_tfpca_eu.RData")
g_tfpca_eu = pred_eu[[2]][,1:674] 
load("fgarch_dpca_eu.RData")
g_dpca_eu = pred_eu[[2]][,1:674] 
load("fgarch_lpca_eu.RData")
g_lpca_eu = pred_eu[[2]][,1:674] 
load("fgarch_mfpca_eu.RData")
g_mfpca_eu = pred_eu[[2]][,1:674] 


load("x_tfpca_eu.RData")
x_tfpca_eu = pred_eu[[2]][,1:674] 
load("x_dpca_eu.RData")
x_dpca_eu = pred_eu[[2]][,1:674] 
load("x_lpca_eu.RData")
x_lpca_eu = pred_eu[[2]][,1:674] 
load("x_mfpca_eu.RData")
x_mfpca_eu = pred_eu[[2]][,1:674] 


eu_rv_mean = rowMeans(eu_rv)
g_tfpca_eu_mean = rowMeans(g_tfpca_eu)
g_dpca_eu_mean = rowMeans(g_dpca_eu)
g_lpca_eu_mean = rowMeans(g_lpca_eu)
g_mfpca_eu_mean = rowMeans(g_mfpca_eu)

x_tfpca_eu_mean = rowMeans(x_tfpca_eu)
x_dpca_eu_mean = rowMeans(x_dpca_eu)
x_lpca_eu_mean = rowMeans(x_lpca_eu)
x_mfpca_eu_mean = rowMeans(x_mfpca_eu)

library(matrixStats)
eu_rv_median = rowMedians(eu_rv)
g_tfpca_eu_median = rowMedians(g_tfpca_eu)
g_dpca_eu_median = rowMedians(g_dpca_eu)
g_lpca_eu_median = rowMedians(g_lpca_eu)
g_mfpca_eu_median = rowMedians(g_mfpca_eu)

x_tfpca_eu_median = rowMedians(x_tfpca_eu)
x_dpca_eu_median = rowMedians(x_dpca_eu)
x_lpca_eu_median = rowMedians(x_lpca_eu)
x_mfpca_eu_median = rowMedians(x_mfpca_eu)


plot(eu_rv_mean, lty=1,type="o", col="black", ylim = c(0, 1.5), lwd = 2, pch = 16 )
lines(g_tfpca_eu_mean, lty=2, type="o", col="chartreuse1", lwd = 2, pch = 17 ) 
lines(g_dpca_eu_mean, lty=2, type="o", col="azure4", lwd = 2, pch = 17 )  
lines(g_lpca_eu_mean, lty=2, type="o", col="cyan1", lwd = 2, pch = 17 ) 
lines(g_mfpca_eu_mean, lty=2, type="o", col="red", lwd = 2, pch = 17 )  
lines(x_tfpca_eu_mean, lty=3, type="o", col="gold1", lwd = 2, pch = 3 ) 
lines(x_dpca_eu_mean, lty=3, type="o", col="seagreen", lwd = 2, pch = 3 )  
lines(x_lpca_eu_mean, lty=3, type="o", col="dodgerblue", lwd = 2, pch = 3 )  
lines(x_mfpca_eu_mean, lty=3, type="o", col="deeppink", lwd = 2, pch = 3 )


####################################
#### inter-daily volatility evaluation
####################################

eu_ocidr = intra.return(euspot)$ocidr
uk_ocidr = intra.return(ukspot)$ocidr
jp_ocidr = intra.return(jpspot)$ocidr

eu_idr = intra.return(euspot)$idr[,2:1275]
uk_idr = intra.return(ukspot)$idr[,2:1275]
jp_idr = intra.return(jpspot)$idr[,2:1275]

full_N = ncol(eu_ocidr)
grid_point = nrow(eu_ocidr)

ip = 193 # 193 #108

realvol_overnight<-function(data,scalar_onr){
  point_grid=nrow(data)
  N=ncol(data)
  real_onvol=matrix(0,N,1)
  for (i in 1:N){
    real_onvol[i,]=sum(data[,i]^2)+scalar_onr[i,]^2
  }
  return(real_onvol)
}

eu_realv = realvol_overnight(eu_idr[1:ip,],as.matrix(eu_ocidr[1,],full_N,1))
uk_realv = realvol_overnight(uk_idr[1:ip,],as.matrix(uk_ocidr[1,],full_N,1))
jp_realv = realvol_overnight(jp_idr[1:ip,],as.matrix(jp_ocidr[1,],full_N,1))

eu_rv = eu_realv[601:1274]
uk_rv = uk_realv[601:1274]
jp_rv = jp_realv[601:1274]

#############################################  load prediction results
load("fgarch_tfpca_eu.RData")
g_tfpca_eu = pred_eu[[2]][ip,1:674]
load("fgarch_tfpca_uk.RData")
g_tfpca_uk = pred_uk[[2]][ip,1:674]
load("fgarch_tfpca_jp.RData")
g_tfpca_jp = pred_jp[[2]][ip,1:674]

load("fgarch_dpca_eu.RData")
g_dpca_eu = pred_eu[[2]][ip,1:674]
load("fgarch_dpca_uk.RData")
g_dpca_uk = pred_uk[[2]][ip,1:674]
load("fgarch_dpca_jp.RData")
g_dpca_jp = pred_jp[[2]][ip,1:674]

load("fgarch_lpca_eu.RData")
g_lpca_eu = pred_eu[[2]][ip,1:674]
load("fgarch_lpca_uk.RData")
g_lpca_uk = pred_uk[[2]][ip,1:674]
load("fgarch_lpca_jp.RData")
g_lpca_jp = pred_jp[[2]][ip,1:674]

load("fgarch_mfpca_eu.RData")
g_mfpca_eu = pred_eu[[2]][ip,1:674]
load("fgarch_mfpca_uk.RData")
g_mfpca_uk = pred_uk[[2]][ip,1:674]
load("fgarch_mfpca_jp.RData")
g_mfpca_jp = pred_jp[[2]][ip,1:674]

load("x_tfpca_eu.RData")
x_tfpca_eu = pred_eu[[2]][ip,1:674]
load("x_tfpca_uk.RData")
x_tfpca_uk = pred_uk[[2]][ip,1:674]
load("x_tfpca_jp.RData")
x_tfpca_jp = pred_jp[[2]][ip,1:674]

load("x_dpca_eu.RData")
x_dpca_eu = pred_eu[[2]][ip,1:674]
load("x_dpca_uk.RData")
x_dpca_uk = pred_uk[[2]][ip,1:674]
load("x_dpca_jp.RData")
x_dpca_jp = pred_jp[[2]][ip,1:674]

load("x_lpca_eu.RData")
x_lpca_eu = pred_eu[[2]][ip,1:674]
load("x_lpca_uk.RData")
x_lpca_uk = pred_uk[[2]][ip,1:674]
load("x_lpca_jp.RData")
x_lpca_jp = pred_jp[[2]][ip,1:674]

load("x_mfpca_eu.RData")
x_mfpca_eu = pred_eu[[2]][ip,1:674]
load("x_mfpca_uk.RData")
x_mfpca_uk = pred_uk[[2]][ip,1:674]
load("x_mfpca_jp.RData")
x_mfpca_jp = pred_jp[[2]][ip,1:674]

### forecasting errors
### mse/qlike
forecast_error<-function(true, pred){
  out_sampleT = length(true)
  maein = abs(true-pred)
  msein = (true-pred)^2
  qlikein = log(pred) + true/pred
 
  mae = sum(maein)/out_sampleT
  mse = sum(msein)/out_sampleT
  qlike = sum(qlikein)/out_sampleT

  return(list(mae=mae, mse=mse, qlike=qlike))
}

eu_rv1=eu_rv
uk_rv1=uk_rv
jp_rv1=jp_rv

fore_err_table = matrix(NA,8,6)
fore_err = forecast_error(eu_rv1,g_tfpca_eu)
fore_err_table[1,1] = fore_err$mse
fore_err_table[1,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_tfpca_uk)
fore_err_table[1,3] = fore_err$mse
fore_err_table[1,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_tfpca_jp)
fore_err_table[1,5] = fore_err$mse
fore_err_table[1,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,g_dpca_eu)
fore_err_table[2,1] = fore_err$mse
fore_err_table[2,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_dpca_uk)
fore_err_table[2,3] = fore_err$mse
fore_err_table[2,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_dpca_jp)
fore_err_table[2,5] = fore_err$mse
fore_err_table[2,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,g_lpca_eu)
fore_err_table[3,1] = fore_err$mse
fore_err_table[3,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_lpca_uk)
fore_err_table[3,3] = fore_err$mse
fore_err_table[3,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_lpca_jp)
fore_err_table[3,5] = fore_err$mse
fore_err_table[3,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,g_mfpca_eu)
fore_err_table[4,1] = fore_err$mse
fore_err_table[4,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,g_mfpca_uk)
fore_err_table[4,3] = fore_err$mse
fore_err_table[4,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,g_mfpca_jp)
fore_err_table[4,5] = fore_err$mse
fore_err_table[4,6] = fore_err$qlike



fore_err = forecast_error(eu_rv1,x_tfpca_eu)
fore_err_table[5,1] = fore_err$mse
fore_err_table[5,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_tfpca_uk)
fore_err_table[5,3] = fore_err$mse
fore_err_table[5,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_tfpca_jp)
fore_err_table[5,5] = fore_err$mse
fore_err_table[5,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,x_dpca_eu)
fore_err_table[6,1] = fore_err$mse
fore_err_table[6,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_dpca_uk)
fore_err_table[6,3] = fore_err$mse
fore_err_table[6,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_dpca_jp)
fore_err_table[6,5] = fore_err$mse
fore_err_table[6,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,x_lpca_eu)
fore_err_table[7,1] = fore_err$mse
fore_err_table[7,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_lpca_uk)
fore_err_table[7,3] = fore_err$mse
fore_err_table[7,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_lpca_jp)
fore_err_table[7,5] = fore_err$mse
fore_err_table[7,6] = fore_err$qlike

fore_err = forecast_error(eu_rv1,x_mfpca_eu)
fore_err_table[8,1] = fore_err$mse
fore_err_table[8,2] = fore_err$qlike

fore_err = forecast_error(uk_rv1,x_mfpca_uk)
fore_err_table[8,3] = fore_err$mse
fore_err_table[8,4] = fore_err$qlike

fore_err = forecast_error(jp_rv1,x_mfpca_jp)
fore_err_table[8,5] = fore_err$mse
fore_err_table[8,6] = fore_err$qlike

xtable(fore_err_table,digits=4)


############################
######### comparasion method
eu_rv1=eu_rv
uk_rv1=uk_rv
jp_rv1=jp_rv

loss_f <- function(true, pred){
  # true = Winsorize(true, quantile(true, probs = c(0.005, 0.995), na.rm = FALSE))
  out_sampleT = length(true)
  mse = (true-pred)^2
  qlike = log(pred) + true/pred

  return(list(mse=mse, qlike=qlike))
}

eu_loss = loss_f(eu_rv1,g_tfpca_eu)
g_tfpca_eu_mse = eu_loss$mse
g_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_dpca_eu)
g_dpca_eu_mse = eu_loss$mse
g_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_lpca_eu)
g_lpca_eu_mse = eu_loss$mse
g_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_mfpca_eu)
g_mfpca_eu_mse = eu_loss$mse
g_mfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_tfpca_eu)
x_tfpca_eu_mse = eu_loss$mse
x_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_dpca_eu)
x_dpca_eu_mse = eu_loss$mse
x_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_lpca_eu)
x_lpca_eu_mse = eu_loss$mse
x_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_mfpca_eu)
x_mfpca_eu_mse = eu_loss$mse
x_mfpca_eu_qlike = eu_loss$qlike
 

loss_eu_mse = cbind(g_tfpca_eu_mse,g_dpca_eu_mse,g_lpca_eu_mse,g_mfpca_eu_mse,x_tfpca_eu_mse,x_dpca_eu_mse,x_lpca_eu_mse,x_mfpca_eu_mse)
loss_eu_qlike = cbind(g_tfpca_eu_qlike,g_dpca_eu_qlike,g_lpca_eu_qlike,g_mfpca_eu_qlike,x_tfpca_eu_qlike,x_dpca_eu_qlike,x_lpca_eu_qlike,x_mfpca_eu_qlike)
 
MCSprocedure(loss_eu_mse, alpha = 0.95, B = 5000, cl = NULL,
ram.allocation = TRUE, statistic = "Tmax", k = NULL, min.k = 3,
verbose = TRUE)
MCSprocedure(loss_eu_qlike, alpha = 0.95, B = 5000, cl = NULL,
ram.allocation = TRUE, statistic = "Tmax", k = NULL, min.k = 3,
verbose = TRUE)

###################
#### DM test
###################

ip=193 

eu_rv1=eu_rv
uk_rv1=uk_rv
jp_rv1=jp_rv

loss_f <- function(true, pred){
  # true = Winsorize(true, quantile(true, probs = c(0.01, 0.99), na.rm = FALSE))
  out_sampleT = length(true)
  mse = (true-pred)^2
  qlike = log(pred) + true/pred

  return(list(mse=mse, qlike=qlike))
}

eu_loss = loss_f(eu_rv1,g_tfpca_eu)
g_tfpca_eu_mse = eu_loss$mse
g_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_dpca_eu)
g_dpca_eu_mse = eu_loss$mse
g_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_lpca_eu)
g_lpca_eu_mse = eu_loss$mse
g_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,g_mfpca_eu)
g_mfpca_eu_mse = eu_loss$mse
g_mfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_tfpca_eu)
x_tfpca_eu_mse = eu_loss$mse
x_tfpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_dpca_eu)
x_dpca_eu_mse = eu_loss$mse
x_dpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_lpca_eu)
x_lpca_eu_mse = eu_loss$mse
x_lpca_eu_qlike = eu_loss$qlike

eu_loss = loss_f(eu_rv1,x_mfpca_eu)
x_mfpca_eu_mse = eu_loss$mse
x_mfpca_eu_qlike = eu_loss$qlike

dm_tab=matrix(NA,16,3)

dm_tab[1,1]=dm.test(g_tfpca_eu_mse,g_dpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[2,1]=dm.test(g_tfpca_eu_mse,g_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[3,1]=dm.test(g_tfpca_eu_mse,g_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[4,1]=dm.test(g_dpca_eu_mse,g_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[5,1]=dm.test(g_dpca_eu_mse,g_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[6,1]=dm.test(g_lpca_eu_mse,g_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[7,1]=dm.test(x_tfpca_eu_mse,x_dpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[8,1]=dm.test(x_tfpca_eu_mse,x_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[9,1]=dm.test(x_tfpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[10,1]=dm.test(x_dpca_eu_mse,x_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[11,1]=dm.test(x_dpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[12,1]=dm.test(x_lpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[13,1]=dm.test(g_tfpca_eu_mse,x_tfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[14,1]=dm.test(g_dpca_eu_mse,x_dpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[15,1]=dm.test(g_lpca_eu_mse,x_lpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[16,1]=dm.test(g_mfpca_eu_mse,x_mfpca_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]



############################################################################
#### classic inter-daily conditional volatility comparasion
############################################################################
ip=193

eu_rv1=eu_rv
uk_rv1=uk_rv
jp_rv1=jp_rv


load("x_mfpca_eu.RData")
x_mfpca_eu = pred_eu[[2]][ip,1:674]
load("x_mfpca_uk.RData")
x_mfpca_uk = pred_uk[[2]][ip,1:674]
load("x_mfpca_jp.RData")
x_mfpca_jp = pred_jp[[2]][ip,1:674]


######## PEER
forecast_error1<-function(true, pred){
  out_sampleT = length(true)
  maein = abs(true-pred)
  msein = (true-pred)^2
  qlikein = log(pred) + true/pred

  mae = sum(maein)/out_sampleT
  mse = sum(msein)/out_sampleT
  qlike = sum(qlikein)/out_sampleT

  return(list(mae=mae, mse=mse, qlike=qlike))
}

peer_table = matrix(NA,9,6)
fore_err = forecast_error1(eu_rv1,eu_garch)
peer_table[1,1] = fore_err$mse
peer_table[1,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_gjr)
peer_table[2,1] = fore_err$mse
peer_table[2,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_x)
peer_table[3,1] = fore_err$mse
peer_table[3,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_fi)
peer_table[4,1] = fore_err$mse
peer_table[4,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_har)
peer_table[5,1] = fore_err$mse
peer_table[5,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_harx)
peer_table[6,1] = fore_err$mse
peer_table[6,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_rvarma)
peer_table[7,1] = fore_err$mse
peer_table[7,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_rgarch)
peer_table[8,1] = fore_err$mse
peer_table[8,2] = fore_err$qlike
fore_err = forecast_error1(eu_rv1,eu_midas)
peer_table[9,1] = fore_err$mse
peer_table[9,2] = fore_err$qlike

fore_err = forecast_error1(uk_rv1,uk_garch)
peer_table[1,3] = fore_err$mse
peer_table[1,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_gjr)
peer_table[2,3] = fore_err$mse
peer_table[2,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_x)
peer_table[3,3] = fore_err$mse
peer_table[3,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_fi)
peer_table[4,3] = fore_err$mse
peer_table[4,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_har)
peer_table[5,3] = fore_err$mse
peer_table[5,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_harx)
peer_table[6,3] = fore_err$mse
peer_table[6,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_rvarma)
peer_table[7,3] = fore_err$mse
peer_table[7,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_rgarch)
peer_table[8,3] = fore_err$mse
peer_table[8,4] = fore_err$qlike
fore_err = forecast_error1(uk_rv1,uk_midas)
peer_table[9,3] = fore_err$mse
peer_table[9,4] = fore_err$qlike

fore_err = forecast_error1(jp_rv1,jp_garch)
peer_table[1,5] = fore_err$mse
peer_table[1,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_gjr)
peer_table[2,5] = fore_err$mse
peer_table[2,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_x)
peer_table[3,5] = fore_err$mse
peer_table[3,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_fi)
peer_table[4,5] = fore_err$mse
peer_table[4,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_har)
peer_table[5,5] = fore_err$mse
peer_table[5,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_harx)
peer_table[6,5] = fore_err$mse
peer_table[6,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_rvarma)
peer_table[7,5] = fore_err$mse
peer_table[7,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_rgarch)
peer_table[8,5] = fore_err$mse
peer_table[8,6] = fore_err$qlike
fore_err = forecast_error1(jp_rv1,jp_midas)
peer_table[9,5] = fore_err$mse
peer_table[9,6] = fore_err$qlike

xtable(peer_table,digits=4)



loss_f <- function(true, pred){  
  # true = Winsorize(true, quantile(true, probs = c(0.01, 0.99), na.rm = FALSE))
  out_sampleT = length(true)
  mse = (true-pred)^2
  qlike = log(pred) + true/pred
  mse = mse[-which(mse>quantile(mse,0.99))]
  return(list(mse=mse, qlike=qlike))
}
###################
#### DM test
###################

dm_tab=matrix(NA,45,3)

dm_tab[1,1]=dm.test(garch_eu_mse,gjr_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[2,1]=dm.test(garch_eu_mse,x_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[3,1]=dm.test(garch_eu_mse,fi_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[4,1]=dm.test(garch_eu_mse,har_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[5,1]=dm.test(garch_eu_mse,harx_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[6,1]=dm.test(garch_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[7,1]=dm.test(garch_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[8,1]=dm.test(garch_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[9,1]=dm.test(gjr_eu_mse,x_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[10,1]=dm.test(gjr_eu_mse,fi_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[11,1]=dm.test(gjr_eu_mse,har_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[12,1]=dm.test(gjr_eu_mse,harx_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[13,1]=dm.test(gjr_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[14,1]=dm.test(gjr_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[15,1]=dm.test(gjr_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[16,1]=dm.test(x_eu_mse,fi_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[17,1]=dm.test(x_eu_mse,har_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[18,1]=dm.test(x_eu_mse,harx_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[19,1]=dm.test(x_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[20,1]=dm.test(x_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[21,1]=dm.test(x_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[22,1]=dm.test(fi_eu_mse,har_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[23,1]=dm.test(fi_eu_mse,harx_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[24,1]=dm.test(fi_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[25,1]=dm.test(fi_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[26,1]=dm.test(fi_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[27,1]=dm.test(har_eu_mse,harx_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[28,1]=dm.test(har_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[29,1]=dm.test(har_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[30,1]=dm.test(har_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[31,1]=dm.test(harx_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[32,1]=dm.test(harx_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[33,1]=dm.test(harx_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[34,1]=dm.test(rvarma_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[35,1]=dm.test(rvarma_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[36,1]=dm.test(rgarch_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[37,1]=dm.test(x_mfpca_eu_mse,garch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[38,1]=dm.test(x_mfpca_eu_mse,gjr_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[39,1]=dm.test(x_mfpca_eu_mse,x_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[40,1]=dm.test(x_mfpca_eu_mse,fi_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[41,1]=dm.test(x_mfpca_eu_mse,har_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[42,1]=dm.test(x_mfpca_eu_mse,harx_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[43,1]=dm.test(x_mfpca_eu_mse,rvarma_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[44,1]=dm.test(x_mfpca_eu_mse,rgarch_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]
dm_tab[45,1]=dm.test(x_mfpca_eu_mse,midas_eu_mse, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)[[1]]


xtable(dm_tab)


load("peer_eu.RData")
load("peer_uk.RData")
load("peer_jp.RData")
load("rgarch_eu.RData")
load("rgarch_uk.RData")
load("rgarch_jp.RData")
x<-1:674
plot(x,eu_rv1,col='black', cex = 1.5, pch=1, ylim=c(0,3.5),main = "", xlab="Date",cex.lab = 1.5,ylab="",xaxt="n")
axis(1,at=ceiling(seq(1,674,134.5)),labels=c("17/05/2017","16/01/2018","06/09/2018","14/05/2019","17/01/2020", "30/09/2020"),cex.axis=1.5)
# lines(g_mfpca_eu,type = "b", col = "magenta1", lwd = 2, pch = 2, cex = 1)
lines(x_mfpca_eu,type = "b", col = "red", lwd = 2, pch = 4, cex = 1)
lines(eu_garch,col='blue',lwd=2)
lines(eu_gjr,col='cyan',lwd=2)
lines(eu_x,col='darkgoldenrod1',lwd=2)
lines(eu_fi,col='chartreuse',lwd=2)
lines(eu_har,col='darkgrey',lwd=2)
lines(eu_harx,col='darkgreen',lwd=2)
lines(eu_rvarma,col='darkorchid1',lwd=2)
lines(eu_rgarch,col='deeppink',lwd=2)
lines(eu_midas,col='deepskyblue',lwd=2)
legend(x = c(20, 20), y = c(3.5, 3.5), legend = c("Realised Volatility","X-MFPCA","GARCH","GJR", "GARCH-X","FIGARCH","HAR","HAR-X","HAR-ARFIMA","RGARCH","GARCH-MIDAS"), lty=c(NA,NA,1,1,1,1,1,1,1,1,1),pch=c(1,4,NA,NA,NA,NA,NA,NA,NA,NA,NA), lwd=c(3,3,3,3,3,3,3,3,3,3,3),col = c('black','red','blue','cyan','darkgoldenrod1','chartreuse','darkgrey','darkgreen','darkorchid1','deeppink','deepskyblue'),y.intersp=1.3, cex = 1.2)


x<-1:674
plot(x,uk_rv1,col='black', cex = 1.5, pch=1, ylim=c(0,3.5),main = "", xlab="Date",cex.lab = 1.5,ylab="",xaxt="n")
axis(1,at=ceiling(seq(1,674,134.5)),labels=c("17/05/2017","16/01/2018","06/09/2018","14/05/2019","17/01/2020", "30/09/2020"),cex.axis=1.5)
# lines(g_mfpca_uk,type = "b", col = "red", lwd = 2, pch = 4, cex = 1)
lines(x_mfpca_uk,type = "b", col = "red", lwd = 2, pch = 4, cex = 1)
lines(uk_garch,col='blue',lwd=2)
lines(uk_gjr,col='cyan',lwd=2)
lines(uk_x,col='darkgoldenrod1',lwd=2)
lines(uk_fi,col='chartreuse',lwd=2)
lines(uk_har,col='darkgrey',lwd=2)
lines(uk_harx,col='darkgreen',lwd=2)
lines(uk_rvarma,col='darkorchid1',lwd=2)
lines(uk_rgarch,col='deeppink',lwd=2)
lines(uk_midas,col='deepskyblue',lwd=2)
legend(x = c(20, 20), y = c(3.5, 3.5), legend = c("Realised Volatility","X-MFPCA","GARCH","GJR", "GARCH-X","FIGARCH","HAR","HAR-X","HAR-ARFIMA","RGARCH","GARCH-MIDAS"), lty=c(NA,NA,1,1,1,1,1,1,1,1,1),pch=c(1,4,NA,NA,NA,NA,NA,NA,NA,NA,NA), lwd=c(3,3,3,3,3,3,3,3,3,3,3),col = c('black','red','blue','cyan','darkgoldenrod1','chartreuse','darkgrey','darkgreen','darkorchid1','deeppink','deepskyblue'),y.intersp=1.3, cex = 1.2)


x<-1:674
plot(x,jp_rv1,col='black', cex = 1.5, pch=1, ylim=c(0,3.5),main = "", xlab="Date",cex.lab = 1.5,ylab="",xaxt="n")
axis(1,at=ceiling(seq(1,674,134.5)),labels=c("17/05/2017","16/01/2018","06/09/2018","14/05/2019","17/01/2020", "30/09/2020"),cex.axis=1.5)
# lines(g_mfpca_jp,type = "b", col = "red", lwd = 2, pch = 4, cex = 1)
lines(x_mfpca_jp,type = "b", col = "red", lwd = 2, pch = 4, cex = 1)
lines(jp_garch,col='blue',lwd=2)
lines(jp_gjr,col='cyan',lwd=2)
lines(jp_x,col='darkgoldenrod1',lwd=2)
lines(jp_fi,col='chartreuse',lwd=2)
lines(jp_har,col='darkgrey',lwd=2)
lines(jp_harx,col='darkgreen',lwd=2)
lines(jp_rvarma,col='darkorchid1',lwd=2)
lines(jp_rgarch,col='deeppink',lwd=2)
lines(jp_midas,col='deepskyblue',lwd=2)
legend(x = c(20, 20), y = c(3.5, 3.5), legend = c("Realised Volatility","X-MFPCA","GARCH","GJR", "GARCH-X","FIGARCH","HAR","HAR-X","HAR-ARFIMA","RGARCH","GARCH-MIDAS"), lty=c(NA,NA,1,1,1,1,1,1,1,1,1),pch=c(1,4,NA,NA,NA,NA,NA,NA,NA,NA,NA), lwd=c(3,3,3,3,3,3,3,3,3,3,3),col = c('black','red','blue','cyan','darkgoldenrod1','chartreuse','darkgrey','darkgreen','darkorchid1','deeppink','deepskyblue'),y.intersp=1.3, cex = 1.2)


