# function to load

#################
### VaR backtest
#################
#' Forecast Daily/intra-day Value-at-Risk
#' var.forecast function forecasts the daily VaR and intra-day VaR curves according to intra-day return curves.
var.forecast <- function(sigma_pred,error_fit,quantile_v,Method){
  # a function to resample curve data with replacement.
  resmp <- function(data,boot_N,boot_in_N){
    bootres=matrix(NA,1,boot_N)
    for (j in 1:boot_N){
      bootres[,j]=data[sample(c(1:boot_in_N),1,replace=TRUE)]
    }
    return(bootres)}

  point_grid=nrow(sigma_pred)
  in_N=ncol(error_fit)
  switch(Method,
         normal = {
           error_var=qnorm(p=quantile_v,mean=0,sd=1)
           var_curves=sigma_pred*error_var
         },
         bootstrap = {
            bootN=1000
            error_var = c(rep(NA,point_grid))
            for (j in 1:point_grid){
               bootres=resmp(error_fit[j,],bootN,in_N)
               error_var[j]=quantile(bootres,quantile_v,na.rm = TRUE)
               if(is.infinite(error_var[j])==TRUE){error_var[j]=error_var[j-1]}
            }
            var_curves=sigma_pred*error_var
         },
         stop("Enter something to switch me!")
  )
  return(var_curves)
}

var.forecast1 <- function(sigma_pred,error_quantile){
  var_curves=sweep(sigma_pred, 1, error_quantile, "*")#sigma_pred*error_quantile
  return(var_curves)
}

#' Violation Process
#' var.vio function returns a violation process for the intra-day VaR curves.
var.vio <- function(yd,var_curve){
  if( ncol(yd)!=ncol(var_curve) ) warning('the number of observations must match for calculating the violation process!')
  point_grid=nrow(yd)
  n=ncol(yd)
  y_vio=matrix(0,point_grid,n)
  for (i in 1:n){
    ind=as.matrix(yd)[,i]<as.matrix(var_curve)[,i]
    y_vio[ind,i]=1
  }
  return(y_vio)
}

var.vio_short <- function(yd,var_curve){
  if( ncol(yd)!=ncol(var_curve) ) warning('the number of observations must match for calculating the violation process!')
  point_grid=nrow(yd)
  n=ncol(yd)
  y_vio=matrix(0,point_grid,n)
  for (i in 1:n){
    ind=as.matrix(yd)[,i]>as.matrix(var_curve)[,i]
    y_vio[ind,i]=1
  }
  return(y_vio)
}

#' Backtest Intra-day VaR forecasts
#' var.backtest function backtests the unbiasedness and the independence hypotheses for the intra-day VaR curve forecasts.
#'
#' vio A (grid_point) x (number of observations) matrix drawn from the violation process curves.
#' tau The nominal/true quantile of the VaR curves.
#' K The maximal lagged autocorrelation considered for the independence test. If it is missing, a default value "K=20" is used.
var.backtest <- function(vio, tau, K=NULL){
  #########################
  ### unbiasedness test ###
  #########################
  # the test T statistics
  H0_Tn<-function(z,tau){
    point_grid=nrow(z)
    n=ncol(z)
    int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}
    Tn=n*int_approx(rowMeans(z)-tau)^2
    return(Tn)
  }
  # get the critical values
  Tn_cv<-function(z,cv_N){
    N=ncol(z)
    point_grid=nrow(z)

    # autocovariance operator from -N to N
    cov_est<-function(x,h){
      N=ncol(x)
      c_est=x[,1:(N-h)]%*%t(x[,(h+1):N])/(N-h)
      return(c_est)
    }

    Gam=0
    Gam=Gam+cov_est(z,0)#long_run_cov2(z, C0 = 3, H = 3)#c

    eig=eigen(Gam)
    eigvals=eig$values
    eigvals=as.vector(eigvals)
    vect=eigvals/point_grid

    int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}

    lim_sum=matrix(0,cv_N,1)

    lim_sum=0
    for (j in 1:length(vect)){
      lim_sum=lim_sum+vect[j]*rnorm(cv_N,mean=0,sd=1)^2}

    cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))
    return(list(cv,lim_sum))
  }

  cv_N=5000 # sample size for approximating critical values
  Tn_stat=H0_Tn(vio,tau)
  limit=Tn_cv(vio,cv_N)
  emplim=limit[[2]]
  unbias_p = 1-ecdf(emplim)(Tn_stat)

  #########################
  ### independent test ###
  #########################
  if(is.null(K) == TRUE) {
    K=20
  }

  T_statistic <- function(fdata, lag_val){
    T = nrow(fdata)
    p = ncol(fdata)
    # calculate autocovariance function
    gamma_l <- function(fdata, lag){
      T = nrow(fdata)
      center_dat = t(scale(fdata, center = TRUE, scale = FALSE))

      gamma_lag_sum = 0
      for(ij in 1:(T - lag)){
        gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij + lag)])))
      }
      return(gamma_lag_sum/T)
    }
    gamma_hat = gamma_l(fdata = fdata, lag = lag_val)
    T_h = T * ((1/p)^2) * sum(gamma_hat^2)
    return(T_h)
  }

  # Portmanteau test statistics
  gaPort.stat <- function(H, datmat){
    vals = rep(0,H)
    for(j in 1:H){
      vals[j] = T_statistic(t(datmat), j)
    }
    return(sum(vals))
  }

  # useful inner functions
  eta <- function(i, j, datmat){
    T = dim(datmat)[2]
    p = dim(datmat)[1]
    l = max(c(i,j))
    sum1 = array(0, c(p, p, p, p))
    for(k in (l+1):T){
      sum1 = sum1 + datmat[,k-i] %o% datmat[,k] %o% datmat[,k-j] %o% datmat[,k]
    }
    return(sum1/T)
  }

  etaM <- function(i, datmat){
    T = dim(datmat)[2]
    p = dim(datmat)[1]
    l = i
    sum1 = array(0, c(p,p))
    for(k in (l+1):T){
      sum1 = sum1 + (datmat[,k-i])^2 %o% (datmat[,k])^2
    }
    return(sum1/T)
  }

  etapopM <- function(H, datmat){
    etal = list()
    for(j in 1:H){
      etal[[j]] = list()
    }
    for(k in 1:H){
      etal[[k]] = etaM(k, datmat)
    }
    return(etal)
  }

  mean.W.2 <- function(x.e, H, datmat){
    sum1 = 0
    p = dim(datmat)[1]
    for(j in 1:H){
      sum1 = sum1 + sum(x.e[[j]])
    }
    return(sum1/p^2)
  }

  etapopNvarMC2 <- function(H, datmat, len1, len2){
    sum1 = 0

    rref = runif(len1, 0, 1)
    rref = sort(rref)
    rrefind = round(rref * dim(datmat)[1])
    rrefind[which(rrefind == 0)] = 1

    rrefG = c(0, rref)
    xd = diff(rrefG)
    gmat = xd %o% xd %o% xd %o% xd

    for(m in 1:H){
      x = eta(m, m, datmat[rrefind,])
      sum1 = sum1 + sum((x*x)*gmat)
    }

    rref2 = runif(len2, 0, 1)
    rref2 = sort(rref2)
    rrefind2 = round(rref2 * dim(datmat)[1])
    rrefind2[which(rrefind2 == 0)] = 1
    rrefG2 = c(0, rref2)
    xd2 = diff(rrefG2)
    gmat2 = xd2 %o% xd2 %o% xd2 %o% xd2

    if(H > 1){
      for(j in 1:(H-1)){
        for(k in (j+1):H){
          if((abs(k-j)) < 3){
            x = eta(j, k, datmat[rrefind,])
            #sum1=sum1+sum( (x.e[[j]][[k]])^2)
            sum1 = sum1 + 2 * sum((x*x)*gmat)
          }#end if

          if((abs(k-j)) >= 3){
            x = eta(j, k, datmat[rrefind2,])
            #sum1=sum1+sum( (x.e[[j]][[k]])^2)
            sum1 = sum1 + 2 * sum((x*x)*gmat2)
          }
        }
      }
    }
    return(2 * sum1)
  }

  vio = vio - rowMeans(vio)
  len1 = 20;len2 = 15; # parameters for Monte Carlo integral, can be altered.
  res = gaPort.stat(K, vio)
  x.e = etapopM(K, vio)
  res2 = mean.W.2(x.e, K, vio)
  res3 = etapopNvarMC2(K, vio, len1, len2)
  beta = res3/(2 * res2)
  nu = (2 * res2^2)/res3
  independent_p = 1 - pchisq(res/beta, df = nu)

  return(list(unbias = unbias_p, independent = independent_p))
}



########################################################
## var backtesting results

load("fgarch_tfpca_eu.RData")
eutrue_out = pred_eu[[3]]
eut_day = pred_eu[[1]] # eu [[1]], uk [[2]], jp[[3]]
eut_intrvol = pred_eu[[2]]
eut_err_quantile025 = pred_eu[[4]] #0.0025
eut_err_quantile01 = pred_eu[[5]] #0.01
eut_err_quantile99 = pred_eu[[6]]

load("fgarch_dpca_eu.RData")
eud_vol_pred = pred_eu[[1]]
eud_intrvol = pred_eu[[2]]
eud_err_quantile025 = pred_eu[[4]]
eud_err_quantile01 = pred_eu[[5]]
eud_err_quantile99 = pred_eu[[6]]

load("fgarch_lpca_eu.RData")
eul_vol_pred = pred_eu[[1]]
eul_intrvol = pred_eu[[2]]
eul_err_quantile025 = pred_eu[[4]]
eul_err_quantile01 = pred_eu[[5]]
eul_err_quantile99 = pred_eu[[6]]

load("fgarch_mfpca_eu.RData")
eum_vol_pred = pred_eu[[1]]
eum_intrvol = pred_eu[[2]]
eum_err_quantile025 = pred_eu[[4]]
eum_err_quantile01 = pred_eu[[5]]
eum_err_quantile99 = pred_eu[[6]]


load("fgarch_tfpca_uk.RData")
uktrue_out = pred_uk[[3]]
ukt_day = pred_uk[[1]] # uk [[1]], uk [[2]], jp[[3]]
ukt_intrvol = pred_uk[[2]]
ukt_err_quantile025 = pred_uk[[4]] #0.0025
ukt_err_quantile01 = pred_uk[[5]] #0.01
ukt_err_quantile99 = pred_uk[[6]]

load("fgarch_dpca_uk.RData")
ukd_vol_pred = pred_uk[[1]]
ukd_intrvol = pred_uk[[2]]
ukd_err_quantile025 = pred_uk[[4]]
ukd_err_quantile01 = pred_uk[[5]]
ukd_err_quantile99 = pred_uk[[6]]

load("fgarch_lpca_uk.RData")
ukl_vol_pred = pred_uk[[1]]
ukl_intrvol = pred_uk[[2]]
ukl_err_quantile025 = pred_uk[[4]]
ukl_err_quantile01 = pred_uk[[5]]
ukl_err_quantile99 = pred_uk[[6]]

load("fgarch_mfpca_uk.RData")
ukm_vol_pred = pred_uk[[1]]
ukm_intrvol = pred_uk[[2]]
ukm_err_quantile025 = pred_uk[[4]]
ukm_err_quantile01 = pred_uk[[5]]
ukm_err_quantile99 = pred_uk[[6]]



load("fgarch_tfpca_jp.RData")
jptrue_out = pred_jp[[3]]
jpt_day = pred_jp[[1]] # jp [[1]], uk [[2]], jp[[3]]
jpt_intrvol = pred_jp[[2]]
jpt_err_quantile025 = pred_jp[[4]] #0.0025
jpt_err_quantile01 = pred_jp[[5]] #0.01
jpt_err_quantile99 = pred_jp[[6]]

load("fgarch_dpca_jp.RData")
jpd_vol_pred = pred_jp[[1]]
jpd_intrvol = pred_jp[[2]]
jpd_err_quantile025 = pred_jp[[4]]
jpd_err_quantile01 = pred_jp[[5]]
jpd_err_quantile99 = pred_jp[[6]]

load("fgarch_lpca_jp.RData")
jpl_vol_pred = pred_jp[[1]]
jpl_intrvol = pred_jp[[2]]
jpl_err_quantile025 = pred_jp[[4]]
jpl_err_quantile01 = pred_jp[[5]]
jpl_err_quantile99 = pred_jp[[6]]

load("fgarch_mfpca_jp.RData")
jpm_vol_pred = pred_jp[[1]]
jpm_intrvol = pred_jp[[2]]
jpm_err_quantile025 = pred_jp[[4]]
jpm_err_quantile01 = pred_jp[[5]]
jpm_err_quantile99 = pred_jp[[6]]





load("x_tfpca_eu.RData")
euxtrue_out = pred_eu[[3]]
euxt_day = pred_eu[[1]] # eu [[1]], uk [[2]], jp[[3]]
euxt_intrvol = pred_eu[[2]]
euxt_err_quantile025 = pred_eu[[4]] #0.0025
euxt_err_quantile01 = pred_eu[[5]] #0.01
euxt_err_quantile99 = pred_eu[[6]]

load("x_dpca_eu.RData")
euxd_vol_pred = pred_eu[[1]]
euxd_intrvol = pred_eu[[2]]
euxd_err_quantile025 = pred_eu[[4]]
euxd_err_quantile01 = pred_eu[[5]]
euxd_err_quantile99 = pred_eu[[6]]

load("x_lpca_eu.RData")
euxl_vol_pred = pred_eu[[1]]
euxl_intrvol = pred_eu[[2]]
euxl_err_quantile025 = pred_eu[[4]]
euxl_err_quantile01 = pred_eu[[5]]
euxl_err_quantile99 = pred_eu[[6]]

load("x_mfpca_eu.RData")
euxm_vol_pred = pred_eu[[1]]
euxm_intrvol = pred_eu[[2]]
euxm_err_quantile025 = pred_eu[[4]]
euxm_err_quantile01 = pred_eu[[5]]
euxm_err_quantile99 = pred_eu[[6]]


load("x_tfpca_uk.RData")
ukxtrue_out = pred_uk[[3]]
ukxt_day = pred_uk[[1]] # uk [[1]], uk [[2]], jp[[3]]
ukxt_intrvol = pred_uk[[2]]
ukxt_err_quantile025 = pred_uk[[4]] #0.0025
ukxt_err_quantile01 = pred_uk[[5]] #0.01
ukxt_err_quantile99 = pred_uk[[6]]

load("x_dpca_uk.RData")
ukxd_vol_pred = pred_uk[[1]]
ukxd_intrvol = pred_uk[[2]]
ukxd_err_quantile025 = pred_uk[[4]]
ukxd_err_quantile01 = pred_uk[[5]]
ukxd_err_quantile99 = pred_uk[[6]]

load("x_lpca_uk.RData")
ukxl_vol_pred = pred_uk[[1]]
ukxl_intrvol = pred_uk[[2]]
ukxl_err_quantile025 = pred_uk[[4]]
ukxl_err_quantile01 = pred_uk[[5]]
ukxl_err_quantile99 = pred_uk[[6]]

load("x_mfpca_uk.RData")
ukxm_vol_pred = pred_uk[[1]]
ukxm_intrvol = pred_uk[[2]]
ukxm_err_quantile025 = pred_uk[[4]]
ukxm_err_quantile01 = pred_uk[[5]]
ukxm_err_quantile99 = pred_uk[[6]]



load("x_tfpca_jp.RData")
jpxtrue_out = pred_jp[[3]]
jpxt_day = pred_jp[[1]] # jp [[1]], uk [[2]], jp[[3]]
jpxt_intrvol = pred_jp[[2]]
jpxt_err_quantile025 = pred_jp[[4]] #0.0025
jpxt_err_quantile01 = pred_jp[[5]] #0.01
jpxt_err_quantile99 = pred_jp[[6]]

load("x_dpca_jp.RData")
jpxd_vol_pred = pred_jp[[1]]
jpxd_intrvol = pred_jp[[2]]
jpxd_err_quantile025 = pred_jp[[4]]
jpxd_err_quantile01 = pred_jp[[5]]
jpxd_err_quantile99 = pred_jp[[6]]

load("x_lpca_jp.RData")
jpxl_vol_pred = pred_jp[[1]]
jpxl_intrvol = pred_jp[[2]]
jpxl_err_quantile025 = pred_jp[[4]]
jpxl_err_quantile01 = pred_jp[[5]]
jpxl_err_quantile99 = pred_jp[[6]]

load("x_mfpca_jp.RData")
jpxm_vol_pred = pred_jp[[1]]
jpxm_intrvol = pred_jp[[2]]
jpxm_err_quantile025 = pred_jp[[4]]
jpxm_err_quantile01 = pred_jp[[5]]
jpxm_err_quantile99 = pred_jp[[6]]



table1 = matrix(NA, 15, 8)

# var025 = var.forecast1(sqrt(eut_intrvol),eut_err_quantile025)
# vio_025 = var.vio(eutrue_out,var025)
# mean(rowMeans(vio_025))
# table1[1,1] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[2,1] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[3,1] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[4,1] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[5,1] = var.backtest(vio_025, tau=0.0025, K=20)$independent

# var025 = var.forecast1(sqrt(eud_intrvol),eud_err_quantile025)
# vio_025 = var.vio(eutrue_out,var025)
# table1[1,2] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[2,2] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[3,2] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[4,2] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[5,2] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(eul_intrvol),eul_err_quantile025)
vio_025 = var.vio(eutrue_out,var025)
table1[1,3] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[2,3] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[3,3] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[4,3] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[5,3] = var.backtest(vio_025, tau=0.0025, K=20)$independent


var025 = var.forecast1(sqrt(eum_intrvol),eum_err_quantile025)
vio_025 = var.vio(eutrue_out,var025)
table1[1,4] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[2,4] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[3,4] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[4,4] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[5,4] = var.backtest(vio_025, tau=0.0025, K=20)$independent


# var025 = var.forecast1(sqrt(euxt_intrvol),euxt_err_quantile025)
# vio_025 = var.vio(eutrue_out,var025)
# table1[1,5] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[2,5] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[3,5] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[4,5] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[5,5] = var.backtest(vio_025, tau=0.0025, K=20)$independent

# var025 = var.forecast1(sqrt(euxd_intrvol),euxd_err_quantile025)
# vio_025 = var.vio(eutrue_out,var025)
# table1[1,6] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[2,6] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[3,6] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[4,6] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[5,6] = var.backtest(vio_025, tau=0.0025, K=20)$independent
# mean(rowMeans(vio_025))
var025 = var.forecast1(sqrt(euxl_intrvol),euxl_err_quantile025)
vio_025 = var.vio(eutrue_out,var025)
table1[1,7] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[2,7] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[3,7] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[4,7] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[5,7] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(euxm_intrvol),euxm_err_quantile025)
vio_025 = var.vio(eutrue_out,var025)
table1[1,8] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[2,8] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[3,8] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[4,8] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[5,8] = var.backtest(vio_025, tau=0.0025, K=20)$independent





# var025 = var.forecast1(sqrt(ukt_intrvol),ukt_err_quantile025)
# vio_025 = var.vio(uktrue_out,var025)
# # mean(rowMeans(vio_025))
# table1[6,1] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[7,1] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[8,1] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[9,1] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[10,1] = var.backtest(vio_025, tau=0.0025, K=20)$independent

# var025 = var.forecast1(sqrt(ukd_intrvol),ukd_err_quantile025)
# vio_025 = var.vio(uktrue_out,var025)
# table1[6,2] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[7,2] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[8,2] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[9,2] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[10,2] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(ukl_intrvol),ukl_err_quantile025)
vio_025 = var.vio(uktrue_out,var025)
table1[6,3] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[7,3] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[8,3] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[9,3] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[10,3] = var.backtest(vio_025, tau=0.0025, K=20)$independent


var025 = var.forecast1(sqrt(ukm_intrvol),ukm_err_quantile025)
vio_025 = var.vio(uktrue_out,var025)
table1[6,4] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[7,4] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[8,4] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[9,4] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[10,4] = var.backtest(vio_025, tau=0.0025, K=20)$independent


# var025 = var.forecast1(sqrt(ukxt_intrvol),ukxt_err_quantile025)
# vio_025 = var.vio(uktrue_out,var025)
# table1[6,5] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[7,5] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[8,5] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[9,5] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[10,5] = var.backtest(vio_025, tau=0.0025, K=20)$independent

# var025 = var.forecast1(sqrt(ukxd_intrvol),ukxd_err_quantile025)
# vio_025 = var.vio(uktrue_out,var025)
# table1[6,6] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[7,6] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[8,6] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[9,6] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[10,6] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(ukxl_intrvol),ukxl_err_quantile025)
vio_025 = var.vio(uktrue_out,var025)
table1[6,7] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[7,7] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[8,7] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[9,7] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[10,7] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(ukxm_intrvol),ukxm_err_quantile025)
vio_025 = var.vio(uktrue_out,var025)
table1[6,8] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[7,8] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[8,8] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[9,8] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[10,8] = var.backtest(vio_025, tau=0.0025, K=20)$independent




# var025 = var.forecast1(sqrt(jpt_intrvol),jpt_err_quantile025)
# vio_025 = var.vio(jptrue_out,var025)
# # mean(rowMeans(vio_025))
# table1[11,1] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[12,1] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[13,1] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[14,1] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[15,1] = var.backtest(vio_025, tau=0.0025, K=20)$independent

# var025 = var.forecast1(sqrt(jpd_intrvol),jpd_err_quantile025)
# vio_025 = var.vio(jptrue_out,var025)
# table1[11,2] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[12,2] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[13,2] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[14,2] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[15,2] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(jpl_intrvol),jpl_err_quantile025)
vio_025 = var.vio(jptrue_out,var025)
table1[11,3] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[12,3] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[13,3] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[14,3] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[15,3] = var.backtest(vio_025, tau=0.0025, K=20)$independent


var025 = var.forecast1(sqrt(jpm_intrvol),jpm_err_quantile025)
vio_025 = var.vio(jptrue_out,var025)
table1[11,4] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[12,4] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[13,4] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[14,4] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[15,4] = var.backtest(vio_025, tau=0.0025, K=20)$independent
#

# var025 = var.forecast1(sqrt(jpxt_intrvol),jpxt_err_quantile025)
# vio_025 = var.vio(jptrue_out,var025)
# table1[11,5] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[12,5] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[13,5] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[14,5] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[15,5] = var.backtest(vio_025, tau=0.0025, K=20)$independent

# var025 = var.forecast1(sqrt(jpxd_intrvol),jpxd_err_quantile025)
# vio_025 = var.vio(jptrue_out,var025)
# table1[11,6] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
# table1[12,6] = var.backtest(vio_025, tau=0.0025, K=1)$independent
# table1[13,6] = var.backtest(vio_025, tau=0.0025, K=5)$independent
# table1[14,6] = var.backtest(vio_025, tau=0.0025, K=10)$independent
# table1[15,6] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(jpxl_intrvol),jpxl_err_quantile025)
vio_025 = var.vio(jptrue_out,var025)
table1[11,7] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[12,7] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[13,7] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[14,7] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[15,7] = var.backtest(vio_025, tau=0.0025, K=20)$independent

var025 = var.forecast1(sqrt(jpxm_intrvol),jpxm_err_quantile025)
vio_025 = var.vio(jptrue_out,var025)
table1[11,8] = var.backtest(vio_025, tau=0.0025, K=1)$unbias
table1[12,8] = var.backtest(vio_025, tau=0.0025, K=1)$independent
table1[13,8] = var.backtest(vio_025, tau=0.0025, K=5)$independent
table1[14,8] = var.backtest(vio_025, tau=0.0025, K=10)$independent
table1[15,8] = var.backtest(vio_025, tau=0.0025, K=20)$independent


xtable(round(table1,2))


#######################################################
####### save the following for the VaR-correct trading strategy
#######################################################
g_var_lfpca_eu01 = var.forecast1(sqrt(eul_intrvol),eul_err_quantile01)
g_var_mfpca_eu01 = var.forecast1(sqrt(eum_intrvol),eum_err_quantile01)
x_var_lfpca_eu01 = var.forecast1(sqrt(euxl_intrvol),euxl_err_quantile01)
x_var_mfpca_eu01 = var.forecast1(sqrt(euxm_intrvol),euxm_err_quantile01)


g_var_lfpca_uk01 = var.forecast1(sqrt(ukl_intrvol),ukl_err_quantile01)
g_var_mfpca_uk01 = var.forecast1(sqrt(ukm_intrvol),ukm_err_quantile01)
x_var_lfpca_uk01 = var.forecast1(sqrt(ukxl_intrvol),ukxl_err_quantile01)
x_var_mfpca_uk01 = var.forecast1(sqrt(ukxm_intrvol),ukxm_err_quantile01)

 
g_var_lfpca_jp01 = var.forecast1(sqrt(jpl_intrvol),jpl_err_quantile01)
g_var_mfpca_jp01 = var.forecast1(sqrt(jpm_intrvol),jpm_err_quantile01)
x_var_lfpca_jp01 = var.forecast1(sqrt(jpxl_intrvol),jpxl_err_quantile01)
x_var_mfpca_jp01 = var.forecast1(sqrt(jpxm_intrvol),jpxm_err_quantile01)


save( g_var_lfpca_eu01, g_var_mfpca_eu01, file = "eu_forecast_intraday_g.RData")
save( x_var_lfpca_eu01, x_var_mfpca_eu01, file = "eu_forecast_intraday_x.RData")

save( g_var_lfpca_uk01, g_var_mfpca_uk01, file = "uk_forecast_intraday_g.RData")
save( x_var_lfpca_uk01, x_var_mfpca_uk01, file = "uk_forecast_intraday_x.RData")

save( g_var_lfpca_jp01, g_var_mfpca_jp01, file = "jp_forecast_intraday_g.RData")
save( x_var_lfpca_jp01, x_var_mfpca_jp01, file = "jp_forecast_intraday_x.RData")




g_var_lfpca_eu99 = var.forecast1(sqrt(eul_intrvol),eul_err_quantile99)
g_var_mfpca_eu99 = var.forecast1(sqrt(eum_intrvol),eum_err_quantile99)
x_var_lfpca_eu99 = var.forecast1(sqrt(euxl_intrvol),euxl_err_quantile99)
x_var_mfpca_eu99 = var.forecast1(sqrt(euxm_intrvol),euxm_err_quantile99)


g_var_lfpca_uk99 = var.forecast1(sqrt(ukl_intrvol),ukl_err_quantile99)
g_var_mfpca_uk99 = var.forecast1(sqrt(ukm_intrvol),ukm_err_quantile99)
x_var_lfpca_uk99 = var.forecast1(sqrt(ukxl_intrvol),ukxl_err_quantile99)
x_var_mfpca_uk99 = var.forecast1(sqrt(ukxm_intrvol),ukxm_err_quantile99)

 
g_var_lfpca_jp99 = var.forecast1(sqrt(jpl_intrvol),jpl_err_quantile99)
g_var_mfpca_jp99 = var.forecast1(sqrt(jpm_intrvol),jpm_err_quantile99)
x_var_lfpca_jp99 = var.forecast1(sqrt(jpxl_intrvol),jpxl_err_quantile99)
x_var_mfpca_jp99 = var.forecast1(sqrt(jpxm_intrvol),jpxm_err_quantile99)


save( g_var_lfpca_eu99, g_var_mfpca_eu99, file = "eu_forecast_intraday_g99.RData")
save( x_var_lfpca_eu99, x_var_mfpca_eu99, file = "eu_forecast_intraday_x99.RData")

save( g_var_lfpca_uk99, g_var_mfpca_uk99, file = "uk_forecast_intraday_g99.RData")
save( x_var_lfpca_uk99, x_var_mfpca_uk99, file = "uk_forecast_intraday_x99.RData")

save( g_var_lfpca_jp99, g_var_mfpca_jp99, file = "jp_forecast_intraday_g99.RData")
save( x_var_lfpca_jp99, x_var_mfpca_jp99, file = "jp_forecast_intraday_x99.RData")

