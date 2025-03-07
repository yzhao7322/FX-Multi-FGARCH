###################
### functions to load
###################
int_approx=function(x){
        temp_n=NROW(x)
        return((1/temp_n)*sum(x))}

norm_fd<-function(x){
    x=as.matrix(x)
    M=ncol(x)
      int_approx <- function(x){
      temp_n = NROW(x)
      return((1/(temp_n)) * sum(x))}
    normx=c(rep(NA,M))
    for (i in 1:M){
      normx[i]=sqrt(int_approx(x[,i]^2))
    }
  return(sum(normx))
}

blockresmp<-function(data,boot_N){
   P=nrow(data)
   N=ncol(data)
   bootres=matrix(NA,P,boot_N)
        for (j in 1:boot_N){
          bootres[,j]=data[,sample(c(1:N),1,replace=TRUE)]
        }
        return(bootres)
}

win_sigma <- function(dy, sigma2fit, delta_hat, alpha_Op_hat, beta_Op_hat){
  point_grid=nrow(dy)
  sigrep = ncol(dy)
  sigma2_fit=matrix(NA,point_grid,(sigrep+1) )
  sigma2_fit[,1]=sigma2fit

  for(i in 1:sigrep){
        for(j in 1:point_grid){ 
            fit_alpha_op = alpha_Op_hat[j,] * (dy[,i])^2
            fit_beta_op = beta_Op_hat[j,] * (sigma2_fit[,i])
            sigma2_fit[j,i+1] = delta_hat[j] + int_approx(fit_alpha_op) + int_approx(fit_beta_op)
        }
  }
  return(sigma2_fit[,2: (sigrep+1) ])
}

win_sigma <- function(dy, sigma2fit, x_fit, delta_hat, alpha_Op_hat, beta_Op_hat, gamma_Op_hat){
  point_grid = nrow(dy)
  sigrep = ncol(dy)
  sigma2_fit=matrix(NA,point_grid,(sigrep+1) )
  sigma2_fit[,1]=sigma2fit

  for(i in 1:sigrep){
        for(j in 1:point_grid){ 
            fit_alpha_op = alpha_Op_hat[j,] * (dy[,i])^2
            fit_beta_op = beta_Op_hat[j,] * sigma2_fit[,i]
            fit_x_op = gamma_Op_hat[j,] * x_fit[,i]
            sigma2_fit[j,i+1] = delta_hat[j] + int_approx(fit_alpha_op) + int_approx(fit_beta_op) + int_approx(fit_x_op)
        }     
  }
  return(sigma2_fit[,2:(sigrep+1) ])
}


###################
### forecasting fgarch(1,1)
###################
# cidr_data : objective intraday return data, with rows presenting the intraday grid points, and columns presenting the sample observations.
# fpcatype : the type of FPCA, either "tfpca", "dfpca" or "lfpca".
forecast_FGARCH<-function(cidr_data, fpcatype){
  point_grid=nrow(cidr_data)
  N=ncol(cidr_data) # the same size - 1274
  win_N = 600 # training sample size
  out_N = N-win_N
  out_win = 4 # forecasting horizon

  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,N,1) 
  intrvol=matrix(NA,point_grid,N)
  err_quantile025=matrix(NA,point_grid,out_N) 
  err_quantile025[point_grid,]=0
  err_quantile01=matrix(NA,point_grid,out_N) 
  err_quantile01[point_grid,]=0
  err_quantile99=matrix(NA,point_grid,out_N) 
  err_quantile99[point_grid,]=0
  true_out=cidr_data[,(win_N+1):N]

  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){###re-estimate the model every 42 observations
    sub_yd = cidr_data[,insample_ind[i,1]:insample_ind[i,2]]
    
    if (fpcatype =="tfpca"){
    	# static FPCA
        basis_est = fpca.est(sub_yd^2)$basis
    } else if (fpcatype =="dfpca"){
        # # dynamic FPCA
        q = floor(log(dim(sub_yd)[2])) #floor(sqrt(dim(sub_yd)[2]))### here I changed
        basis_est = dfpca(sub_yd^2, q, weight = 'flat')$basis
    } else if (fpcatype =="lfpca"){
        # # long-memory FPCA
        basis_est = LW_estimate(sub_yd^2, choice_nbasis = "ratio")$basis
    }
    
    for (j in 1:ncol(basis_est)){
      basis_est[,j]=basis_est[,j]/norm_fd(basis_est[,j])
    }

    # FGARCH(1,1)
    fgarch_fit=fgarch_est(sub_yd,basis_est)

    delta_hat = fgarch_fit[[2]]
    alpha_Op_hat = fgarch_fit[[3]]
    beta_Op_hat = fgarch_fit[[4]]

    error_fit = fgarch_fit[[6]]
    bootres = blockresmp(error_fit,10000)
    sigma500 = fgarch_fit[[1]][,win_N]

    yfit=fgarch_fit[[7]]

    sample_42 = cidr_data[,( outsample_ind[i]:(outsample_ind[i+1]-1) )]
    sigma2_42 = win_sigma(sample_42, sigma500, delta_hat, alpha_Op_hat, beta_Op_hat)

    for (j in 1:out_winv[i]){
      intrvol[,(insample_ind[i,2]+j - win_N) ] = sigma2_42[,j]
      vol_pred[ (insample_ind[i,2]+j - win_N) ] = sigma2_42[193,j]
    
      for (u in 1:(point_grid-1)){
        err_quantile025[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.0025,na.rm = TRUE)
        err_quantile01[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.01,na.rm = TRUE)
        err_quantile99[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.99,na.rm = TRUE)
      }
    }
  }
  return(list(vol_pred,intrvol,true_out,err_quantile025,err_quantile01,err_quantile99))
}
# output - vol_pred: predicted inter-daily conditional volatility
#          intrvol: predicted intraday conditional volatility
#          true_out: true OCIDR curves from out-of-sample
#          err_quantile025: 0.25% quantile of bootstrapped error
#          err_quantile01: 1% quantile of bootstrapped error
#          err_quantile99: 99% quantile of bootstrapped error 


# cidr_data : objective intraday return data, with rows presenting the intraday grid points, and columns presenting the sample observations.
# cidr_data2 : other intraday return data that helps to derive the common factor, with rows presenting the intraday grid points, and columns presenting the sample observations.
# cidr_data3 : other intraday return data that helps to derive the common factor, with rows presenting the intraday grid points, and columns presenting the sample observations.
forecast_multi_FGARCH<-function(cidr_data,cidr_data2,cidr_data3){
  point_grid=nrow(cidr_data)
  N=ncol(cidr_data) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,N,1) 
  intrvol=matrix(NA,point_grid,N)
  err_quantile025=matrix(NA,point_grid,out_N) 
  err_quantile025[point_grid,]=0
  err_quantile01=matrix(NA,point_grid,out_N) 
  err_quantile01[point_grid,]=0
  err_quantile99=matrix(NA,point_grid,out_N) 
  err_quantile99[point_grid,]=0
  true_out = cidr_data2[,(win_N+1):N]
  
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){###re-estimate the model every 42 observations
    sub_yd = cidr_data[,insample_ind[i,1]:insample_ind[i,2]]
    sub_yd2 = cidr_data2[,insample_ind[i,1]:insample_ind[i,2]]
    sub_yd3 = cidr_data3[,insample_ind[i,1]:insample_ind[i,2]]

    # multi-level FPCA
    sub_N = ncol(sub_yd)
    obj_dat = array(NA, dim = c(point_grid, sub_N, 3))
    obj_dat[,,1] = sub_yd^2
    obj_dat[,,2] = sub_yd2^2
    obj_dat[,,3] = sub_yd3^2

    basis_est = Multilevel_FPCA(obj_dat)
    basis_com = basis_est$basis_com_R 
    basis_res = basis_est$basis_en_U
    ncomp_R = basis_est$com_R_ncomp
    ncomp_U = basis_est$en_U_ncomp
    
    basis_res1 = as.matrix(basis_res[,1:1,1]) # eu 1:ncomp_en_U[1]
    basis_res2 = as.matrix(basis_res[,1:1,2]) # uk
    basis_res3 = as.matrix(basis_res[,1:1,3]) # jp
    
    basis_res_obj = basis_res2

    for (j in 1:ncol(basis_com)){
      basis_com[,j]=basis_com[,j]/norm_fd(basis_com[,j])
    }
    for (j in 1:ncol(basis_res_obj)){
      basis_res_obj[,j]=basis_res_obj[,j]/norm_fd(basis_res_obj[,j])
    }

    basis = cbind(basis_com[,1:1],basis_res_obj)
    # FGARCH(1,1)
    fgarch_fit=fgarch_est(sub_yd2,basis)
  
    delta_hat = fgarch_fit[[2]]
    alpha_Op_hat = fgarch_fit[[3]]
    beta_Op_hat = fgarch_fit[[4]]

    error_fit = fgarch_fit[[6]]
    bootres = blockresmp(error_fit,10000)
    sigma500 = fgarch_fit[[1]][,win_N]

    yfit=fgarch_fit[[7]]
    
    sample_42 = cidr_data2[,( outsample_ind[i]:(outsample_ind[i+1]-1) )]
    sigma2_42 = win_sigma(sample_42, sigma500, delta_hat, alpha_Op_hat, beta_Op_hat)

    for (j in 1:out_winv[i]){
      intrvol[,(insample_ind[i,2]+j - win_N) ] = sigma2_42[,j]
      vol_pred[ (insample_ind[i,2]+j - win_N) ] = sigma2_42[point_grid,j]
    
      for (u in 1:(point_grid-1)){
        err_quantile025[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.0025,na.rm = TRUE)
        err_quantile01[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.01,na.rm = TRUE)
        err_quantile99[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.99,na.rm = TRUE)
      }
    }
  }
  return(list(vol_pred,intrvol,true_out,err_quantile025,err_quantile01,err_quantile99))
}


############################
## forecasting with fgarch-x
############################

forecast_X<-function(cidr_data, x_data, fpcatype){
  point_grid=nrow(cidr_data)
  N=ncol(cidr_data) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,N,1) 
  intrvol=matrix(NA,point_grid,N)
  err_quantile025=matrix(NA,point_grid,out_N) 
  err_quantile025[point_grid,]=0
  err_quantile01=matrix(NA,point_grid,out_N) 
  err_quantile01[point_grid,]=0
  err_quantile99=matrix(NA,point_grid,out_N) 
  err_quantile99[point_grid,]=0
  true_out=cidr_data[,(win_N+1):N]

  out_winv = diff(outsample_ind)
  
  for (i in 1:mdl_rep){###re-estimate the model every 42 observations
    sub_yd = cidr_data[,insample_ind[i,1]:insample_ind[i,2]]
    sub_xd = x_data[,insample_ind[i,1]:insample_ind[i,2]]
    
    if (fpcatype =="tfpca"){
    	# static FPCA
        basis_est = fpca.est(sub_yd^2)$basis
    } else if (fpcatype =="dfpca"){
        # # dynamic FPCA
        q = floor(log(dim(sub_yd)[2]))
        basis_est = dfpca(sub_yd^2, q, weight = 'flat')$basis
    } else if (fpcatype =="lfpca"){
        # # long-memory FPCA
        basis_est = LW_estimate(sub_yd^2, choice_nbasis = "ratio")$basis
    }

    for (j in 1:ncol(basis_est)){
      basis_est[,j]=basis_est[,j]/norm_fd(basis_est[,j])
    }
  
    x_fit=fgarchx_est(sub_yd,sub_xd,basis_est)
  
    delta_hat = x_fit[[2]]
    Op_hat = x_fit[[3]]
    alpha_Op_hat = Op_hat[[1]]
    beta_Op_hat = Op_hat[[2]]
    gamma_Op_hat = Op_hat[[2]]
  
    error_fit = x_fit[[5]]
    bootres = blockresmp(error_fit,10000)
    sigma500 = x_fit[[1]][,win_N]

    yfit=x_fit[[6]]

    # window sample input, e.g., first window: 500 - 541, to forecast 501 - 542.
    sample_42 = cidr_data[, ( outsample_ind[i]:(outsample_ind[i+1]-1) ) ]
    x_42      = x_data[, ( outsample_ind[i]:(outsample_ind[i+1]-1) ) ]
    sigma2_42 = win_sigma(sample_42, sigma500, x_42, delta_hat, alpha_Op_hat, beta_Op_hat, gamma_Op_hat)

    for (j in 1:out_winv[i]){
      intrvol[,(insample_ind[i,2]+j - win_N) ] = sigma2_42[,j]
      vol_pred[ (insample_ind[i,2]+j - win_N) ] = sigma2_42[point_grid,j]
    
      for (u in 1:(point_grid-1)){
        err_quantile025[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.0025,na.rm = TRUE)
        err_quantile01[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.01,na.rm = TRUE)
        err_quantile99[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.99,na.rm = TRUE)
      }
    }

  }
  return(list(vol_pred,intrvol,true_out,err_quantile025,err_quantile01,err_quantile99))
}



forecast_multi_X<-function(cidr_data, cidr_data2, cidr_data3, x_data){
  point_grid=nrow(cidr_data)
  N=ncol(cidr_data) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,N,1) 
  intrvol=matrix(NA,point_grid,N)
  err_quantile025=matrix(NA,point_grid,out_N) 
  err_quantile025[point_grid,]=0
  err_quantile01=matrix(NA,point_grid,out_N) 
  err_quantile01[point_grid,]=0
  err_quantile99=matrix(NA,point_grid,out_N) 
  err_quantile99[point_grid,]=0

  true_out=cidr_data[,(win_N+1):N]

  out_winv = diff(outsample_ind)
  
  for (i in 1:mdl_rep){###re-estimate the model every 42 observations
    sub_yd = cidr_data[,insample_ind[i,1]:insample_ind[i,2]]
    sub_yd2 = cidr_data2[,insample_ind[i,1]:insample_ind[i,2]]
    sub_yd3 = cidr_data3[,insample_ind[i,1]:insample_ind[i,2]]
    sub_xd = x_data[,insample_ind[i,1]:insample_ind[i,2]]

    # multi-level FPCA
    sub_N = ncol(sub_yd)
    obj_dat = array(NA, dim = c(point_grid, sub_N, 3))
    obj_dat[,,1] = sub_yd^2
    obj_dat[,,2] = sub_yd2^2
    obj_dat[,,3] = sub_yd3^2

    basis_est = Multilevel_FPCA(obj_dat)
    basis_com = basis_est$basis_com_R 
    basis_res = basis_est$basis_en_U
    
    basis_res1 = as.matrix(basis_res[,1:1,1]) # eu
    basis_res2 = as.matrix(basis_res[,1:1,2]) # uk
    basis_res3 = as.matrix(basis_res[,1:1,3]) # jp
    basis_res_obj = basis_res1

    for (j in 1:ncol(basis_com)){
      basis_com[,j]=basis_com[,j]/norm_fd(basis_com[,j])
    }
    for (j in 1:ncol(basis_res_obj)){
      basis_res_obj[,j]=basis_res_obj[,j]/norm_fd(basis_res_obj[,j])
    }

    basis = cbind(basis_com[,1:1],basis_res_obj)

    x_fit=fgarchx_est(sub_yd,sub_xd,basis)
  
    delta_hat = x_fit[[2]]
    Op_hat = x_fit[[3]]
    alpha_Op_hat = Op_hat[[1]]
    beta_Op_hat = Op_hat[[2]]
    gamma_Op_hat = Op_hat[[2]]
  
    error_fit = x_fit[[5]]
    bootres = blockresmp(error_fit,10000)
    sigma500 = x_fit[[1]][,win_N]

    yfit=x_fit[[6]]
  
    # window sample input, e.g., first window: 500 - 541, to forecast 501 - 542.
    sample_42 = cidr_data[, ( outsample_ind[i]:(outsample_ind[i+1]-1) ) ] ##### change here for eu, jp, uk!!!
    x_42      = x_data[, ( outsample_ind[i]:(outsample_ind[i+1]-1) ) ]
    sigma2_42 = win_sigma(sample_42, sigma500, x_42, delta_hat, alpha_Op_hat, beta_Op_hat, gamma_Op_hat)

    for (j in 1:out_winv[i]){
      intrvol[,(insample_ind[i,2]+j - win_N) ] = sigma2_42[,j]
      vol_pred[ (insample_ind[i,2]+j - win_N) ] = sigma2_42[point_grid,j]
    
      for (u in 1:(point_grid-1)){
        err_quantile025[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.0025,na.rm = TRUE)
        err_quantile01[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.01,na.rm = TRUE)
        err_quantile99[u,(insample_ind[i,2]+j - win_N) ]=quantile(bootres[u,],0.99,na.rm = TRUE)
      }
    }

  }
  return(list(vol_pred,intrvol,true_out,err_quantile025,err_quantile01,err_quantile99))
}

