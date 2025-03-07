#####################
#### basis methods
#####################

##################################################
### tfpca -Truncated functional principal component analysis 
##################################################
fpca.est <- function(yd){
  grid_point=nrow(yd)
  N=ncol(yd)
  yd=yd-rowMeans(yd)

    times=rep(0,grid_point)
    for(i in 1:grid_point){times[i]=i/grid_point}
    basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=32,norder=4)
    y_sq_fdsmooth=smooth.basis(argvals=times,y=yd,fdParobj=basis_obj)
    y_sq_fd=y_sq_fdsmooth$fd
    pca_obj=pca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd), centerfns=TRUE)  #centerfns=FALSE.
    eigen_functs=pca_obj$harmonics
    eigen_v=pca_obj$values
    tve=eigen_v/sum(eigen_v)
    ortho_basis_matrix=matrix(0,nrow=grid_point,ncol=10)
    for(j in 1:10){indi_basis=eval.fd(times,eigen_functs[j])
    ortho_basis_matrix[,j]=indi_basis}

    ncomp = 1; gap = 1
    while(gap<0.001)
    {
      gap = 1-(eigen_v[ncomp+1]/eigen_v[1])/(eigen_v[ncomp]/eigen_v[1])
      ncomp = ncomp+1
    }

    basis1 = as.matrix(ortho_basis_matrix[,1:ncomp])
    for (j in 1:ncomp){
      if(all(basis1[,j]<=0)==TRUE){
        basis1[,j] = -basis1[,j]
    }}
    basis1[basis1<0]<-0
  return(list(basis = basis1, tve = tve))
}

##################################################
## dfpca - dynamic functional principal component analysis
##################################################

# weight function (flat top)
wt.flat <- function(x, c){
    if ( -1< x & x <= -c)
    {
        return(x/(1-c) + 1/(1-c))
    }
    if (-c< x & x < c)
    {
        return(1)
    }
    if ( c< x & x <1)
    {
        return(x/(c-1) - 1/(c-1))
    }
    else return(0)
}

# weight function (point)
wt.point <- function(x){
    if ( -1< x & x < 1)
    {
        return(1-abs(x))
    }
    else return(0)
}

dfpca <- function(x, q, weight = c('flat','point'), c = 0.5){
  # xs <- scale(x, scale = FALSE)
  xs=x-rowMeans(x)
  n <- ncol(x) # sample size
  m <- nrow(x) # function dimension
  c1 <- array(0, dim = c(n, m, m)) # c1[h,,] is the covariance at h-1
  c2 <- array(0, dim = c(n, m, m)) # c2[h,,] is the covariance at 1-h
  for(h in 1:n)
  {
    for( k in h: n)
    {
      c1[h,,] <- c1[h,,] + xs[,k]%*%t(xs[,k-h+1])
      c2[h,,] <- c2[h,,] + xs[,k-h+1]%*%t(xs[,k])
    }
    c1[h,,] <- c1[h,,]/n
    c2[h,,] <- c2[h,,]/n
  }
  
  f <- array(0, dim = c(m, m)) # f is weighted sum of the covariance at all lags
  if (weight == 'flat'){
    for (h in 2:n)
    {
      f <- f + wt.flat(h/q, c)*c1[h,,] + wt.flat(-h/q, c)*c2[h,,]
    }
  }
  if (weight == 'point')
  {
    for (h in 2:n)
    {
      f <- f + wt.point(h/q)*c1[h,,] + wt.point(-h/q)*c2[h,,]
    }
  }
  f <- f + c1[1,,]
  md <- eigen(f, symmetric = T)

  md$values[md$values<0] <- 0
  eigen_value = md$values
  ncomp = 1; gap = 1
  while(gap<0.001)
    {
      gap = 1-(eigen_value[ncomp+1]/eigen_value[1])/(eigen_value[ncomp]/eigen_value[1])
      ncomp = ncomp+1
    }

  score <- t(xs)%*%md$vectors[, 1:ncomp]
  # fitted <- t(score%*%t(md$vectors[, 1:ncomp])) + apply(x, 2, mean)
  varprop <- eigen_value/sum(eigen_value) # cumulated

  basis <- as.matrix(md$vectors[, 1:ncomp])
  for (j in 1:ncomp){
    if(all(basis[,j]<=0)==TRUE){
      basis[,j] = -basis[,j]
  }}

  # basis <- cbind(apply(x, 2, mean), md$vectors[, 1:ncomp])
  # non-negativity
  basis[basis<0]<-0

  return(list(basis = basis, varprop =varprop))
}

##########################################################
## lfpca - long-range functional principal component analysis
##########################################################

LW_estimate <- function(fun_dat, band_const = 1, choice_nbasis = c("tve", "ratio"), 
                       est_method = c("classical", "modified"))
{
    est_method = match.arg(est_method)
    T = ncol(fun_dat)
    m = min(T - 1, nrow(fun_dat)) * band_const

    X_bar = rowMeans(fun_dat, na.rm = TRUE)
    center_dat = sweep(fun_dat, 1, X_bar)
   
    # calculating long-run covariance for a given lag
    
    gamma_l <- function(lag, T)
    {
        gamma_lag_sum = 0
        if(lag >= 0)
        {
            for(ij in 1:(T-lag))
            {
                gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij+lag)])))
            }
        }
        else
        {
            for(ij in 1:(T - abs(lag)))
            {
                gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij + abs(lag)]) %*%
                                t(as.matrix(center_dat[,ij])))
            }
        }
        return(gamma_lag_sum/T)
    }
   
    # calculating the sum of long-run covariance for all lags (m equals the total number of grid points)
    C = 0
    index = seq(-m, m, by = 1)
    for(k in 1:length(index))
    {
        C = C + (m - abs(index[k])) * gamma_l(lag = index[k], T = T)
    }
    
    # eigen decomposition
    eigen_decomp = eigen(C)
    
    # select the number of component based on 90% total amount of variation
    eigen_value = eigen_decomp$values
    
    if(choice_nbasis == "tve")
    {
        prop_eigen_value = eigen_value/sum(eigen_value)
        ncomp = head(which(prop_eigen_value >= 0.90), 1)
    }
    if(choice_nbasis == "ratio")
    {   
        prop_eigen_value = eigen_value/sum(eigen_value)
        ncomp = 1; gap = 1
        while(gap<0.001)
        {
           gap = 1-(eigen_value[ncomp+1]/eigen_value[1])/(eigen_value[ncomp]/eigen_value[1])
           ncomp = ncomp+1
        }
    }

    basis = as.matrix(eigen_decomp$vectors[,1:ncomp])
    for (j in 1:ncomp){
      if(all(basis[,j]<=0)==TRUE){
        basis[,j] = -basis[,j]
    }}

    # first functional principal component
    eigen_function = matrix(eigen_decomp$vectors[,1], nrow = 1)
    score = as.numeric(eigen_function %*% fun_dat)
    
    # local whittle estimator
    # implementing the exact local Whittle estimator of Shimotsu and Phillips (2005) that is consistent and asymptotically normal as long as the optimization range is less than 9/2, so that it is possible to estimate the memory of stationary as well as non-stationary processes.
    delta = 0.7 # a tunning parameter set 0< delta < 1.
    d_est = ELW(score, floor(1+T^delta))$d # floor(1+T^delta) - bandwith parameter
    
    return(list(d_est = d_est, ncomp = ncomp, varprop = prop_eigen_value, basis = basis))
}


########################################################
## mfpca - Multilevel functional principal component analysis
########################################################

Multilevel_FPCA <-function(data){
  grid_point = dim(data)[1]
  sample_size = dim(data)[2]
  multi_d = dim(data)[3]

  # common mean mu_hat
    com_mean = rowMeans(apply(data, c(1,3),mean)) #rowMeans(data[,,1])

  # entity-specific mean eta_hat_i
    en_mean = apply(data, c(1,3),mean) - com_mean

  # common trend R_hat
    com_R = apply(data, c(1,2),mean) - com_mean

  # entity-specific trend U_hat_i
    en_U = array(NA, dim = c(grid_point, sample_size, multi_d))
    for (i in 1:multi_d){
      en_U[,,i] = data[,,i] - rowMeans(data[,,i]) - apply(data, c(1,2),mean) + com_mean
    }

  # covariance and eigen-decomposition analysis
    eigen_decomp_en_mean = eigen(cov(t(en_mean)))
    lambda_en_mean = eigen_decomp_en_mean$values
    tve_en_mean = lambda_en_mean/sum(lambda_en_mean)
    basis_en_com = eigen_decomp_en_mean$vectors[,1:10]

  # long-run covariance and eigen-decomposition analysis
    eigen_decomp_com_R = eigen(long_run_covariance_estimation(com_R))
    lambda_com_R = eigen_decomp_com_R$values
    tve_com_R = lambda_com_R/sum(lambda_com_R)
    basis_com_R = eigen_decomp_com_R$vectors[,1:10]
    
    tve_en_U = matrix(NA, 10, multi_d)
    basis_en_U = array(NA, dim = c(grid_point, 10, multi_d))


    for (i in 1:multi_d){
      eigen_decomp_en_U = eigen(long_run_covariance_estimation(en_U[,,i]))
      lambda_en_U = eigen_decomp_en_U$values
      tve = lambda_en_U/sum(lambda_en_U)
      tve_en_U[,i] = tve[1:10]

      enuv = eigen_decomp_en_U$vectors[,1:10]
      for (j in 1:10){
        if(all(enuv[,j]<=0)==TRUE){
          enuv[,j] = -enuv[,j]
      }}
      enuv[enuv<0]<-0

      basis_en_U[,,i] = enuv
    }
    
    # number of basis
    # entity-specific mean
    en_mean_ncomp = 1; gap = 1
    while(gap<0.001)
    {
      gap = 1-(tve_en_mean[en_mean_ncomp+1]/tve_en_mean[1])/(tve_en_mean[en_mean_ncomp]/tve_en_mean[1])
      en_mean_ncomp = en_mean_ncomp+1
    }
    
    # common trend R_hat
    com_R_ncomp = 1; gap = 1
    while(gap<0.001)
    {
      gap = 1-(tve_com_R[com_R_ncomp+1]/tve_com_R[1])/(tve_com_R[com_R_ncomp]/tve_com_R[1])
      com_R_ncomp = com_R_ncomp+1
    }

    # entity-specific trend
    en_U_ncomp = 1; gap = 1
    # should be checked by 3 entities. here use the default 1.
    

    # non-negative
    for (j in 1:10){
      if(all(basis_en_com[,j]<=0)==TRUE){
        basis_en_com[,j] = -basis_en_com[,j]
    }}
    basis_en_com[basis_en_com<0]<-0

    for (j in 1:10){
      if(all(basis_com_R[,j]<=0)==TRUE){
        basis_com_R[,j] = -basis_com_R[,j]
    }}
    basis_com_R[basis_com_R<0]<-0

  return(list(basis_en_com = basis_en_com, basis_com_R = basis_com_R, basis_en_U = basis_en_U, tve_en_mean = tve_en_mean, tve_com_R = tve_com_R, tve_en_U = tve_en_U, en_mean_ncomp = en_mean_ncomp, com_R_ncomp = com_R_ncomp, en_U_ncomp = en_U_ncomp))
}


