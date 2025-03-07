
project_values<-function(xt,M){
  U=nrow(xt)
    N=ncol(xt)
    xt=xt-rowMeans(xt)

    times=rep(0,U)
    for(i in 1:U){times[i]=i/U}
    basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=20,norder=4)  #original: nbasis=64, norder=6
    x_fdsmooth=smooth.basis(argvals=times,y=xt,fdParobj=basis_obj)
    x_fd=x_fdsmooth$fd
    pca_obj=pca.fd(x_fd,nharm=10,harmfdPar=fdPar(x_fd), centerfns=FALSE)  #centerfns=TRUE may be better
    eigen_values=pca_obj$values
    return(eigen_values)
}

project<-function(xt,M){
  U=nrow(xt)
    N=ncol(xt)
    xt=xt-rowMeans(xt)

    times=rep(0,U)
    for(i in 1:U){times[i]=i/U}
    basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=20,norder=4)  #original: nbasis=64, norder=6
    x_fdsmooth=smooth.basis(argvals=times,y=xt,fdParobj=basis_obj)
    x_fd=x_fdsmooth$fd
    pca_obj=pca.fd(x_fd,nharm=10,harmfdPar=fdPar(x_fd), centerfns=FALSE)  #centerfns=TRUE may be better
    eigen_functs=pca_obj$harmonics
    ortho_basis_matrix=matrix(0,nrow=U,ncol=10)
    for(j in 1:10){ortho_basis_matrix[,j]=eval.fd(times,eigen_functs[j])}
    
    int_approx=function(x){
        temp_n=NROW(x)
        return((1/temp_n)*sum(x))}
    # This will hold the y^2 projections onto phi_1,...,phi_M
    x_proj_coefs1=matrix(0,nrow=M,ncol=N)
    xi=matrix(NA,M,N)
    for(i in 1:N){
            for(j in 1:M){xi[j,i]=int_approx(xt[,i]*ortho_basis_matrix[,j])}}
  return(xi)
}


# normality test
normality_test<-function(data,M){
    N=ncol(data)
    xi=t(project(data,M))

    ##### T2 test
    lambda=project_values(data,M)
    Z_rv=matrix(NA,N,M)
    U_rv=matrix(NA,N,1)
    for (j in 1:M){
        Z_rv[,j]=xi[,j]/sqrt(lambda[j])
        U_rv=cbind(U_rv,Z_rv[,j]^3)
        U_rv=cbind(U_rv,Z_rv[,j]^4)
    }
    U_rv=U_rv[,2:(2*M+1)]
    U_mu=t(as.matrix(c(rep(c(0,3),M))))

    U_demean=matrix(NA,N,2*M)
    for (i in 1:N){
    U_demean[i,]=U_rv[i,]-U_mu}
    
    lrcov=lrvar(U_demean,type="Newey-West",prewhite=FALSE)
    
    T2=1/N*t(colSums(U_demean))%*%inv(lrcov)%*%colSums(U_demean)
    P_T2=1-pchisq(T2,df=2*M)
    return(list(T2,P_T2))
}

# stationary test
T_stationary1 <- function (sample, L = 49, J = 500, MC_rep = 1000, cumulative_var = 0.9, 
    Ker1 = FALSE, Ker2 = TRUE, h = ncol(sample)^0.5, pivotal = FALSE) 
{
    xrefine = N = ncol(sample)
    trefine = nrow(sample)
    if (Ker1) {
        K = function(x) {
            output = min(1, max(1.1 - abs(x), 0))
            return(output)
        }
    }
    if (Ker2) {
        K = function(x) {
            output = min(1, max(2 - 2 * abs(x), 0))
            return(output)
        }
    }
    basis = create.fourier.basis(c(0, 1), L)
    ld = length(cumulative_var)
    h1 = h
    X1_bar = rowMeans(sample)
    mean_subtracted_X1 = sample - X1_bar
    gamma_hat = list()
    for (i in 0:(N - 1)) {
        temp_matrix = matrix(rep(0, trefine^2), trefine, trefine)
        for (j in (i + 1):N) {
            temp_matrix = temp_matrix + (mean_subtracted_X1[, 
                j] %*% t(mean_subtracted_X1[, j - i]))
        }
        gamma_hat[i + 1] = list(temp_matrix/N)
    }
    cov_sample1 = gamma_hat[[1]]
    for (index in 1:(N - 1)) {
        cov_sample1 = cov_sample1 + K(index/h1) * (gamma_hat[[index + 
            1]] + t(gamma_hat[[index + 1]]))
    }
    Z_matrix = cov_sample1
    e1 = list()
    for (index in 1:L) {
        e1[index] = list(as.matrix(eval.basis(evalarg = (1:trefine)/trefine, 
            basisobj = basis, Lfdobj = 0)[, index]))
    }
    eigenvalues1 = (eigen(Z_matrix)$values)/trefine
    D = matrix(0, L, L)
    for (k in 1:L) {
        for (ell in 1:L) {
            Integrand_matrix = Z_matrix * (e1[[k]] %*% t(e1[[ell]]))
            D[k, ell] = 1/(trefine^2) * sum(Integrand_matrix)
        }
    }
    eigenpairs = eigen(D)
    eigenvectors = eigenpairs$vec
    eigenvectors2 = eigen(Z_matrix)$vectors
    eigenvalues = eigenpairs$val
    evals = eigenvalues
    if (pivotal) {
        d = c(1:ld)
        switch = 0
        stoper = 1
        spot = 1
        while (switch == 0) {
            while ((sum(eigenvalues[c(1:spot)])/sum(eigenvalues)) < 
                cumulative_var[stoper]) {
                spot = spot + 1
            }
            d[stoper] = spot
            stoper = stoper + 1
            if (stoper == (length(d) + 1)) {
                switch = 1
            }
        }
        T_N0 = 1:ld
        for (r in 1:ld) {
            ds = d[r]
            inp.matrix = matrix(0, ds, N)
            eig.v.norm = ((trefine)^0.5) * eigenvectors2
            for (j in (1:ds)) {
                for (k in (1:N)) {
                  inp.matrix[j, k] = t(sample[, k]) %*% (eig.v.norm[, 
                    j])/trefine
                }
            }
            T_Nsum = rep(0, ds)
            for (j in (1:ds)) {
                s.0 = sum(inp.matrix[j, (1:xrefine)])
                for (x in (1:xrefine)) {
                  T_Nsum[j] = T_Nsum[j] + (1/xrefine) * ((1/N^0.5) * 
                    (sum(inp.matrix[j, (1:x)]) - (x/xrefine) * 
                      s.0))^2
                }
            }
            T_N0[r] = sum(T_Nsum/eigenvalues[1:ds])
        }
        T = vector(, MC_rep)
        T_array = matrix(0, length(d), MC_rep)
        lambda = eigenvalues
        for (dd in 1:length(d)) {
            for (k in 1:MC_rep) {
                z = rnorm(d[dd] * J)
                tot = 0
                for (n in c(1:d[dd])) {
                  sum1 = 0
                  sum1 = sum((z[c(((n - 1) * d[dd] + 1):((n - 
                    1) * d[dd] + J))]/(pi * c(1:J)))^2)
                  tot = tot + sum1
                }
                T_array[dd, k] = T[k] = tot
            }
        }
        p_values = (1:ld)
        for (dd in 1:length(d)) {
            p_values[dd] = round(1 - ecdf(T_array[dd, ])(T_N0[dd]), 
                4)
        }
    }
    if (pivotal == FALSE) {
        T = vector(, MC_rep)
        T_array = (1:MC_rep)
        lambda = eigenvalues
        d = min(c(length(which(lambda > 0)), 15))
        for (k in 1:MC_rep) {
            z = rnorm(d * J)
            tot = 0
            for (n in c(1:d)) {
                sum1 = 0
                sum1 = sum((z[c(((n - 1) * d + 1):((n - 1) * 
                  d + J))]/(pi * c(1:J)))^2)
                tot = tot + lambda[n] * sum1
            }
            T_array[k] = T[k] = tot
        }
        int = sum(((1/sqrt(N)) * ((sample[, 1]) - (1/xrefine) * 
            rowSums(sample)))^2/(xrefine * trefine))
        for (x in 2:xrefine) {
            int = int + sum(((1/sqrt(N)) * (rowSums(sample[, 
                1:x]) - (x/xrefine) * rowSums(sample)))^2/(xrefine * 
                trefine))
        }
        statT_N = int
        p_values = (1:ld)
        for (dd in 1:length(d)) {
            p_values = round(1 - ecdf(T_array)(statT_N), 4)
        }
    }
    return(p_values)
}


# Portmanteau test statistics 
##################
# Test statistics
##################

T_statistic <- function(fdata, lag_val)
{
    # calculate sample size and dimensionality
  
    T = nrow(fdata)
    p = ncol(fdata)
  
    # calculate autocovariance function
  
    gamma_l <- function(fdata, lag)
    {
        T = nrow(fdata)
        center_dat = t(scale(fdata, center = TRUE, scale = FALSE))
    
        gamma_lag_sum = 0
        for(ij in 1:(T - lag))
        {
            gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij + lag)])))
        }
        return(gamma_lag_sum/T)
    }
    gamma_hat = gamma_l(fdata = fdata, lag = lag_val)
    T_h = T * ((1/p)^2) * sum(gamma_hat^2)
    return(T_h)
}    

################### idependent tests
# Portmanteau test statistics 

gaPort.stat <- function(H, datmat)
{
    vals = rep(0,H)
    for(j in 1:H)
    {
        vals[j] = T_statistic(t(datmat), j)
    }
    return(sum(vals))
}

eta <- function(i, j, datmat)
{
    T = dim(datmat)[2]
    p = dim(datmat)[1]
    l = max(c(i,j))
    sum1 = array(0, c(p, p, p, p))
    for(k in (l+1):T)
    {
        sum1 = sum1 + datmat[,k-i] %o% datmat[,k] %o% datmat[,k-j] %o% datmat[,k]
    }
    #if(i<=j){
    #sum2=(datmat[,(1+j-i):(T-i)]%*%t(datmat[,(1+j):T]))%o%(datmat[,(1+j):T]%*%t(datmat[,1:(T-j)]))
    #}
    return(sum1/T)
}

etaM <- function(i, datmat)
{
    T = dim(datmat)[2]
    p = dim(datmat)[1]
    l = i
    sum1 = array(0, c(p,p))
    for(k in (l+1):T)
    {
        sum1 = sum1 + (datmat[,k-i])^2 %o% (datmat[,k])^2
    }
    #if(i<=j){
    #sum2=(datmat[,(1+i):(T-i)]%*%t(datmat[,(1+j):T]))%o%(datmat[,(1+j):T]%*%t(datmat[,1:(T-j)]))
    #}
    return(sum1/T)
}


etapop <- function(H, datmat){
    etal = list()
    for(j in 1:H)
    {
        etal[[j]] = list()
    }
    for(k in 1:H)
    {
        for(i in 1:H)
        {
            etal[[k]][[i]] = eta(k, i, datmat)
            #print(k)
        }
    }
    return(etal)
}

etapopN <- function(H, datmat)
{
    etal = list()
    for(j in 1:H)
    {
        etal[[j]]=list()
    }
    for(k in 1:H)
    {
        for(i in k:H)
        {
            etal[[k]][[i]] = eta(k, i, datmat)
            #print(k)
        }
    }

    for(k in 2:H)
    {
        for(i in (1:(k-1)))
        {
            etal[[k]][[i]] = etal[[i]][[k]]
        }
    }
    return(etal)
}

etapopNvar <- function(H, datmat, len)
{
    p = dim(datmat)[1]
    xe = floor((p/len)*(1:len))
    sum1 = 0

    for(j in 1:H)
    {
        x = eta(j, j, datmat[xe,])
        sum1 = sum1 + sum((x*x))
    }

    if(H > 1)
    {
        for(j in 1:(H-1))
        {
            for(k in (j+1):H)
            {
                x = eta(j, k, datmat[xe,])
                #sum1=sum1+sum( (x.e[[j]][[k]])^2)
                sum1 = sum1 + 2 * sum((x*x))
            }
        }
    }#end if H>1

    sum1 = 2*sum1/(len+1)^4
    return(sum1)
}


etapopM <- function(H, datmat)
{
    etal = list()
    for(j in 1:H)
    {
        etal[[j]] = list()
    }
    for(k in 1:H)
    {
        etal[[k]] = etaM(k, datmat)
        #print(k)
    }
    return(etal)
}




mean.W <- function(x.e, H, datmat)
{
    sum1 = 0
    p = dim(datmat)[1]
    for(j in 1:H)
    {
        for(k in 1:p)
        {
            sum1 = sum1 + sum(diag(x.e[[j]][[j]][,,k,k]))
        }
    }
    return(sum1/p^2)
}


mean.W.2 <- function(x.e, H, datmat)
{
    sum1 = 0
    p = dim(datmat)[1]
    for(j in 1:H)
    {
        sum1 = sum1 + sum(x.e[[j]])
    }
    return(sum1/p^2)
}

var.W <- function(x.e, H, datmat)
{
    sum1 = 0
    p = dim(datmat)[1]

    for(j in 1:H)
    {
        for(k in 1:H)
        {
            sum1 = sum1 + sum( (x.e[[j]][[k]])^2)
            #sum1=sum1+sum( (x.e[[j]][[k]])*(x.e[[j]][[k]]))
        }
    }
    sum1 = 2 * sum1/p^4
    return(sum1)
}

var.WN <- function(x.e, H, datmat)
{
    sum1 = 0
    p = dim(datmat)[1]
    for(j in 1:H)
    {
        sum1 = sum1 + sum((x.e[[j]][[j]])*(x.e[[j]][[j]]))
    }

    for(j in 1:(H-1))
    {
        for(k in (j+1):H)
        {
            #sum1=sum1+sum( (x.e[[j]][[k]])^2)
            sum1 = sum1 + 2 * sum((x.e[[j]][[k]])*(x.e[[j]][[k]]))
        }
    }
    sum1 = 2 * sum1/p^4
    return(sum1)
}

etapopNvarMC <- function(H, datmat, len)
{
    sum1 = 0

    rref = runif(len, 0, 1)
    rref = c(sort(rref), 1)
    rrefind = round(rref * dim(datmat)[1])
    rrefind[which(rrefind==0)] = 1

    rrefG = c(0, rref)
    xd = diff(rrefG)
    gmat = xd %o% xd %o% xd %o% xd

    p = dim(datmat)[1]
    for(j in 1:H)
    {
        x = eta(j, j, datmat[rrefind,])
        sum1 = sum1 + sum((x*x)*gmat)
    }

    if(H > 1)
    {
        for(j in 1:(H-1))
        {
            for(k in (j+1):H)
            {
                x = eta(j, k, datmat[rrefind,])
                #sum1=sum1+sum( (x.e[[j]][[k]])^2)
                sum1 = sum1 + 2 * sum((x*x)*gmat)
            }
        }
    }#end H if
    #sum1=2*sum1/p^4
    return(2 * sum1)
}


etapopNvarMC2 <- function(H, datmat, len1, len2)
{
    sum1 = 0

    rref = runif(len1, 0, 1)
    #rref=c(sort(rref),1)
    rref = sort(rref)
    rrefind = round(rref * dim(datmat)[1])
    rrefind[which(rrefind == 0)] = 1

    rrefG = c(0, rref)
    xd = diff(rrefG)
    gmat = xd %o% xd %o% xd %o% xd

    for(m in 1:H)
    {
        x = eta(m, m, datmat[rrefind,])
        sum1 = sum1 + sum((x*x)*gmat)
    }

    rref2 = runif(len2, 0, 1)
    #rref2=c(sort(rref2),1)
    rref2 = sort(rref2)
    rrefind2 = round(rref2 * dim(datmat)[1])
    rrefind2[which(rrefind2 == 0)] = 1
    rrefG2 = c(0, rref2)
    xd2 = diff(rrefG2)
    gmat2 = xd2 %o% xd2 %o% xd2 %o% xd2

    if(H > 1)
    {
        for(j in 1:(H-1))
        {
            for(k in (j+1):H)
            {
                if((abs(k-j)) < 3)
                {
                    x = eta(j, k, datmat[rrefind,])
                    #sum1=sum1+sum( (x.e[[j]][[k]])^2)
                    sum1 = sum1 + 2 * sum((x*x)*gmat)
                }#end if

                if((abs(k-j)) >= 3)
                {
                    x = eta(j, k, datmat[rrefind2,])
                    #sum1=sum1+sum( (x.e[[j]][[k]])^2)
                    sum1 = sum1 + 2 * sum((x*x)*gmat2)
                }
            }
        }
    }#end if H>1
    return(2 * sum1)
}

port.test.ga <- function(datmat, H, len1, len2)
{
    datmat = datmat - rowMeans(datmat)

    res = gaPort.stat(H, datmat)
    x.e = etapopM(H, datmat)
    res2 = mean.W.2(x.e, H, datmat)
    res3 = etapopNvarMC2(H, datmat, len1, len2)
    beta = res3/(2 * res2)
    nu = (2 * res2^2)/res3
    
    return(1 - pchisq(res/beta, df = nu))
}


#######################################0.05528557
# The Function^2-based Portmenteau Test
#######################################

func1<-function(fdata) ## fdata - each column contains a functional object
{ 
  TT=ncol(fdata)
  p=nrow(fdata)
  fdata=fdata-rowMeans(fdata)
  # fdata=scale(fdata,center = TRUE, scale = FALSE)
  y2=matrix(NA,p,TT)
  for(i in 1:TT){
     y2[,i]=fdata[,i]^2
   }
   y2=y2-rowMeans(y2)
return(y2)
}

func1_emp<-function(fdata) ## fdata - each column contains a functional object
{ 
  TT=ncol(fdata)
  p=nrow(fdata)
  fdata=fdata-rowMeans(fdata)
  # fdata=scale(fdata,center = TRUE, scale = FALSE)
  y2=matrix(NA,p,TT)
  for(i in 1:TT){
     y2[,i]=fdata[,i]^2
   }
   y2=y2-rowMeans(y2)
return(y2)
}

gamma2<-function(data, lag_val) ## data should be y2
{
  TT=ncol(data)
  p=nrow(data)
  # data=data-rowMeans(data)
  gam2_sum=0
  for (j in 1:(TT-lag_val)){
    gam2_sum=gam2_sum+data[,j]%*%t(data[,(j+lag_val)])
  }
  gam2=gam2_sum/TT
  return(gam2)
}

Port2.stat<-function(H,fdata) ## the portmenteau test based on the function2
{
  TT=ncol(fdata)
  p=nrow(fdata)
  vals=rep(0,H)
  #fdata=scale(fdata, center = TRUE, scale = FALSE)
  y2=func1_emp(fdata)
  # sigma2=gamma2(y2,0)
  for(h in 1:H){
    gam2hat=gamma2(y2,h)
    #rho=gam2hat/sigma2### rho - autocorrelation function
    vals[h]=TT*((1/p)^2)*sum(gam2hat^2)}
  stat=sum(vals)
  return(stat)
} 

etaM <- function(datmat)
{
  T = dim(datmat)[2]
  p = dim(datmat)[1]
  datmat=func1(datmat)
  sum1 = 0#array(0, c(p,p,p,p))#############
  for(k in 1:T)
  {
    sum1 = sum1 + (datmat[,k])^2# %o% (datmat[,k]) %o% (datmat[,k]) %o% (datmat[,k])
  }
  return(sum1/T)
}


mean.W.2 <- function(H, datmat,len1)
{
  sum1 = 0
  p = dim(datmat)[1]
  
  etal=etaM(datmat)
  sum1 = (sum(etal)/p)^2
  sum1=sum1*H

  # sum1 = H*trace.tensor(etal,c("I1","I2"),c("I3","I4"))
  return(sum1)
}

etaM_cov <- function(datmat)
{
  T = dim(datmat)[2]
  p = dim(datmat)[1]
  datmat=func1(datmat)
  sum1 = 0#array(0, c(p,p,p,p))#############
  for(k in 1:T)
  {
    sum1 = sum1 + (datmat[,k])%*%t(datmat[,k])# %o% (datmat[,k]) %o% (datmat[,k]) %o% (datmat[,k])
  }
  return(sum1/T)
}

etapopNvarMC2 <- function(H, datmat)
{
  T = dim(datmat)[2]
  p = dim(datmat)[1]
  
    x = etaM_cov(datmat)^2
    sum1=sum(x)/(p^2)
    sum1=2*H*(sum1^2)
    # x=etaM(datmat[rrefind,])^2
    # x=as.tensor(x^2)
    # sum1=2*H^2*trace.tensor(x,c("I1","I2"),c("I3","I4"))
  
  return(sum1)
}


port.test.2 <- function(datmat, H)
{ 
  res = Port2.stat(H, datmat)

  res2 = mean.W.2(H, datmat)
  res3 = etapopNvarMC2(H, datmat)
  beta = res3/(2 * res2)
  nu = (2 * res2^2)/res3
  
  return(1 - pchisq(res/beta, df = nu))
}



#  Local Whittle estimator

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
        prop_eigen_value = cumsum(eigen_value/sum(eigen_value))
        ncomp = head(which(prop_eigen_value >= 0.90), 1)
    }
    if(choice_nbasis == "ratio")
    {   
        prop_eigen_value = cumsum(eigen_value/sum(eigen_value))
        ncomp = 1; gap = 1
        while(gap<0.001)
        {
           gap = 1-(eigen_value[ncomp+1]/eigen_value[1])/(eigen_value[ncomp]/eigen_value[1])
           ncomp = ncomp+1
        }
    }
    
    basis = eigen_decomp$vectors[,1:ncomp]

    # first functional principal component
    eigen_function = matrix(eigen_decomp$vectors[,1], nrow = 1)
    score = as.numeric(eigen_function %*% fun_dat)
    
    # local whittle estimator
    # implementing the exact local Whittle estimator of Shimotsu and Phillips (2005) that is consistent and asymptotically normal as long as the optimization range is less than 9/2, so that it is possible to estimate the memory of stationary as well as non-stationary processes.
    delta = 0.7 # a tunning parameter set 0< delta < 1.
    d_est = ELW(score, floor(1+T^delta))$d # floor(1+T^delta) - bandwith parameter
    
    return(list(C = C, score = score, d_est = d_est,
                ncomp = ncomp, varprop = prop_eigen_value, 
                eigen_function = as.numeric(eigen_function),basis = basis))
}

