# functions to load
int_approx=function(mat){
    temp_n=NROW(mat)
    return((1/temp_n)*sum(mat))
} 

positive_c <- function(x) { x[x<0] <- 0 ; x}
    
error_sim <- function(grid_point, samplesize){
  ti=1:grid_point/grid_point
  comat=matrix(NA,grid_point,grid_point)
  for (i in 1:grid_point){
    comat[i,]=exp(1)^(-ti[i]/2-ti/2)*pmin(exp(1)^(ti[i]),exp(1)^(ti))
  }
  epsilon=mvrnorm(n = samplesize, mu = c(rep(0,grid_point)), Sigma = comat, empirical = FALSE)
  return(t(epsilon))
}

####################
## fgarch model
####################
# input - y: objective data, with rows present the grid points, and columns present the sample observations
#         basis: the chosen bases for dimention reduction
fgarch_est<-function(y,basis){
  basis=as.matrix(basis)
  M=ncol(basis)
  point_grid=nrow(y)
  T=ncol(y)

  function_to_minimize=function(x,data1){
   # data ranges from y1 to yT  ; 
   # basis (a [point_grid X M] matrix)
   # x the parameters need to be optimised
   undat=as.matrix(data1[[1]])
   basis=as.matrix(data1[[2]])
   
   T=ncol(undat)
   point_grid=nrow(undat)
   M=ncol(basis)
   
   # standardize data
   undat=scale(undat,center=TRUE,scale=FALSE)
   data=matrix(0,point_grid,T)
   for (i in 1:point_grid){
     data[i,]=(undat[i,])/sd(undat[i,])
   }
   
   ### set scalar coefficients
   cons=x[1:M]
   
   alpha=matrix(as.vector(x[(M+1):(M+M^2)]),M,M)
   beta=matrix(as.vector(x[(M+M^2+1):(M+2*M^2)]),M,M)

   int_approx=function(mat){
        temp_n=NROW(mat)
        return((1/temp_n)*sum(mat))}    
   
   ### set kernel coefficients
   times = 1:point_grid/point_grid
   delta=0
   for (k in 1:M){
       delta=delta+cons[k]*basis[,k]
   }
    
   alpha_Op=matrix(0,NROW(times),NROW(times))
      for(i in 1:NROW(times)){
          for(j in 1:NROW(times)){
          temp_alpha_m=0
          for(t in 1:M){
          for(s in 1:M){
              temp_alpha_m=temp_alpha_m+alpha[t,s]*basis[i,t]*basis[j,s]
               }
          }
      alpha_Op[i,j]=temp_alpha_m
      }}

   beta_Op=matrix(0,NROW(times),NROW(times))
      for(i in 1:NROW(times)){
          for(j in 1:NROW(times)){
          temp_beta_m=0
          for(t in 1:M){
          for(s in 1:M){    
              temp_beta_m=temp_beta_m+beta[t,s]*basis[i,t]*basis[j,s]
               }
          }
        beta_Op[i,j]=temp_beta_m
      }}
  
   ### set initial values of y and sigma
   est_sigma2=matrix(NA, point_grid, T) 
   est_sigma2[,1]=matrix(delta)
                    
   ### function to minimise   
   elle=as.vector(rep(0,T-1))

   for (i in 2:T){
       for (t in 1:point_grid){
                   arch_part=alpha_Op[t,]*(data[,i-1]^2)
                   garch_part=beta_Op[t,]*est_sigma2[,i-1]
           est_sigma2[t,i]=delta[t]+int_approx(arch_part)+int_approx(garch_part)
       }

       elle_t=0
       for (k in 1:M){
       elle_t=elle_t+int_approx(data[,i]^2*basis[,k])/int_approx(est_sigma2[,i]*basis[,k])+log(int_approx(est_sigma2[,i]*basis[,k]))}
     
       elle[i-1]=elle_t}

       elle=na.omit(elle)
       m_elle=mean(elle)
       if(is.nan(m_elle)==TRUE){
         m_elle=99
       }
 
      return(m_elle)
  }

  stav=c(rep(0.1,M),rep(0.1,M^2), rep(0.3,M^2))######Caution !!
  lower1=c(rep(0.01,M),rep(0.01,M^2), rep(0.01,M^2))
  upper1=c(rep(1,M),rep(0.95,2*M^2)) ######Caution !! 

  y1=list(y,basis)
  para=optimx(par=as.vector(stav),fn=function_to_minimize,data1=y1,lower=lower1,upper=upper1,itnmax=1000,method="L-BFGS-B")#,gR=NULL,lower=lower1,upper=upper1,itnmax=1000,method="L-BFGS-B")bobyqa
  para=as.vector(as.numeric(para[1:(M+2*M^2)]))

  d_hat=para[1:M]
  alpha_hat=matrix(para[(M+1):(M+M^2)],M,M)
  beta_hat=matrix(para[(M+M^2+1):(M+2*M^2)],M,M)

    times = 1:point_grid/point_grid
    delta_hat=0
    for (i in 1:M){
      delta_hat=delta_hat+d_hat[i]*basis[,i]
    }
    delta_hat=positive_c(delta_hat)

    alpha_Op_hat=matrix(0,NROW(times),NROW(times))
      for(i in 1:NROW(times)){
          for(j in 1:NROW(times)){
          temp_alpha_m=0
          for(t in 1:M){
          for(s in 1:M){
              temp_alpha_m=temp_alpha_m+alpha_hat[t,s]*basis[i,t]*basis[j,s]
               }
          }
      alpha_Op_hat[i,j]=temp_alpha_m
      }}
    alpha_Op_hat=positive_c(alpha_Op_hat)

    beta_Op_hat=matrix(0,NROW(times),NROW(times))
      for(i in 1:NROW(times)){
            for(j in 1:NROW(times)){
            temp_beta_m=0
            for(t in 1:M){
            for(s in 1:M){    
              temp_beta_m=temp_beta_m+beta_hat[t,s]*basis[i,t]*basis[j,s]
               }
           }
        beta_Op_hat[i,j]=temp_beta_m
      }}
    beta_Op_hat=positive_c(beta_Op_hat)


    sigma2_fit=matrix(1,point_grid,T+1)
    sigma2_fit[,1]=delta_hat
    for(j in 2:(T+1)){
        #first fill in sigma2 column:
        for(i in 1:point_grid)
        { 
            fit_alpha_op = alpha_Op_hat[i,] * (y[,(j-1)])^2
            fit_beta_op = beta_Op_hat[i,] * (sigma2_fit[,j-1])
            
            sigma2_fit[i,j] = delta_hat[i] + int_approx(fit_alpha_op) + int_approx(fit_beta_op)
        }     }
    error_fit=y/sqrt(abs(sigma2_fit[,1:T]))
    error_fit[,1]=0

    y_fit=sqrt(abs(sigma2_fit[,1:T]))*error_sim(point_grid,T)

    return(list(sigma2_fit,delta_hat,alpha_Op_hat,beta_Op_hat,para,error_fit,y_fit)) #######error_fit for bootstrapping
}
# output - sigma2_fit: predicted conditional volatility
#          delta_hat: estimated unconditional volatility 
#          alpha_Op_hat: estimated ARCH coefficient operator
#          beta_Op_hat: estimated GARCH coefficient operator
#          para: estimated scalar parameters
#          error_fit: estimated innovations with 
#          y_fit: fitted/predicted y with the assumption that innovation following an OU process.

###################
### FGARCH-X model
###################

fgarchx_est<-function(y,xc,basis){
  basis=as.matrix(basis)
  M=ncol(basis)
  point_grid=nrow(y)
  T=ncol(y)

  function_to_minimize=function(x,data1){
   # x the parameters need to be optimised
   data=as.matrix(data1[[1]])
   x_data=as.matrix(data1[[2]])
   basis=as.matrix(data1[[3]])
   q=3

   T=ncol(data)
   point_grid=nrow(data)
   basis=as.matrix(basis)
   M=ncol(basis)
   
   ### set scalar coefficients
   cons=x[1:M]

   alpha=list()
   for (h in 1:q){
   alpha_c=matrix(as.vector(x[(M+(h-1)*M^2+1):(M+h*M^2)]),M,M)
   alpha[[h]]<-alpha_c}

   int_approx=function(mat){
        temp_n=NROW(mat)
        return((1/temp_n)*sum(mat))}    
   positive_c <- function(x) { x[x<0] <- 0 ; x}
   ### set kernel coefficients
   times = 1:point_grid/point_grid
   delta=0
   for (k in 1:M){
       delta=delta+cons[k]*basis[,k]
   }
   
   alpha_Op=list()
   for (h in 1:q){   
       alpha_M=matrix(0,NROW(times),NROW(times))
       A_mat=as.matrix(alpha[[h]])
       for(i in 1:NROW(times)){
          for(j in 1:NROW(times)){
          temp_alpha_m=0
          for(t in 1:M){
          for(s in 1:M){
              temp_alpha_m=temp_alpha_m+A_mat[t,s]*basis[i,t]*basis[j,s]
               }
          }
      alpha_M[i,j]=temp_alpha_m
      }}
      alpha_Op[[h]]=positive_c(alpha_M)
   }
  
   ### set initial values of y and sigma
   est_sigma2=matrix(0, point_grid, T)   
   est_sigma2[,1]=delta

   ### function to minimise   
   elle=as.vector(rep(0,T-1))
   for (i in 2:T){
       for (t in 1:point_grid){
                   arch_part = garch_part = x_part = 0

                   arch_part=alpha_Op[[1]][t,] * data[,i-1]^2
                   garch_part=alpha_Op[[2]][t,] * est_sigma2[,i-1]
                   x_part=alpha_Op[[3]][t,] * x_data[,i-1]
                   
           est_sigma2[t,i]=delta[t]+int_approx(arch_part)+int_approx(garch_part)+int_approx(x_part)
       }
       elle_t=0
       for (k in 1:M){
       elle_t=elle_t+(int_approx(data[,i]^2*basis[,k])/int_approx(est_sigma2[,i]*basis[,k]))+log(int_approx(est_sigma2[,i]*basis[,k]))}
     
      elle[i-1]=elle_t}
      elle=na.omit(elle)
      m_elle=mean(elle)
       if(is.nan(m_elle)==TRUE){
         m_elle=99
       }
      return(m_elle)
  }

  q=3
  stav=c(rep(0.1,M),rep(0.1,M^2),rep(0.5,M^2),rep(0.05,M^2))######Caution !!
  lower1=c(rep(0.01,M), rep(0.01,q*M^2))
  upper1=c(rep(1,M),rep(0.95,(q*M^2))) ######Caution !! 

  y1=list(y,xc,basis)

  para=optimx(par=as.vector(stav),fn=function_to_minimize,data1=y1,lower=lower1,upper=upper1,method="L-BFGS-B")#,gR=NULL,lower=lower1,upper=upper1,itnmax=1000,method="L-BFGS-B")
  para=as.vector(as.numeric(para[1:(M+q*(M^2))]))

  d_hat=para[1:M]
  alpha_hat=list()
  for (h in 1:q){
    alpha_hat[[h]]=matrix(para[(M+(h-1)*M^2+1):(M+h*M^2)],M,M)
  }

  times = 1:point_grid/point_grid
  delta_hat=0
  for (i in 1:M){
    delta_hat=delta_hat+d_hat[i]*basis[,i]
  }

  alpha_Op_hat=list()
  for (h in 1:q){
       alpha_M=matrix(0,NROW(times),NROW(times))
       A_mat=as.matrix(alpha_hat[[h]])
       for(i in 1:NROW(times)){
          for(j in 1:NROW(times)){
          temp_alpha_m=0
          for(t in 1:M){
          for(s in 1:M){
              temp_alpha_m=temp_alpha_m+A_mat[t,s]*basis[i,t]*basis[j,s]
               }
          }
      alpha_M[i,j]=temp_alpha_m
      }}
      alpha_Op_hat[[h]]=positive_c(alpha_M)
  }

  sigma2_fit=matrix(NA,point_grid,T+1)
  sigma2_fit[,1]=delta_hat

  for(i in 2:(T+1)){
        #first fill in sigma2 column:
        for(j in 1:point_grid)
        {   fit_arch_op = fit_garch_op = fit_x_op = 0
            fit_arch_op = alpha_Op_hat[[1]][j,] * y[,i-1]^2
            fit_garch_op = alpha_Op_hat[[2]][j,] * (sigma2_fit[,i-1])
            fit_x_op = alpha_Op_hat[[3]][j,] * xc[,i-1]
            
            sigma2_fit[j,i] = delta_hat[j] + int_approx(fit_arch_op) + int_approx(fit_garch_op) + int_approx(fit_x_op)
        }}
  error_fit=y/sqrt(abs(sigma2_fit[,1:T]))
  error_fit[,1]=0
  y_fit=sqrt(abs(sigma2_fit[,1:T]))*error_sim(point_grid,T)
  return(list(sigma2_fit,delta_hat,alpha_Op_hat,para,error_fit,y_fit))
}
