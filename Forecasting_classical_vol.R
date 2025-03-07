library('sde')
library('xtable')
library('mvtnorm')
library('xtable')
library('tseries')
library('tictoc')
library('optimx')
library('compiler')
library('parallel')
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
library('highfrequency')
library('HARModel')
library('reshape2')
library('tidyr')
library('zoo')
library('scoringRules')
library('forecast')
library('MCS')
library('garchx')
library('arfima')

setwd('your working directory/')
euspot = read.csv('euspot_24.csv',sep=",",header=FALSE)
ukspot = read.csv('ukspot_24.csv',sep=",",header=FALSE)
jpspot = read.csv('jpspot_24.csv',sep=",",header=FALSE)

eubid = read.csv('eubid_24.csv',sep=",",header=FALSE)
euask = read.csv('euask_24.csv',sep=",",header=FALSE)
ukbid = read.csv('ukbid_24.csv',sep=",",header=FALSE)
ukask = read.csv('ukask_24.csv',sep=",",header=FALSE)
jpbid = read.csv('jpbid_24.csv',sep=",",header=FALSE)
jpask = read.csv('jpask_24.csv',sep=",",header=FALSE)

realvol_overnight<-function(intra_ret,scalar_onr){
  point_grid=nrow(intra_ret)
  N=ncol(intra_ret)
  real_onvol=c(rep(0,N))
  for (i in 1:N){
    real_onvol[i]=sqrt(sum(intra_ret[,i]^2)+scalar_onr[i]^2)
  }
  return(real_onvol)
}

# daily return
day.return <- function(yd, ip){
  grid_point=nrow(yd)
  N=ncol(yd)
  day_re=c(rep(0,N))
    
  for(i in c(2:N)){
    day_re[i]=log(yd[ip,i])-log(yd[ip,i-1])
  }
  return( day_re*100)
}

# overnight return
on.return <- function(yd, ip){
  grid_point=nrow(yd)
  N=ncol(yd)
  day_re=c(rep(0,N))
    
  for(i in c(2:N)){
    day_re[i]=log(yd[1,i])-log(yd[ip,i-1])
  }
  return( day_re*100)
}

# intraday log return
idr.return <- function(yd, ip){
  # grid_point=nrow(yd)
  grid_point = ip
  N=ncol(yd)
  idr_re=matrix(0,grid_point-1,N)
  
  for(j in c(2:grid_point)){
    for(i in c(1:N)){
      idr_re[j-1,i]=log(yd[j,i])-log(yd[j-1,i])
    }
  }
  return( idr_re*100 )
}

# daily bid ask spread
bid.ask <- function(bid, ask, ip){
   grid_point=nrow(bid) 
   N = ncol(bid)
   spread = c(rep(NA, N))
   
   for (i in c(1:N)){
     ask_d = ask[ip,i]
     bid_d = bid[ip,i]
     mid_d = (ask_d+bid_d)/2
     spread[i] = (ask_d-bid_d)/mid_d
   }
   return(spread)
}


N = ncol(euspot)
ip = 193 ##
euspot = euspot[,2:N]
ukspot = ukspot[,2:N]
jpspot = jpspot[,2:N]

eubid = eubid[,2:N]
ukbid = ukbid[,2:N]
jpbid = jpbid[,2:N]

euask = euask[,2:N]
ukask = ukask[,2:N]
jpask = jpask[,2:N]

eu_idr = idr.return(euspot, ip)
uk_idr = idr.return(ukspot, ip)
jp_idr = idr.return(jpspot, ip)

# demean
eu_idr=eu_idr-rowMeans(eu_idr)
uk_idr=uk_idr-rowMeans(uk_idr)
jp_idr=jp_idr-rowMeans(jp_idr)

dbidask_eu = bid.ask(eubid, euask, ip)
dbidask_uk = bid.ask(ukbid, ukask, ip)
dbidask_jp = bid.ask(jpbid, jpask, ip)

eu_clre = day.return(euspot, ip)
uk_clre = day.return(ukspot, ip)
jp_clre = day.return(jpspot, ip)

# demean
eu_clre = eu_clre-mean(eu_clre)
uk_clre = uk_clre-mean(uk_clre)
jp_clre = jp_clre-mean(jp_clre)

eu_onre = on.return(euspot, ip)
uk_onre = on.return(ukspot, ip)
jp_onre = on.return(jpspot, ip)

## GARCH

forecast_garch<-function(day_re){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){## 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    
    G1=garchFit(~garch(1,1),sub_yd,trace=FALSE)
 
    preds = c(predict(G1,n.ahead=out_winv[i] )$standardDeviation)

    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = preds
  }
  
  return( vol_pred )
}


eu_garch = forecast_garch(eu_clre)
uk_garch = forecast_garch(uk_clre)
jp_garch = forecast_garch(jp_clre)


## GJR

forecast_gjr<-function(day_re){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){# 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]

    spec <- ugarchspec(
      variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)
    )
    G1 <- ugarchfit(spec = spec, data = sub_yd, solver = "hybrid", trace = FALSE)
    preds <- sigma(ugarchforecast(G1, n.ahead = out_winv[i]))
    vol_pred[(insample_ind[i, 2] + 1 - win_N):(insample_ind[i, 2] + out_winv[i] - win_N)] <- as.numeric(preds)
  }
  
  return( vol_pred )
}


eu_gjr = forecast_gjr(eu_clre)
uk_gjr = forecast_gjr(uk_clre)
jp_gjr = forecast_gjr(jp_clre)


## GARCH-X


forecast_garchx<-function(day_re, dbidask){
  N=length(day_re) #1274  
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){# 
    #od = outsample_ind[i+1]-1
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    sub_x = abs(dbidask[insample_ind[i,1]:insample_ind[i,2]])

    specx = ugarchspec(mean.model=list(armaOrder=c(0,0)),
                  variance.model = list(external.regressors = as.matrix(sub_x) ) )
    G1 = ugarchfit(spec=specx, data=sub_yd,  out.sample=out_winv[i] )
 
    preds = sigma(ugarchforecast(G1, n.ahead=out_winv[i], externalforecasts = as.matrix(sub_x) ) )#n.roll = out_winv[i], 
 
    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = c(preds)
  }
  
  return( vol_pred )
}

eu_x = forecast_garchx(eu_clre, dbidask_eu)
uk_x = forecast_garchx(uk_clre, dbidask_uk)
jp_x = forecast_garchx(jp_clre, dbidask_jp)


## FIGARCH
forecast_figarch<-function(day_re){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)
  
  specFigarch = ugarchspec(mean.model=list(armaOrder=c(0,0)),
                  variance.model = list(model = "fiGARCH",submodel="GARCH", garchOrder = c(1,1)),
                  distribution="norm")

  for (i in 1:mdl_rep){ 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    
    G1 = ugarchfit(spec=specFigarch, data=sub_yd, fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all',trunclag = 500))#, solver.control=list(trace = 1))

    preds = c(sigma(ugarchforecast(G1, n.ahead=out_winv[i] ) ))
 
    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = preds
  }
  
  return( vol_pred )
}

eu_fi = forecast_figarch(eu_clre)
uk_fi = forecast_figarch(uk_clre) 
jp_fi = forecast_figarch(jp_clre)


## HAR
forecast_har<-function(day_re, idr, onr){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){ 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    sub_idr = idr[,insample_ind[i,1]:insample_ind[i,2]]
    sub_onr = onr[insample_ind[i,1]:insample_ind[i,2]]

    datesim=seq(as.Date("2000/1/1"), by = "day", length.out = length(sub_yd))
    real_feed=realvol_overnight(sub_idr,sub_onr)
    real_feed=as.xts(real_feed,datesim)

    G1=HARModel::HARForecast(real_feed, periods = c(1,5,22), nRoll=out_winv[i], nAhead=out_winv[i], type = "HAR")
    preds = c(getForc(G1))
 

    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = preds
  }
  
  return( vol_pred )
}

eu_har = forecast_har(eu_onre,eu_idr,eu_onre)
uk_har = forecast_har(uk_onre,uk_idr,uk_onre)
jp_har = forecast_har(jp_onre,jp_idr,jp_onre)


## HAR-X
forecast_harx<-function(day_re, idr, dbidask, onr){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){ 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    sub_idr = idr[,insample_ind[i,1]:insample_ind[i,2]]
    sub_x = abs(dbidask[insample_ind[i,1]:insample_ind[i,2]])
    sub_onr = onr[insample_ind[i,1]:insample_ind[i,2]]

    datesim=seq(as.Date("2000/1/1"), by = "day", length.out = length(sub_yd))
    real_feed1=realvol_overnight(sub_idr,sub_onr)
    real_feed=as.xts(real_feed1,datesim)
    x_feed = as.xts(sub_x,datesim)
    
    real_feedc = real_feed
    preds = c()
    for(j in 1:out_winv[i]){
      G1=highfrequency::HARmodel(data = real_feedc, periods = c(1,5,22), RVest = c("rCov"),type = "HAR", inputType = "RM", externalRegressor = x_feed)
      preds =  c(preds, predict(G1))

      real_feedc = as.xts(c(real_feed1,preds)[(i+1):(i+win_N)],datesim)
    }
 
    
    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = preds
  }
  
  return( vol_pred )
}

eu_harx = forecast_harx(eu_clre,eu_idr,dbidask_eu,eu_onre)
uk_harx = forecast_harx(uk_clre,uk_idr,dbidask_uk,uk_onre)
jp_harx = forecast_harx(jp_clre,jp_idr,dbidask_jp,jp_onre)



## HAR-ARFIMA
forecast_rvarma<-function(day_re, idr, onr){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){ 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    sub_idr = idr[,insample_ind[i,1]:insample_ind[i,2]]
    sub_onr = onr[insample_ind[i,1]:insample_ind[i,2]]

    real_feed=realvol_overnight(sub_idr,sub_onr)

    G1=arfima::arfima(real_feed, order=c(1,0,1))
    preds = predict(G1, n.ahead = out_winv[i], predint = 0.95)[[1]]$Forecast
 

    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = preds
  }
  
  return( vol_pred )
}

eu_rvarma = forecast_rvarma(eu_onre,eu_idr,eu_onre)
uk_rvarma = forecast_rvarma(uk_onre,uk_idr,uk_onre)
jp_rvarma = forecast_rvarma(jp_onre,jp_idr,jp_onre)




## RGARCH
forecast_rgarch<-function(day_re, idr, onr){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){ 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    sub_idr = idr[,insample_ind[i,1]:insample_ind[i,2]]
    sub_onr = onr[insample_ind[i,1]:insample_ind[i,2]]

    real_feed=realvol_overnight(sub_idr,sub_onr)
    sub_yd = as.xts(sub_yd,seq(as.Date("2014/3/1"), by = "day", length.out = length(sub_yd)))
    real_feed = as.xts(real_feed,seq(as.Date("2014/3/1"), by = "day", length.out = length(real_feed)))


    spec = ugarchspec(mean.model = list(armaOrder = c(1, 1), include.mean = FALSE), variance.model = list(model = 'realGARCH', garchOrder = c(1, 1)))
    setbounds(spec)<-list(alpha2=c(-1,1))
    fit = ugarchfit(spec, sub_yd, solver = 'hybrid', solver.control = list(trace = TRUE), realizedVol = real_feed)#'hybrid'
    rgarch_pred = ugarchforecast(fit, n.ahead=out_winv[i])
    sigma_pred=sigma(rgarch_pred)
    
    if(sigma_pred[1]>2){print(i)}

    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = sigma_pred
  }
  
  return( vol_pred )
}

eu_rgarch = forecast_rgarch(eu_clre,eu_idr,eu_onre)
uk_rgarch = forecast_rgarch(uk_clre,uk_idr,uk_onre)
jp_rgarch = forecast_rgarch(jp_clre,jp_idr,jp_onre)



# GARCH - MIDAS
library('rumidas')

eumidas = read.csv('euspot_midas.csv',sep=",",header=FALSE)
ukmidas = read.csv('ukspot_midas.csv',sep=",",header=FALSE)
jpmidas = read.csv('jpspot_midas.csv',sep=",",header=FALSE)
ip = 193

realvol<-function(data){
  point_grid=nrow(data); N=ncol(data)
  real_onvol=c(rep(0,N))
  for (i in 1:N){
    real_onvol[i]=sqrt(sum(data[,i]^2))
  }
  return(real_onvol)
}

real_eu=realvol(eumidas)
real_uk=realvol(ukmidas)
real_jp=realvol(jpmidas)

eu_mrv = colMeans(matrix(real_eu, 16)) # rollmean  rollsum
uk_mrv = colMeans(matrix(real_uk, 16))
jp_mrv = colMeans(matrix(real_jp, 16))

eu_mrv = as.xts(eu_mrv,seq(as.Date("2003/1/1"), by = "months", length.out = length(eu_mrv)))
uk_mrv = as.xts(uk_mrv,seq(as.Date("2003/1/1"), by = "months", length.out = length(uk_mrv)))
jp_mrv = as.xts(jp_mrv,seq(as.Date("2003/1/1"), by = "months", length.out = length(jp_mrv)))


forecast_midas<-function(day_re, mx){
  N=length(day_re) #1274
  win_N = 600
  out_N = N-win_N
  out_win = 4
  
  insample_ind = as.matrix(cbind(c(seq(1,N,out_win)),c(seq(win_N,N,out_win)))[1:169,])
  outsample_ind = c(seq(win_N,N,out_win),N)
  mdl_rep = nrow(insample_ind)

  # containers for predictors
  vol_pred=matrix(NA,out_N,1)
  out_winv = diff(outsample_ind)

  for (i in 1:mdl_rep){ 
    sub_yd = day_re[insample_ind[i,1]:insample_ind[i,2]]
    trun_y = sub_yd[421:600]

    trun_y = as.xts(trun_y,seq(as.Date("2014/3/1"), by = "day", length.out = length(trun_y)))
    midas_rv = mv_into_mat(trun_y,mx,K=12,type="monthly")

    spec = ugmfit(model="GM", skew="NO", distribution="norm", trun_y, midas_rv, K=12)
    midas_pred = multi_step_ahead_pred(spec, out_winv[i])

    vol_pred[ (insample_ind[i,2]+1 - win_N):(insample_ind[i,2]+out_winv[i] - win_N) ] = c(midas_pred)
  }  
  return( vol_pred )
}

eu_midas = forecast_midas(eu_clre,eu_mrv)
uk_midas = forecast_midas(uk_clre,uk_mrv)
jp_midas = forecast_midas(jp_clre,jp_mrv*10)



save(eu_garch,eu_gjr, eu_x,eu_fi,eu_har,eu_harx,eu_rvarma, eu_rgarch, eu_midas, file = "peer_eu.RData" )
save(uk_garch,uk_gjr, uk_x,uk_fi,uk_har,uk_harx,uk_rvarma, uk_rgarch, uk_midas, file = "peer_uk.RData" )
save(jp_garch,jp_gjr, jp_x,jp_fi,jp_har,jp_harx,jp_rvarma, jp_rgarch, jp_midas, file = "peer_jp.RData" )

# setwd('your directory')
# save(eu_rgarch, file = "rgarch_eu.RData" )
# save(uk_rgarch, file = "rgarch_uk.RData" )
# save(jp_rgarch, file = "rgarch_jp.RData" )