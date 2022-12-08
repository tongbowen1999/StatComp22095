#' @title Compute volatility
#' @description  compute the dynamic volatility of a stock by SMA/GMA/EWMA
#' @param data the close price of the stock
#' @param SMA return the SMA result or not
#' @param GMA return the GMA result or not
#' @param EWMA return the EWMA result or not
#' @return the dynamic volatility and some figures
#' @examples
#' \dontrun{
#' library(StatComp22095)
#' data(data)
#' volatility(data, EWMA = TRUE, GMA = FALSE, SMA = FALSE)
#' }
#' @importFrom graphics axis lines par
#' @importFrom stats dnorm qnorm
#' @importFrom utils head
#' @useDynLib StatComp22095
#' @export
volatility <- function(data, SMA = TRUE, GMA = TRUE, EWMA = TRUE){
  ##data
  shjc_data <- subset(data,select = c("date","close"))
  
  ##log(1+r)
  shjcnum <- dim(shjc_data)[1]
  SH <- data.frame(date=shjc_data$date[-1],rt= log(shjc_data$close[-1]/shjc_data$close[-shjcnum]))
  
  ##log plot
  par(mfrow = c(1,1))
  plot(SH$rt,type = "l",main = "log(1+r)",col = 3,ylim = c(-0.15,0.15),xlab ="",ylab = "")
  
  
  ## annual dynamic volatility
  #SMA
  N1 <- which(SH$date == "2016-01-04")
  N2 <- which(SH$date == "2021-12-31")
  Num <- N2 - N1 + 1
  SMA_SH <- rep(NA,Num)
  for(i in 1:Num){
    fInd <- N1 + i - 1
    mean_SH <- mean(SH$rt[(fInd - 4):(fInd - 1)])
    SMA_SH[i] <- sqrt(sum((SH$rt[(fInd - 4):(fInd - 1)] - mean_SH)^2/4))
  }
  if (SMA == TRUE){
    print('SMA:')
    print(head(SMA_SH))
  }
  #GMA
  GMA_SH <- rep(NA,Num)
  wt <- c(0.1,0.2,0.3,0.4)
  for(i in 1:Num){
    fInd <- N1 + i - 1
    mean_SH <- mean(SH$rt[(fInd - 4):(fInd - 1)])
    GMA_SH[i] <- sqrt(sum(wt*(SH$rt[(fInd - 4):(fInd - 1)] - mean_SH)^2))
  }
  if (GMA == TRUE){
    print('GMA:')
    print(head(GMA_SH))
  }
  
  #EWMA
  EWMA_SH <- rep(NA,Num);M=4
  lambda <- 0.94;
  vlambda <- c(1,lambda,lambda^2,lambda^3)
  wt <- sort(vlambda/sum(vlambda))
  for(i in 1:Num){
    fInd <- N1 + i - 1
    wtmean_SH <- mean(wt*SH$rt[(fInd-4):(fInd-1)])
    EWMA_SH[i] <- sqrt(sum(wt*(SH$rt[(fInd-4):(fInd-1)]-mean_SH)^2))
  }
  if (EWMA == TRUE){
    print('EWMA:')
    print(EWMA_SH)
  }
  
  ##plot
  Dynavol_SH <- data.frame(date = SH$date[N1:N2], SMA = SMA_SH,GMA = GMA_SH, EWMA = EWMA_SH)
  years <- c("2016","2017","2018","2019","2020","2021") 
  ind_m <- rep(NA,length(years))
  for(i in 1:length(years)){  
    ind_m[i] <- min(which(substr(Dynavol_SH$date,1,4) == years[i])) 
  }
  par(mfrow = c(1,1),pin = c(2,1.2))
  plot(Dynavol_SH$SMA,type = "l",main = "SMA",col = 1,ylim = c(0,0.1),ylab = "",xaxt = "n",xlab ="") 
  axis(side = 1,at = ind_m,labels = years) 
  
  par(mfrow = c(1,1),pin = c(2,1.2))
  plot(Dynavol_SH$GMA,type = "l",main = "GMA",col = 2,ylim = c(0,0.1),ylab = "",xaxt = "n",xlab ="")
  axis(side = 1,at = ind_m,labels = years)
  
  par(mfrow = c(1,1),pin = c(2,1.2))
  plot(Dynavol_SH$EWMA,type = "l",main = "EWMA",col = 3,ylim = c(0,0.1),ylab = "",xaxt = "n",xlab ="")
  axis(side = 1,at = ind_m,labels = years)
  
  print(summary(SH$rt))
}



#' @title Compute VAR and ES
#' @description use garch model to compute the VAR and ES of one stock
#' @param data the close price of the stock
#' @return VAR and ES for this stock and some figures
#' @examples
#' \dontrun{
#' library(StatComp22095)
#' data(data)
#' VAR_and_ES(data)
#' }
#' @importFrom fGarch garchFit
#' @export
VAR_and_ES <- function(data){
  zuhe_data = subset(data,select = c("date","close"))
  zuhenum = dim(zuhe_data)[1]
  ZUHE=data.frame(date=zuhe_data$date[-1],rt=log(zuhe_data$close[-1]/zuhe_data$close[-zuhenum]))#对数收益率
  
  ##fGarch to compute 1.4(day) VaR
  hT = 800;c = 0.99
  endN = which(ZUHE$date == "2020-12-31")
  staN = endN - hT
  
  
  fit_ZUHE = fGarch::garchFit(~ garch(1,1),data = ZUHE$rt[staN:endN])
  pred_ZUHE = fGarch::predict(fit_ZUHE,n.ahead = 1)##forward one day
  VaR_ZUHE = -(pred_ZUHE$meanForecast + pred_ZUHE$standardDeviation*qnorm(p = 1-c))
  
  ##VaR、ES(year)
  N3 = which(ZUHE$date == "2021-01-04")
  N4 = which(ZUHE$date == "2021-12-31")
  predict_Num = N4 - N3 + 1
  VaR_ZUHE = rep(NA,predict_Num)
  ES_ZUHE = rep(NA,predict_Num)
  for(j in 1:predict_Num){
    pDate = endN + j; pInd = (pDate - hT):(pDate - 1)
    fit_ZUHE = fGarch::garchFit(~ garch(1,1), data = ZUHE$rt[pInd])
    pred_ZUHE = fGarch::predict(fit_ZUHE,n.ahead = 1)
    VaR_ZUHE[j] = -(pred_ZUHE$meanForecast + pred_ZUHE$standardDeviation*qnorm(p = 1-c))
    ES_ZUHE[j]=-(pred_ZUHE$meanForecast+pred_ZUHE$standardDeviation*(-dnorm(qnorm(p = 1-c))/(1-c)))
  }
  
  Risk_ZUHE = data.frame(date = ZUHE$date[N3:N4],VaR = VaR_ZUHE,ES = ES_ZUHE)
  mon = c("2021-01","2021-03","2021-05","2021-07","2021-09","2021-12")
  ind_m1 = rep(NA,length(mon))
  for(i in 1:length(mon)){   
    ind_m1[i] = min(which(substr(Risk_ZUHE$date,1,7) == mon[i])) 
  }
  print(Risk_ZUHE)
  
  par(mfrow = c(1,1),pin = c(2,1.2))
  plot(-Risk_ZUHE$VaR,type = "l",main = "VaR",col = 1,ylim = c(-0.15,0.15),ylab = "",xaxt = "n",xlab = "")
  lines(ZUHE$rt[N3:N4],col = 4)
  axis(side = 1,at = ind_m1,labels = mon)
  
  par(mfrow = c(1,1),pin = c(2,1.2))
  plot(-Risk_ZUHE$ES,type = "l",main = "ES",col = 1,ylim = c(-0.15,0.15),ylab = "",xaxt = "n",xlab = "")
  lines(ZUHE$rt[N3:N4],col = 2)
  axis(side = 1,at = ind_m1,labels = mon)
}

