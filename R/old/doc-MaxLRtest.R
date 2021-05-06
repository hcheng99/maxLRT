#' @title MaxLRtest
#' @description perform the maximum weighted logrank test
#' @param dat the dataset containing time, event and group in the 1st-3rd column
#' @param Wlist list containing the weight functions (see example below)
#' @param base the base for weight function, the default is pooled KM estimate.
#'  options: c("Combined", "KM", "N")
#' @param alpha significance level, Default: 0.05
#' @param side two-sided or one-sided(less or greater) test, Default: c("two.sided")
#' @return stat: the maximum statistic
#' @return stat_mat: statistics from each weight function and covariance matrix
#' @return critV: the critical value
#' @return details: the dataset containing information at each event time
#' point
#' @return  pvalue: p-value
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname MaxLRtest
#' @export
#' @import dplyr
#' @import survival
MaxLRtest <- function(dat
                        ,Wlist
                        ,base=c("Combined","KM","N")
                        ,alpha=0.05
                        ,side=c("two.sided")

){
  ## transform the input data to the matrix ready for analysis
  datM <- load.mat(dat)
  ## only keep the timepoints with event
  datM <- datM[datM$event>0,]
  if (base=="Combined"){

    tnew0 <- c(1,
              cumprod(chkV(datM$n.risk.x,1-datM$n.event.x/datM$n.risk.x))*table(dat$V3)[2]+
      cumprod(chkV(datM$n.risk.y,1-datM$n.event.y/datM$n.risk.y))*table(dat$V3)[1])
    tnew <- 1-tnew0[1:nrow(datM)]

  }else if (base=="KM"){
    # based on the pooled survival
    s_fit<- survfit(Surv(dat[,1], dat[,2])~1 , data = dat)
    f_s <- stepfun(s_fit$time,y=c(0,1-s_fit$surv),right = TRUE)
    tnew <- f_s(datM$time)
  }else if (base=="N"){
    # based on time/person at risk
    # tnew <- datM$time/max(datM$time)
    tnew <- 1-datM$risk/nrow(dat)
  }
  ## number of weight functions
  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(datM),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)
  for (j in 1:wn){
    W[,j] <- Wlist[[j]](tnew)
  }
  for (k1 in 1:wn){
    for (k2 in 1:wn){
      Vmat[k1,k2] <-  t(W[,k1]*W[,k2]) %*%datM$V
    }
  }
  Zstat <- apply(W,2,function(x) {CalZ(x,data=datM)})
  ## correlation matrix
  rho_est <- cov2cor(Vmat)
  if (side=="two.sided"){
    crit <- qmvnorm(1-alpha,tail="both.tails",mean=rep(0,wn),sigma = rho_est)$quantile
    statistic <- max(abs(Zstat))
    p.value <- 1-pmvnorm(-1*statistic,statistic,mean=rep(0,wn),sigma = rho_est)[1]
  }else if (side=="greater"){
    crit <- qmvnorm(alpha,tail="upper.tail",mean=rep(0,wn),sigma = rho_est)$quantile
    statistic <- max(abs(Zstat))*sign(Zstat[1])
    p.value <- 1-pmvnorm(-Inf,statistic,mean=rep(0,wn),sigma = rho_est)[1]
  }else if (side=="less"){
    crit <- qmvnorm(alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
    statistic <- max(abs(Zstat))*sign(Zstat[1])
    p.value <- 1-pmvnorm(statistic,Inf,mean=rep(0,wn),sigma = rho_est)[1]

  }

  res <- statistic>=crit
  stat.mat <- data.frame(Zstat,rho_est)
  names(stat.mat) <- c("Zstat",paste0("W",1:wn))
  Rdat <- data.frame(datM[,c("time","event","risk","n.event.x","exp.x","n.risk.x","V")],W)
  names(Rdat)[-c(1:7)] <- paste0("W",1:wn)
  list <-list(stat=statistic
              ,stat.mat=stat.mat
              ,critV=crit
              ,details=Rdat
              ,pvalue=p.value)
  class(list) <- 'maxLR'

  return(list)

}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname plot.maxLR
#' @export
plot.maxLR <- function(x) {
    datM <- x$details
    xtime <- datM$time
    wn <- ncol(datM)-7
    plot(xtime,datM$n.event.x,cex=1,col=1,xlab="time",
         ylab="weight",pch=1,ylim=c(-1,1))
    for (i in 1:wn){
      points(xtime,datM[,7+i],cex=0.3,col=1+i,pch=1+i)
    }

    legend("bottomright",legend=c("orginal data",paste("weight",1:wn)),
           col=1:(wn+1),pch = 1:(wn+1),cex=0.8)
    mtext(paste(paste0("Statistic",1:wn,sep=":"),round(x$stat.mat[,1],3),
                collapse = "  "),side=3,cex=0.7)
}
