

###input parameters
n2pwr.maxLR<- function(entry   = 1
                       ,fup      = 1
                       ,CtrlHaz
                       ,hazR
                       ,transP1
                       ,transP0
                       ,Wlist
                       ,entry_pdf0=function(x){(1/entry)*(x>=0&x<=entry)}
                       ,entry_pdf1=entry_pdf0
                       ,eventN
                       ,totalN
                       ,ratio    = 1
                       ,alpha    = 0.05
                       ,alternative=c("two.sided","less","greater")
                       ,k        = 100
                       ,criteria = 100
){
  if (missing(eventN)&missing(totalN)){
    stop("At least one of eventN/totalN must be provided")
  }
  tot_time <- entry+fup
  num <- k*tot_time
  # create the subintervals
  x <- seq(0,tot_time,by=1/k)
  ctrlRate <- CtrlHaz(x)
  haz_val <- hazR(x)*ctrlRate
  haz_point <- x*k
  ## load the transition matrix
  load <- trans.mat(num=num,x=x,ctrlRate=ctrlRate,haz_val=haz_val,
                    haz_point=haz_point,ratio=ratio,
                    transP1=transP1,transP0=transP0,k=k,
                    fup=fup,entry=entry,entry_pdf0=entry_pdf0,
                    entry_pdf1=entry_pdf1)

  pdat <- load$pdat
  eprob <- stats::weighted.mean(c(pdat$C_E[num],pdat$E_E[num]),w=c(1,ratio))
  if (missing(eventN)){ eventN <- round(totalN*eprob)}
  if (missing(totalN)){ totalN <- round(eventN/eprob)}
  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(pdat),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)

  for (j in 1:wn){
    W[,j] <- Wlist[[j]](pdat$S)
  }
  dnum <- eventN
  for (k1 in 1:wn){
    for (k2 in 1:wn){
      Vmat[k1,k2] <-  dnum*t(W[,k1]*W[,k2]) %*%(pdat$rho*pdat$eta)
    }
  }
  rho_est <- stats::cov2cor(Vmat)
  mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))/sqrt(diag(Vmat))
  if (alternative=="two.sided"){
    ftwo <- function(i){
      crit <- mvtnorm::qmvnorm(1-alpha,tail="both.tails",
                               mean=rep(0,wn),sigma = rho_est)$quantile
      power <- 1-mvtnorm::pmvnorm(-crit,crit,mean=mu,sigma = rho_est)
      return(c(crit,power))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    power <- ftwod[2]
  }else if (alternative=="less"){ # l1 <l0
    ftwo <- function(i){
      crit <- mvtnorm::qmvnorm(alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
      power <- 1-mvtnorm::pmvnorm(crit,Inf,mean=mu,sigma = rho_est)[1]
      return(c(crit,power))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    crit <- ftwod[1];  power <- ftwod[2]
  }else if (alternative=="greater"){
    ftwo <- function(i){
      crit <- qmvnorm(alpha,tail="upper.tail",mean=rep(0,wn),sigma = rho_est)$quantile
      power <- 1-pmvnorm(-Inf,crit,mean=mu,sigma = rho_est)[1]
      return(c(crit,power))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    crit <- ftwod[1];  power <- ftwod[2]
  }


  listall <- list( power = as.numeric(power)
                   ,eventN  = eventN
                   ,totalN = totalN
                   ,prob_event =eprob
                   ,L_trans = load$L_trans
                   ,pdat = pdat
                   ,studytime=c(entry,fup)
                   ,RandomizationRatio=ratio

  )

  return(listall)


}
