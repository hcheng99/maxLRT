###input parameters
#' @title Power Calculation with Maximum Weighted Logrank Test
#' @description  \code{n2pwr.maxLR} calcualtes the power given either the
#' number of events or number of subjects
#' @param entry a numeric value indicating the enrollment time, Default: 1
#' @param fup a numeric value indicating the minimum follow-up time for subjects.
#'  , Default: 1
#' @param CtrlHaz a function,  specifying the hazard function for control group.
#' @param hazR a function, specifying the hazard ratio function between
#' treatment and control group
#' @param transP1 a numeric vector of length 2, consisting of the transition
#' probability from
#' receiving treatment to drop-out (drop-out rate) and
#' from receiving treatment to receiving control (drop-in rate) per time unit.
#' @param transP0 a numeric vector of length 2, consisting of the transition
#' probability from
#' receiving control to drop-out (drop-out rate) and
#' from receiving control to receiving treatment (drop-in rate) per time unit.
#' @param Wlist a list, consisting of weight functions applied to the test.
#' The element of the list must be functions. Default is a list of one constant
#' function, corresponding to the logrank test.
#' @param eventN the number of events
#' @param totalN the number of subjects
#' @param ratio allocation ratoi, Default: 1
#' @param alpha type i error, Default: 0.05
#' @param alternative alternative hypothesis, Default: c("two.sided", "less", "greater")
#' @param k an integer, indicating number of sub-intervals per time unit, Default: 100
#' @param nocensor PARAM_DESCRIPTION, Default: FALSE
#' @param criteria PARAM_DESCRIPTION, Default: 100
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{weighted.mean}},\code{\link[stats]{cor}}
#'  \code{\link[mvtnorm]{qmvnorm}},\code{\link[mvtnorm]{pmvnorm}}
#' @rdname n2pwr.maxLR
#' @export
#' @importFrom stats weighted.mean cov2cor
#' @importFrom mvtnorm qmvnorm pmvnorm
n2pwr.maxLR<- function(entry   = 1
                       ,fup      = 1
                       ,CtrlHaz
                       ,hazR
                       ,transP1
                       ,transP0
                       ,Wlist
                       ,eventN
                       ,totalN
                       ,ratio    = 1
                       ,alpha    = 0.05
                       ,alternative=c("two.sided","less","greater")
                       ,k        = 100
                       ,nocensor = FALSE
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
  load <- trans.mat(num,x,ctrlRate,haz_val,haz_point,ratio,
                    transP1,transP0,nocensor)
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
      set.seed(i)
      crit <- mvtnorm::qmvnorm(1-alpha,tail="both.tails",
                               mean=rep(0,wn),sigma = rho_est)$quantile
      power <- 1-mvtnorm::pmvnorm(-crit,crit,mean=mu,sigma = rho_est)
      return(c(crit,power))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    power <- ftwod[2]
  }else if (alternative=="less"){ # l1 <l0
    ftwo <- function(i){
      set.seed(i)
      crit <- mvtnorm::qmvnorm(alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
      power <- 1-mvtnorm::pmvnorm(crit,Inf,mean=mu,sigma = rho_est)[1]
      return(c(crit,power))
    }
    ftwod <- apply(do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) ,2,mean)
    crit <- ftwod[1];  power <- ftwod[2]
  }else if (alternative=="greater"){
    ftwo <- function(i){
      set.seed(i)
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
