###input parameters
#' @title Power Calculation with Maximum Weighted Logrank Test
#' @description  \code{n2pwr.NPH} calculates the power given either the
#' number of events or number of subjects
#' @param method a text specifying the calculation method, either
#' \code{"MaxLR"} or \code{"Projection"}. Maximum weighted
#' logrank test is used if \code{"MaxLR"} is specified; otherwise,
#' projection test is used.
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
#' @param entry_pdf0 a function, indicating the probability density function (pdf)
#' of enrollment time for control group. The default assumes a uniform distribution
#' corresponding to the constant enrollment rate.
#' Default: function(x) {
#'    (1/entry) * (x >= 0 & x <= entry)
#'}
#' @param entry_pdf1 a pdf functionof enrollment time for treatment
#' @param eventN the number of events
#' @param totalN the number of subjects
#' @param ratio allocation ratio, Default: 1
#' @param alpha type i error, Default: 0.05
#' @param alternative alternative hypothesis, Default: c("two.sided", "less", "greater")
#' @param k an integer, indicating number of sub-intervals per time unit, Default: 100
#' @param criteria an integer indicating the maximum iteration allowed in
#' obtaining the number of events. See details , Default: 500
#' @return
#' a list of components:
#' \item{power}{asymptotic power }
#' \item{eventN}{number of events. If this is the input, it is the input
#' number. If totalN is specified, eventN is calculated based on provided
#' parameters}
#' \item{totalN}{number of subjects. If eventN is specified, totalN is calculated
#' based on provided parameters. Otherwise, it the same as input.}
#' \item{prob_event}{event probability at the end of trial}
#' \item{L_trans}{a list, consisting of tranisition matrix at each inteval}
#' \item{pdat}{ a data frame including all the intermediate variables in the calculation.
#' }
#'  \item{studytime}{a vector of length 2, including the entry and follow-up time as input}
#'  \item{RandomizationRatio}{as input}
#' @details
#' Function \code{npwr.maxLR} calculate the asymptotic power given number
#' of events or number of subjects using maximum weighted logrank test
#' iteratively. Check function \code{pwr2n.maxLR} for more details.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  t_enrl <- 12
#'t_fup <- 18
#'lmd0 <- -log(0.2)/10
#'#delayed treatment effects
#'f_hr_delay <- function(x){(x<=6)+(x>6)*0.75}
#'ss <- 1
#'ratio <-1
#'alpha <- 0.05
#'beta <- 0.1
#'maxc <- gen.wgt(method="Maxcombo")
#'pwr1 <- n2pwr.NPH(entry   = t_enrl
#'                    ,fup      = t_fup
#'                  ,CtrlHaz = function(x){x^0*lmd0}
#'                  ,hazR = f_hr_delay
#'                  ,transP1 = c(0,0)
#'                  ,transP0 = c(0,0)
#'                  ,Wlist = maxc
#'                  ,eventN = 50 # targeted number of events
#'                  ,ratio    = 1
#'                  ,alpha    = 0.05
#'                  ,alternative=c("two.sided")
#'                  ,k        = 100
#'                  ,criteria = 100
#')
#'  }
#' }
#' @seealso
#'  \code{\link{pwr2n.NPH}}
#' @rdname n2pwr.NPH
#' @export
#' @importFrom mvtnorm qmvnorm pmvnorm
n2pwr.NPH<- function(method = "MaxLR"
                       ,entry   = 1
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
  load <- trans.mat(numN=num,x=x,ctrlRate=ctrlRate,haz_val=haz_val,
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
  if (method == "MaxLR"){
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
  }
  else if (method == "Projection"){
    if (alternative!="two.sided"){
      cat(c("note: only two-sided is supported for projection test."))
    }
    ## get the rank of the variance matrix
    mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))
    vr <- qr(Vmat)$rank
    crit <- stats::qchisq(1-alpha,df=vr)
    ## get the noncentral parameter
    lmd <- t(mu)%*%MASS::ginv(Vmat)%*%mu
    power <- stats::pchisq(crit,df=vr,ncp=lmd,lower.tail = FALSE)

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
