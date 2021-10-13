###input parameters
#' @title Sample Size Calculation with Maximum Weighted Logrank Test
#'
#' @description \code{pwr2n.maxLR} calculates the number of events and
#' subjects required to achieve pre-specified power in the setup of two groups.
#' The method extends the calculation in the framework of the Markov model by Lakatos, allowing
#' for using the maximum weighted logrank tests with an arbitrary number of weight
#' functions. If only one weight function is provided, the test is essentially
#' the classic (weighted) logrank test.
#'
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
#' @param entry_pdf1 a pdf of enrollment time for treatment
#' group. See \code{entry_pdf0}, Default: assume same pdf as control group.
#' @param ratio an integer, indicating the randomization ratio between treatment
#' and control group, Default: 1
#' @param alpha type I error rate, Default: 0.05
#' @param beta type II error rate, Default: 0.1
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "\code{two.sided}", "\code{greater}","\code{less}". See details.
#' Default: c("\code{two.sided}")
#' @param criteria an integer indicating the maximum iteration allowed in
#' obtaning the number of events. See details , Default: 500
#' @param k an integer, indicating number of sub-intervals per time unit,
#'  Default: 100
#' @return
#' An object of class "\code{maxLRT}" with corresponding \code{plot} function.
#' The object is a list containing the following components:
#'  \item{eventN}{total number of events}
#'  \item{totalN}{total number of subjects}
#'  \item{pwr}{actual power given the number of events}
#'  \item{prob_event}{event probability at the end of trial}
#'  \item{L_trans}{a list, consisting of tranisition matrix at each inteval}
#'  \item{pdat}{ a data frame including all the intermediate variables in the calculation.
#'  see Details.}
#'  \item{studytime}{a vector of length 2, including the entry and follow-up time as input}
#'  \item{RandomizationRatio}{as input}
#' @details
#' The detailed calculation procedure can be found in the reference paper. The number
#' of subjects is determined by several factors, including the control hazard, hazard
#' ratio function, entry time and distribution, follow-up time, etc. The number of
#' events is mainly determined by the control hazard and hazard ratio function.
#' The study design assumes \code{entry} time units of
#' enrollment and at least \code{fup} time units of follow-up. If enrollment
#' time \code{entry} is set to zero, all subjects are enrolled simultaneously,
#' so there is no staggered entry. Otherwise, if
#' \code{entry} is greater than 0, administrative censoring is considered. The user-defined
#'enrollment time function, hazard function for the control group and hazard ratio function can be either discrete or continuous.
#'Various non-proportional hazards types are accommodated. See examples below.
#'If multiple weight functions are provided in \code{Wlist}, a maximum weighted logrank
#'test or combination test is implemented. An iterative procedure
#'is used to obtain the event number based on the multivariate normal distribution.  Package
#'\pkg{mvtnorm} is used to calculate the quantiles. Because the algorithm is slightly
#'seed dependent, the quantiles are mean values of ten replicates.
#'
#'The "\code{alternative}" option supports both two-sided and one-sided test.
#' Let \eqn{\Lambda_1} and \eqn{\Lambda_0} denote the cumulative hazard of
#' treatment and control group. The \code{less} option tests
#' \eqn{H_0: \Lambda_1 > \Lambda_0} against
#' \eqn{H_a: \Lambda_1 <= \Lambda_0}. The \code{greater} option tests
#'\eqn{H_0: \Lambda_1 < \Lambda_0} against \eqn{H_a: \Lambda_1 >= \Lambda_0}.
#'
#' @references
#' Bender, R., Augustin, T., & Blettner, M. (2005). Generating survival times to simulate Cox proportional
#' hazards models. Statistics in medicine, 24(11), 1713-1723.
#' @examples
#' \dontrun{
#' if(interactive()){
#'### Define weight functions weight functions
#'timef1 <- function(x){1}
#'timef2 <- function(x){(x)}
#'timef3 <- function(x){1-x}
#'timef5 <- function(x,pp=1/2){(x<=pp)*(-1/pp*x+1)+(x>pp)*(-1/(1-pp)*(x-pp))}
#'W_4e <- list(timef1,timef2,timef3,timef5)
#'t_enrl <- 5
#'t_fup <- 5
#'lmd0 <- -log(0.2)/10
#'f_hr <- function(x){0.5*x^0}
#'ratio <-1
#'alpha <- 0.05
#'beta <- 0.1
#'## uniform enrollment, proportional hazards and logrank test
#'ef <- function(x){(1/t_enrl)*(x>0&x<=t_enrl)}
#'size1 <- maxLRT::pwr2n.maxLR(entry     = t_enrl
#'                             ,fup      = t_fup
#'                             ,k        = 100
#'                             ,ratio    = ratio
#'                             ,Wlist  = list(timef1)
#'                             ,CtrlHaz=function(x){lmd0*x^0}
#'                             ,transP1=c(0,0)
#'                             ,transP0=c(0,0)
#'                             ,hazR     = f_hr
#'                             ,alpha    = alpha
#'                             ,beta     = beta
#'                             ,entry_pdf0=ef
#'                             ,entry_pdf1 =ef
#')
#'
#'target_N <- round(size1$totalN,digits=0)
#'target_E <- round(size1$eventN,digits=0)
#'c(target_E,target_N)
#'
#'## delayed response with maximum weighted logrank test
#'
#'t_enry <- 12
#'t_fup <- 18
#'lmd0 <- log(2)/12
#'f_hr_delay <- function(x){(x<=6)+(x>6)*0.75}
#'size <- maxLRT::pwr2n.maxLR(entry     = t_enrl
#'                            ,fup      = t_fup
#'                            ,k        = 100
#'                            ,ratio    = ratio
#'                            ,Wlist  = W_4e
#'                            ,CtrlHaz=function(x){lmd0*x^0}
#'                            ,transP1=c(0,0)
#'                            ,transP0=c(0,0)
#'                            ,hazR     = f_hr_delay
#'                            ,alpha    = alpha
#'                            ,beta     = beta
#'
#'
#')
#'target_N <- round(size$totalN,digits=0)
#'target_E <- round(size$eventN,digits=0)
#'#'c(target_E,target_N)#'
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{cor}},\code{\link[stats]{weighted.mean}}
#'  \code{\link[mvtnorm]{qmvnorm}},\code{\link[mvtnorm]{pmvnorm}}
#' @rdname pwr2n.maxLR
#' @export
#' @importFrom stats cov2cor weighted.mean
#' @importFrom mvtnorm qmvnorm pmvnorm
pwr2n.maxLR<- function(entry   = 1
                       ,fup      = 1
                       ,CtrlHaz
                       ,hazR
                       ,transP1
                       ,transP0
                       ,Wlist = list(function(x){ x^0})
                       ,entry_pdf0=function(x){(1/entry)*(x>=0&x<=entry)}
                       ,entry_pdf1=entry_pdf0
                       ,ratio    = 1
                       ,alpha    = 0.05
                       ,beta     = 0.1
                       ,alternative = c("two.sided")
                       ,criteria = 500
                       ,k        = 100
){

  if (!alternative %in% c("two.sided","greater","less")){
    stop("The alternative must be one of 'two.sided','greater','less'.")
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
  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(pdat),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)
  event <- c()
  if (alternative=="two.sided"){
    part2=(qnorm(1-alpha/2)+qnorm(1-beta))^2
  }else {
    part2=(qnorm(1-alpha)+qnorm(1-beta))^2
  }
  for (j in 1:wn){
    W[,j] <- Wlist[[j]](pdat$S)
    wgt <- W[,j]
    part1=t(pdat$rho)%*%(pdat$eta*wgt^2)
    part3=t(pdat$rho)%*%(pdat$gamma*wgt)
    part3=part3^2
    event[j] <- part1*part2/part3
  }
  #print(c(part1,part2,part3))
  ## get the min and max sample size
  E_min=min(event)
  E_max=max(event)
  eta=-1
  dnum <- E_min
  count <- 0
  while(eta <0){
    count <- count+1
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
      ftwod <- do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) %>% apply(2,mean)
      power <- ftwod[2]
      crit <- ftwod[1]
    }else if (alternative=="less"){ # l1 <l0
      ftwo <- function(i){
        set.seed(i)
        crit <- mvtnorm::qmvnorm(alpha,tail="lower.tail",mean=rep(0,wn),sigma = rho_est)$quantile
        power <- 1-mvtnorm::pmvnorm(crit,Inf,mean=mu,sigma = rho_est)[1]
        return(c(crit,power))
      }
      ftwod <- do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) %>% apply(2,mean)
      crit <- ftwod[1];  power <- ftwod[2]
    }else if (alternative=="greater"){
      ftwo <- function(i){
        set.seed(i)
        crit <- qmvnorm(alpha,tail="upper.tail",mean=rep(0,wn),sigma = rho_est)$quantile
        power <- 1-pmvnorm(-Inf,crit,mean=mu,sigma = rho_est)[1]
        return(c(crit,power))
      }
      ftwod <- do.call(rbind,sapply(1:10,ftwo,simplify=FALSE)) %>% apply(2,mean)
      crit <- ftwod[1];  power <- ftwod[2]
    }

    if (power<1-beta&count<criteria) {dnum<-dnum+1}
    else if (count>=criteria){
      warning(paste0("the algoritm doesn't converge within ",criteria," iterations;
                     the current power is ",power,"; event size: ",dnum))
      ;break}
    else {break}

  }

  eprob <- stats::weighted.mean(c(pdat$C_E[num],pdat$E_E[num]),w=c(1,ratio))
  Nsize <- dnum/eprob


  listall <- list( eventN  = dnum
                   ,totalN = Nsize
                   ,pwr = as.numeric(power)
                   ,prob_event =eprob
                   ,L_trans = load$L_trans
                   ,pdat = pdat
                   ,studytime=c(entry,fup)
                   ,RandomizationRatio=ratio

  )
  class(listall) <-"MaxLRpwr"
  return(listall)


}

#*********************************************
#*show the survival plot/ hazards plots
#*********************************************
#' @title Graphical Display of Design Parameters in Sample Size Calculation
#' @description Displays graphs of survival, hazards, drop-out and censor over time
#' as specified in the calculation.
#' @param x object of the \code{pwr2n.maxLR} function
#' @param type a vector of string, specifying the graphs to display. The options
#' include "\code{hazard}","\code{survival}","\code{dropout}","\code{event}", and
#' "\code{censor}". If \code{type} is not provided, all the available graphs are
#' generated.
#'
#' @param ... additional graphical arguments passed to the plot function
#' @return
#' plots are produced on the current graphics device
#' @details
#'The \code{type} argument provides five options to visualize the trial in design.
#'Option \code{survival} shows the survival probabilities of treatment and control
#'group over time. \code{hazard} option provides the hazard rates and hazard ratio
#'over time.\code{dropout} shows the proportion of drop-out subjects across the trial duration.
#'\code{censor} shows the proportion of censored subjects over time.
#'
#' @examples
#' \dontrun{
#' if(interactive()){

#'timef1 <- function(x){1}
#'timef2 <- function(x){(x)}
#'timef5 <- function(x,pp=1/2){(x<=pp)*(-1/pp*x+1)+(x>pp)*(-1/(1-pp)*(x-pp))}
#'W_4e <- list(timef1,timef2,timef3,timef5)
#'t_enrl <- 5
#'t_fup <- 5
#'lmd0 <- -log(0.2)/10
#'f_hr <- function(x){0.5*x^0}
#'ratio <-1
#'alpha <- 0.05
#'beta <- 0.1
#'### crossing hazards with maximum weighted logrank test
#'f_hr_cross <- function(x){(x<=6)*1.2+(x>6)*0.75}
#'size <- maxLRT::pwr2n.maxLR(entry     = t_enrl
#'                            ,fup      = t_fup
#'                            ,k        = 100
#'                            ,ratio    = ratio
#'                            ,Wlist  = W_4e
#'                            ,CtrlHaz=function(x){lmd0*x^0}
#'                            ,transP1=c(0.00,0)
#'                            ,transP0=c(0.00,0)
#'                            ,hazR     = f_hr_cross
#'                            ,alpha    = alpha
#'                            ,beta     = beta
#'
#'
#')
#'## only plot hazards functions
#'plot(size,type="hazard")
#'## display all graphs
#'plot(size)
#'  }
#' }
#' @seealso
#'  \code{\link[maxLRT]{pwr2n.maxLR}}
#' @rdname plot.MaxLRpwr
#' @export
#' @importFrom graphics lines legend par
plot.MaxLRpwr<- function(x,type=c("hazard","survival","dropout","event","censor"),...) {
  datM <- x$pdat
  totalN <- x$totalN
  ratio <- x$RandomizationRatio
  tval <- 1
  if( missing(type)){ tval <- 0}
  ## draw the survival curves
  if (tval==0|"survival" %in% type ){
    with(datM,{
      plot(ti,S1,cex=0.1,lty=1,ylim=c(0,1),col=1,xlab="Time",
           ylab="Survival Probability",main="Survival Curves",...)
      graphics::lines(ti,S0,col=2,lty=2,...)
      graphics::legend("bottomleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)

    })
  }

  ## hazard functions
  if (tval==0|"hazard" %in% type ){
    graphics::par(mfrow=c(1,2))
    ymax <- max(c(datM$hazard_C,datM$hazard_E))*1.1
    with(datM,{
      plot(ti,hazard_E,cex=0.1,lty=1,col=1,xlab="Time",ylab="Hazard Rate",
           ylim=c(0,ymax),main="Hazard Curves",...)
      graphics::lines(ti,hazard_C,col=2,lty=2,...)
      graphics::legend("bottomright",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
    ## hazard ratio
    with(datM,{
      plot(ti,theta,cex=0.1,lty=1,col=1,xlab="Time",
           ylim=c(min(theta)*0.9,max(theta)*1.1),
           ylab="Hazard Ratio (treatment over control)",
           main="Hazard Ratio over Time",...)


    })

  }
  ## drop out
  if (tval==0|"dropout" %in% type){
    ymax <- max(c(datM$E_L,datM$C_L))*1.1
    with(datM,{
      plot(ti,E_L,cex=0.1,lty=1,col=1,xlab="Time",ylab="proporion of drop out",
           ylim=c(0,ymax),main="Drop-out overtime",...)
      graphics::lines(ti,C_L,col=2,lty=2,...)
      graphics::legend("topleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
  }
  ## censor
  if (tval==0|"censor" %in% type){
    ymax <- max(c(datM$E_C,datM$C_C))*1.1
    with(datM,{
      plot(ti,E_C,cex=0.1,lty=1,col=1,xlab="Time",
           ylab="proporion of censor",
           ylim=c(0,ymax),main="Administrative censoring overtime",...)
      graphics::lines(ti,C_C,col=2,lty=2,...)
      graphics::legend("topleft",legend=c("treatment","control"),
                       col=1:2,lty=1:2,cex=0.8)
    })
  }
  ##event number
  if (tval==0|"event" %in% type ){
    with(datM,{
      plot(ti,round(eprob*totalN,digits=0),cex=0.1,lty=1,col=1,xlab="Time",
           ylab="number of events",
           main="Events over time",...)
      graphics::lines(ti,round(E_E*totalN*ratio/(ratio+1),digits=0),col=2,lty=2,...)
      graphics::lines(ti,round(C_E*totalN/(ratio+1),digits=0),col=3,lty=3,...)
      graphics::legend("topleft",legend=c("overall","treatment","control"),
                       col=1:3,lty=1:3,cex=0.8)
    })
  }


  #
  #   graphics::mtext(paste(paste0("Statistic",1:wn,sep=":"),round(x$stat.mat[,1],3),
  #                         collapse = "  "),side=3,cex=0.7)
}
