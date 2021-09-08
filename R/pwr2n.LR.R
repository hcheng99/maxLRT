pwr2n.LR <- function( method    = c("schoenfeld","freedman")
                          ,lambda0
                          ,lambda1
                          ,phi
                          ,t_enrl    =0
                          ,t_fup
                          ,alpha     = 0.05
                          ,beta      = 0.1
                          ,two.side  = TRUE
                          ,Lparam=NULL
)
  # the calculation assumes the exponential distribution for control and
  # treatment group
  # phi: randomization ratio
  #
{
  if (missing(method)) {stop("method must be specified")}
  if (!(match.arg(method)  %in% (c("schoenfeld","freedman")))){
    stop("Error: method must be schoenfeld or freedman")
  }
  if (!is.null(Lparam)){
    l_shape=Lparam[1]; l_scale=Lparam[2]
    calp <- function(lambda1){
      f1 <- function(x){lambda1*exp(-x*lambda1-(x/l_scale)^l_shape)}
      F1 <- Vectorize(function(u){integrate(f1,lower=0,upper=u)$value})
      if (t_enrl==0){
        e1 <-F1(t_fup)
      }else {
        e1 <- stats::integrate(F1,lower=t_fup,upper=t_enrl+t_fup)$value/t_enrl
      }

      return(e1)
    }
    }
  HR <- lambda1/lambda0
  if (t_enrl==0){
    if (is.null(Lparam)){
      e0 <- 1-exp(-lambda0*t_fup)
      e1 <- 1-exp(-lambda1*t_fup)
    }else {
      e0 <- calp(lambda0)
      e1 <- calp(lambda1)
    }


  }
  else {
    if (is.null(Lparam)){
      intef <- function(x,l){(1-exp(-l*x))}
      e0 <-stats::integrate(function(x){intef(x,l=lambda0)},lower=t_fup,upper=t_enrl+t_fup)
      e1 <-stats::integrate(function(x){intef(x,l=lambda1)},lower=t_fup,upper=t_enrl+t_fup)
      e0 <- e0$value/t_enrl
      e1 <- e1$value/t_enrl
    }else {
      e0 <- calp(lambda0)
      e1 <- calp(lambda1)
    }



  }
  ## calculate the event rate
  erate <- stats::weighted.mean(c(e0,e1),
                         w=c(1/(1+phi),phi/(1+phi)))
  ## calculate the number of events
  numerator=(stats::qnorm(1-alpha/2)+stats::qnorm(1-beta))^2

  if (method == "schoenfeld"){
    Dnum=numerator*(1+phi)^2/phi/log(HR)^2

  }else if (method=="freedman"){
    # numerator=(qnorm(1-alpha/2)+2*sqrt(HR)*qnorm(1-beta)/(HR+1))^2
    Dnum=numerator*(1+HR*phi)^2/phi/(1-HR)^2
    # Dnum=numerator*4*HR/(1-HR)^2

  }

  N=Dnum /erate
  return(list(eventN=Dnum,totalN=N))
}

cal_event <- function(
  phi
  ,lambda1
  ,lambda0
  ,t_enrl
  ,t_fup
  ,l_shape
  ,l_scale
){
  calp <- function(lambda1){
    f1 <- function(x){lambda1*exp(-x*lambda1-(x/l_scale)^l_shape)}
    F1 <- Vectorize(function(u){integrate(f1,lower=0,upper=u)$value})
    e1 <- integrate(F1,lower=t_fup,upper=t_enrl+t_fup)$value/t_enrl
    return(e1)
  }

  ## treatment group
  e1 <- calp(lambda1)
  ## placebo
  e0 <- calp(lambda0)
  e <- phi/(1+phi)*e1+e0/(1+phi)
  return(e)
}

