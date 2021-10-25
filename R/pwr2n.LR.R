pwr2n.LR <- function( method    = c("schoenfeld","freedman")
                          ,lambda0
                          ,lambda1
                          ,ratio  =  1
                          ,entry  = 0
                          ,fup
                          ,alpha = 0.05
                          ,beta  = 0.1
                          ,alternative = c("two.sided")
                          ,Lparam=NULL
                          ,summary = TRUE
)
  # the calculation assumes the exponential distribution for control and
  # treatment group
  # ratio: randomization ratio
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
      if (entry==0){
        e1 <-F1(fup)
      }else {
        e1 <- stats::integrate(F1,lower=fup,upper=entry+fup)$value/entry
      }

      return(e1)
    }
    }
  HR <- lambda1/lambda0
  if (entry==0){
    if (is.null(Lparam)){
      e0 <- 1-exp(-lambda0*fup)
      e1 <- 1-exp(-lambda1*fup)
    }else {
      e0 <- calp(lambda0)
      e1 <- calp(lambda1)
    }


  }
  else {
    if (is.null(Lparam)){
      intef <- function(x,l){(1-exp(-l*x))}
      e0 <-stats::integrate(function(x){intef(x,l=lambda0)},lower=fup,upper=entry+fup)
      e1 <-stats::integrate(function(x){intef(x,l=lambda1)},lower=fup,upper=entry+fup)
      e0 <- e0$value/entry
      e1 <- e1$value/entry
    }else {
      e0 <- calp(lambda0)
      e1 <- calp(lambda1)
    }



  }
  ## calculate the event rate
  erate <- stats::weighted.mean(c(e0,e1),
                         w=c(1/(1+ratio),ratio/(1+ratio)))
  ## calculate the number of events
 if (alternative == "two.sided"){
   numerator=(stats::qnorm(1-alpha/2)+stats::qnorm(1-beta))^2
 } else if (alternative == "one.sided"){
   numerator=(stats::qnorm(1-alpha)+stats::qnorm(1-beta))^2
 } else {
   stop ("alternative must be either 'two-sided' or 'one-sided'!")
 }
  if (method == "schoenfeld"){
    Dnum=numerator*(1+ratio)^2/ratio/log(HR)^2

  }else if (method=="freedman"){
    Dnum=numerator*(1+HR*ratio)^2/ratio/(1-HR)^2
  }

  N=Dnum /erate

  if(summary ==TURE){
    cat("-----Summary of the Input Parameters----- \n")
    inparam <- c("Method", "Lambda1/Lambda0","Entry Time", "Follow-up Time",
                 "Allocation Ratio", "Type I Error", "Type II Error",
                 "Alternative","Drop-out Parameter")
    if (is.null(Lparam)) {Lapram <- NA}
    inval <- c(method, paste0(lambda1,"/",lambda0),fup,ratio, alpha, beta,
               alternative,Lparam)
    inputdata <- data.frame(parameter=inparam, value=inval)
    print(inputdata, row.names = FALSE)
    cat("-----Summary of the Output Parameters----- \n ")
    outparam <- c("Number of Events", "Number of Total Sampe Size",
                  "Overall Event Rate")
    outval <- c(Dnum, N, Dnum/N)
    outputdata <- data.frame(parameter=outparam, value=outval)
    print(outputdata, row.names = FALSE)
  }
  return(list(eventN=Dnum,totalN=N))
}

cal_event <- function(
  ratio
  ,lambda1
  ,lambda0
  ,entry
  ,fup
  ,l_shape
  ,l_scale
){
  calp <- function(lambda1){
    f1 <- function(x){lambda1*exp(-x*lambda1-(x/l_scale)^l_shape)}
    F1 <- Vectorize(function(u){integrate(f1,lower=0,upper=u)$value})
    e1 <- integrate(F1,lower=fup,upper=entry+fup)$value/entry
    return(e1)
  }

  ## treatment group
  e1 <- calp(lambda1)
  ## placebo
  e0 <- calp(lambda0)
  e <- ratio/(1+ratio)*e1+e0/(1+ratio)
  return(e)
}

