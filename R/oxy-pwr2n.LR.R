#' @title Sample Size Calculation under Proportional Hazards
#' @description \code{pwr2n.LR} calculates the total number of events and total
#' number of subjects required given the provided design parameters based on either
#' schoenfeld or freedman formula.
#' @param method calculation formula, Default: c("schoenfeld", "freedman")
#' @param lambda0 hazard rate for the control group
#' @param lambda1 hazard rate for the treatment group
#' @param phi randomization ratio between treatment and control. For example,
#'  phi=2 if randomization ratio is 2:1 to treatment and control group.
#' @param t_enrl enrollment time. A constant enrollment rate is assumed,
#'  Default: 0
#' @param t_fup follow-up time.
#' @param alpha type I error rate, Default: 0.05
#' @param beta type II error rate. For example,if the target power is 80%, beta is 0.2.
#'  Default: 0.1
#' @param two.side a logical scalar; if TRUE then a two-sided test is used.If
#'  FALSE, a one-sided test is used. Default: TRUE
#' @param Lparam a vector of parameters for the drop-out Weibull distribution,
#'  See Details below. Default: NULL
#' @return
#'  \item{eventN}{The total number of events}
#' \item{totalN}{The total number of subjects}
#' @details The total event number is determined only by \eqn{\alpha, \beta} and
#'  hazard ratio, i.e., \eqn{\lambda_1/\lambda_0}. Other design parameters such as
#'  enrollment period affects the event probability and thus the total sample size.
#'  A fixed duration design is assumed in the calculation. ALl patients are enrolled
#'  at a constant rate within \code{t_enrl} time and have at least \code{t_fup}
#'  time of follow-up. So the total study duration is \code{t_enrl}+\code{t_fup}.
#'  If drop-out is expected, a Weibull distribution with shape parameter -\eqn{\alpha}
#'  and scale parameter - \eqn{\beta} is considered. The CDF is
#'  \eqn{F(x)=1-exp(-(x/\beta)^\alpha)}.
#' @author   Huan Cheng (hcheng1118@@hotmail.com)
#' @references
#' Schoenfeld, D. (1981) The asymptotic properties of nonparametric
#'  tests for comparing survival distributions. Biometrika, 68,
#' 316–319.
#'
#' Freedman, L. S. (1982) Tables of the number of patients required
#' in clinical trials using the logrank test. Statistics in medicine, 1, 121–129.

#' @examples
#'l0 <- log(2)/14; HR <- 0.8; RR <- 2; t_enrl <- 12; t_fup <- 12;
#' eg1 <- pwr2n.LR( method    = c("schoenfeld")
#'                  ,l0
#'                  ,l0*HR
#'                  ,phi=RR
#'                  ,t_enrl
#'                  ,t_fup
#'                  ,alpha     = 0.05
#'                  ,beta      = 0.1
#'                  ,two.side  = TRUE
#' )
#' # event number, total subjects, event probability
#' c(eg1$eventN,eg1$totalN,eg1$eventN/eg1$totalN)
#'
#' # example 2: drop-out from an exponential with median time is 30
#' eg2 <- pwr2n.LR( method    = c("schoenfeld")
#'                  ,l0
#'                  ,l0*HR
#'                  ,phi=RR
#'                  ,t_enrl
#'                  ,t_fup
#'                  ,alpha     = 0.05
#'                  ,beta      = 0.1
#'                  ,two.side  = TRUE
#'                  ,Lparam = c(1,30/log(2))
#' )
#' # event number, total subjects, event probability
#' c(eg2$eventN,eg2$totalN,eg2$eventN/eg2$totalN)
#'
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{integrate}},\code{\link[stats]{weighted.mean}},\code{\link[stats]{Normal}}
#' @rdname pwr2n.LR
#' @export
#' @importFrom stats integrate weighted.mean qnorm
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
  print(erate)
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phi PARAM_DESCRIPTION
#' @param lambda1 PARAM_DESCRIPTION
#' @param lambda0 PARAM_DESCRIPTION
#' @param t_enrl PARAM_DESCRIPTION
#' @param t_fup PARAM_DESCRIPTION
#' @param l_shape PARAM_DESCRIPTION
#' @param l_scale PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cal_event
#' @export
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

