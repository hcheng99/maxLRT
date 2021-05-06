#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: c("schoenfeld", "freedman")
#' @param lambda0 PARAM_DESCRIPTION
#' @param lambda1 PARAM_DESCRIPTION
#' @param phi PARAM_DESCRIPTION
#' @param t_enrl PARAM_DESCRIPTION, Default: 0
#' @param t_fup PARAM_DESCRIPTION
#' @param alpha PARAM_DESCRIPTION, Default: 0.05
#' @param beta PARAM_DESCRIPTION, Default: 0.1
#' @param two.side PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
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
  HR <- lambda1/lambda0
  if (t_enrl==0){
    e0 <- list(value=1-exp(-lambda0))
    e1 <- list(value=1-exp(-lambda1))

  }
  else {
    intef <- function(x,l){(1-exp(-l*x))}
    e0 <-stats::integrate(function(x){intef(x,l=lambda0)},lower=t_fup,upper=t_enrl+t_fup)
    e1 <-stats::integrate(function(x){intef(x,l=lambda1)},lower=t_fup,upper=t_enrl+t_fup)
  }
  ## calculate the event rate
  erate <- stats::weighted.mean(c(e0$value,e1$value),
                         w=c(1/(1+phi),phi/(1+phi)))/(t_enrl)
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
