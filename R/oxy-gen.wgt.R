#' @title Weight Functions for \code{MaxLRtest}
#' @description Generate commonly used weight functions for \code{MaxLRtest}
#' function
#' @param method a vector of text specifying the method(s). The method(s)
#' must be one or some of c(\code{"LR"}, \code{"FH"}, \code{"Wilcoxon"},
#'  \code{"Tarone"},\code{"Maxcombo"},\code{"Maxcross"}).  Default: c("LR")
#' @param param a vector of length 2. If \code{FH} method is used, the
#' \eqn{\rho} and \eqn{\gamma} parameters must be provided, Default: 1
#' @param theta a value within (0,1). If method \code{Maxcross} is selected,
#' \code{theta} should be specified.See details. Default: 0.5
#' @return
#'  a list of weight functions
#' @details
#' The weight function for Fleming-Harrington (FH) test is \eqn{S(t)^\rho(1-S(t)^\gamma)}
#' If \code{FH} test is used, both \eqn{\rho} and \eqn{\gamma} should be provided.
#' The weight for Tarone and Ware test is \eqn{y(t)^{1/2}}, where \eqn{y(t)} is number
#' of subjects at risk. The weight for Wilcoxon test is \eqn{y(t)}.See Klein (2003) for
#' more details. Both Maxcombo test and test proposed by Cheng and He (2021)
#' needs four weight functions. Cheng's method is more sensitive in detecting
#' crossing hazards. A nuisance parameter \code{theta} is required to be
#' specified. \code{theta} represents the CDF at crossing time point. If the hazards
#' occur in the early phase of the study, a smaller value should be chosen. The
#' default value is 0.5.
#'
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#' #logrank test
#' gen.wgt(method="LR")
#' # FH and logrank test
#' fn <- gen.wgt(method=c("FH","LR"),param=c(1,1))
#' fn[[1]](0.4)
#'  }
#' }
#' @rdname gen.wgt
#' @export
gen.wgt <- function(method=c("LR")
                    , param
                    , theta=0.5 )
  {

  if (missing(method)){stop("method must be specified.")}
  if (sum(method=="FH")>1){stop("one FH method a time.")}

  lst <- list()
  nlen <- length(method)
  for (i in 1:nlen){
    if (!method[i] %in% c("LR","FH","Wilcoxon","Tarone","Maxcombo",
                          "Maxcross")){
      stop("method must be one of 'LR','FH','Wilcoxon','Tarone',
           'Maxcombo','Maxcross'")
    }
    if ("LR" %in% method[i]){
      lst <- c(lst,LR=function(x){x^0})
    }
    if ("FH" %in% method[i]){
      if (missing(param)){stop("param must be specified for FH method")}
      if (length(param)!=2){stop("both rho and gamma must be
                              specified for FH method")}
      lst <- c(lst, FH=function(x){(1-x)^param[1]*x^(param[2])})
    }
    if ("Wilcoxon" %in% method[i]){
      lst <- c(lst,Wil=function(x){1-x})
    }
    if ("Tarone" %in% method[i]){
      lst <- c(lst,tar= function(x){(1-x)^0.5})
    }
    if ("Maxcombo" %in% method[i]){
      lst <- c(lst,list(
        conf=function(x){x^0},
        f01=function(x){x},
        f10=function(x){1-x},
        f11=function(x){(1-x)*x}
      ))
    }
    if ("Maxcross" %in% method[i]){
      if (missing(theta)){stop(
        "theta must be specified for Maxcross test."
      )
      }
      lst <- c(lst,list(
        conf=function(x){x^0},
        f01=function(x){x},
        f10=function(x){1-x},
        fcross=function(x,pp=theta){
            (x<=pp)*(-1/pp*x+1)+(x>pp)*(-1/(1-pp)*(x-pp))}
      ))

    }



  }

 return(lst)

}
