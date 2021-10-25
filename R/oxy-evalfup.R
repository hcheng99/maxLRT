##################
## evaluate the relationship between follow-up time and sample size
###################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param lower.time PARAM_DESCRIPTION
#' @param upper.time PARAM_DESCRIPTION
#' @param size PARAM_DESCRIPTION
#' @param increment PARAM_DESCRIPTION, Default: 0.5
#' @param xlabel PARAM_DESCRIPTION, Default: 'Follow-up Time'
#' @param ylabl PARAM_DESCRIPTION, Default: 'Total Sample Size'
#' @param title PARAM_DESCRIPTION, Default: 'Relationship between Follow-up and
#' Total Sample Size'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname evalfup
#' @export
evalfup <- function(object, lower.time, upper.time, size,
                    increment=0.5, xlabel = "Follow-up Time",
                    ylabl = "Total Sample Size",
                    title = "Relationship between Follow-up and \n Total Sample Size"
                    ){
  fupseq <- seq(from=lower.time, to=upper.time, by=increment)
  N <- c()
  D <- c()
  for (i in 1: length(fupseq)){
    tmp <- pwr2n.maxLR(entry     = object$studytime[1]
                ,fup      = fupseq[i]
                ,k        = object$inputfun$k
                ,ratio    = object$RandomizationRatio
                ,Wlist  = object$inputfun$Weightfunctions
                ,CtrlHaz=object$inputfun$controalhazard
                ,transP1=object$inputfun$transP1
                ,transP0=object$inputfun$transP0
                ,hazR     = object$inputfun$hazardratio
                ,alpha    = object$inputfun$alpha
                ,beta     = object$inputfun$beta
                ,entry_pdf0=object$inputfun$entrypdf0
                ,entry_pdf1 =object$inputfun$entry_pdf1
                ,summary=FALSE
    )
    N <-c(N,tmp$totalN)
    D <- c(D,tmp$eventN)
  }
  Nint <- ceiling(N/(object$RandomizationRatio+1))*(object$RandomizationRatio+1)
  ymax <- max(size,Nint)
  ymin <- min(size,Nint)
  if (ymax<size){
    cat("The largest required sample size (",ymax,") is less than the target
        sampe size", size, ". Consider reduce the follow-up time")
  }else if (ymin >size){
    cat("The smallest required sample size (",ymin,") is greater than the target
        sampe size", size, ". Consider increase the follow-up time")
  }
  interp <- approx(fupseq,y=Nint,method="linear")
  y1 <- interp$y
  x1 <- interp$x
  p1 <- max(which(y1 >= size))
  p2 <- min(which(y1  <= size))
  xsize <- mean(c(x1[p1:p2]))
  move <- (upper.time-lower.time)/8
  plot(x=fupseq,y=Nint,ylim=c(ymin,ymax),cex=0.8,col="red",
       xlab= xlabel, ylab=ylabel,
       main = title)
  points(interp$x,interp$y,pch=16,cex=0.3)
  segments(x0=xsize,x1=xsize,y0=0,y1=size,lty=2)
  segments(x0=0,x1=xsize,y0=size,y1=size,lty=2)
  text(paste("(",round(xsize,digits = 1),",",size,")"),x=xsize+move,
  y=size,cex=0.7)
  legend("topright",legend=c("Orignal","Interpolated"),pch=c(1,16),
         col=c("red","black"),cex=c(0.8))
 # return values
  original <- list(x=fuptime,y=N)
  retrun(list(
    approx.time = xsize,
    original = original,
    interp = interp,
    Esize = D

  ))
}



