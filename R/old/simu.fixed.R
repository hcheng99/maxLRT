# fitdat <- sim_time_fixed(bsl_dis="loglogistic"
#                           ,param=c(2,12)
#                           ,totN=100
#
#                           ,fHR=fun_list[[1]]#the non proportion hazard function
#                           ,prop=0.5 #proportionof  treatment group
#                           ,T0 =18 #entry time
#                           ,Tfup=24
#                           ,seed=100
#                           ,upInt=2500
# )
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param bsl_dist PARAM_DESCRIPTION, Default: c("weibull", "loglogistic")
#' @param param PARAM_DESCRIPTION
#' @param totN PARAM_DESCRIPTION
#' @param fHR PARAM_DESCRIPTION
#' @param prop PARAM_DESCRIPTION
#' @param T0 PARAM_DESCRIPTION, Default: 18
#' @param Tfup PARAM_DESCRIPTION, Default: 24
#' @param seed PARAM_DESCRIPTION
#' @param upInt PARAM_DESCRIPTION, Default: 100
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname simu.fixed
#' @export
simu.fixed <- function(bsl_dist=c("weibull","loglogistic")
                            ,param
                            # alpha=1 corresponds to exponential
                            ,totN
                            ,fHR #the non proportion hazard function
                            ,prop #proportion of treatment group
                            ,T0 =18 #entry time
                            ,Tfup=24
                            ,seed
                            ,upInt=100


){

  # to note: n0 is the # of subjects in treatment group, which is assumed to
  # have fewer events

  n0 <- totN*prop
  n1 <- totN-n0
  trt <- c(rep(1,n0),rep(0,n1))
  a=param[1]
  b=param[2]
  if (bsl_dist=="weibull"){
    set.seed(seed)
    T_0 <- rweibull(n0,a,1/b)
    #-- get the cummulative hazard function----#
    Hf <- function(t){exp(-1* a*b*integrate( function(x){(x*b)^(a-1)*fHR(x)},0,t)$value)}
  }else if (bsl_dist=="loglogistic"){
    set.seed(seed)
    T_0 <- exp(rlogis(n0,log(b),1/a))
    Hf <- function(t){exp(-1* a/b*integrate( function(x){(x/b)^{a-1}/(1+(x/b)^a)*fHR(x)},0,t)$value)}
  }
  gen_t <- function(y){uniroot(function(x){Hf(x)-y},interval = c(0,upInt),extendInt="yes")$root}
  set.seed(seed*10)
  U1 <- runif(n1)
  T_1 <-lapply(U1, gen_t) %>%unlist%>%as.vector()
  Tm<-c(T_0,T_1)

  ## in extreme case, time is 0, add 0.1;
  Tm <- ifelse(Tm==0,0.1,Tm)

  ## for entry time
  set.seed(seed+1)
  t0 <- runif(totN,0,T0)
  ot <- t0+Tm

  dat <- data.frame(id=1:totN,ent=t0,time=Tm,trt=trt,ot=ot)
  dat$t_cnsr <- with(dat,ifelse(ot>T0+Tfup,0,1))
  dat$t_val <- with(dat,ifelse(ot>T0+Tfup,T0+Tfup-ent,Tm))


  return(dat)
}
