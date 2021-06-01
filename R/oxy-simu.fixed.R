# fitdat <- sim_time_fixed(bsl_dis="loglogistic"
#                           ,param=c(2,12)
#                           ,totN=100
#
#                           ,fHR=fun_list[[1]]#the non proportion hazard function
#                           ,prop=0.5 #proportionof  treatment group
#                           ,entry =18 #entry time
#                           ,fup=24
#                           ,seed=100
#                           ,upInt=2500
# )
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param bsl_dist PARAM_DESCRIPTION, Default: c("weibull", "loglogistic")
#' @param param PARAM_DESCRIPTION
#' @param totN PARAM_DESCRIPTION
#' @param fHR PARAM_DESCRIPTION
#' @param ratio PARAM_DESCRIPTION
#' @param entry PARAM_DESCRIPTION, Default: 18
#' @param fup PARAM_DESCRIPTION, Default: 24
#' @param seed PARAM_DESCRIPTION
#' @param upInt PARAM_DESCRIPTION, Default: 100
#' @param drop_param PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{Weibull}},\code{\link[stats]{integrate}},\code{\link[stats]{Logistic}},\code{\link[stats]{uniroot}},\code{\link[stats]{Uniform}}
#' @rdname simu.fixed
#' @export 
#' @importFrom stats rweibull integrate rlogis uniroot runif
simu.fixed <- function(bsl_dist=c("weibull","loglogistic")
                            ,param
                            # alpha=1 corresponds to exponential
                            ,totN
                            ,fHR #the non proportion hazard function
                            ,ratio # # of trt/# of placebo
                            ,entry =18 #entry time
                            ,fup=24
                            ,seed
                            ,upInt=100
                            ,drop_param


){

  # to note: n1 is the # of subjects in treatment group, which is assumed to
  # have fewer events
  prop <- ratio/(ratio+1)
  n1 <- ceiling(totN*prop)
  n0 <- totN-n1
  if (ceiling(n1)!=n1|ceiling(n0)!=n0){
    stop("The number of subjects in each group must be interger. Check!")
  }
  trt <- c(rep(0,n0),rep(1,n1))
  a=param[1]
  b=param[2]

  if (bsl_dist=="weibull"){
    set.seed(seed)
    T_0 <- stats::rweibull(n0,a,1/b)
    #-- get the cummulative hazard and survival function----#
    Hf <- function(t){exp(-1* a*b*stats::integrate( function(x){(x*b)^(a-1)*fHR(x)},0,t)$value)}
  }else if (bsl_dist=="loglogistic"){
    set.seed(seed)
    T_0 <- exp(stats::rlogis(n0,log(b),1/a))
    Hf <- function(t){exp(-1* a/b*integrate( function(x){(x/b)^{a-1}/(1+(x/b)^a)*fHR(x)},0,t)$value)}
  }
  gen_t <- function(y){stats::uniroot(function(x){Hf(x)-y},interval = c(0,upInt),extendInt="yes")$root}
  set.seed(seed*10)
  U1 <- stats::runif(n1)
  T_1 <-lapply(U1, gen_t) %>%unlist%>%as.vector()
  Tm<-c(T_0,T_1)

  ## in extreme case, time is 0, add 0.1;
  Tm <- ifelse(Tm==0,0.1,Tm)
  if (!missing(drop_param)){
    a_drop <- drop_param[1]
    b_drop <- drop_param[2]
    set.seed(seed*99)
    drop_time <- stats::rweibull(totN,a_drop,1/b_drop)
  }else{
    drop_time <- rep(Inf,totN)
  }

  ## for entry time
  set.seed(seed+1)
  t0 <- runif(totN,0,entry)
  ot <- t0+Tm
  Dur <- entry+fup # the length of study
  dat <- data.frame(id=1:totN,ent=t0,time=Tm,trt=trt,ot=ot,
                    drop_time=drop_time,ot_drop=t0+drop_time)
  # find the smallest time
  min_ind <- apply(cbind(dat$ot,Dur,dat$ot_drop),1,which.min)
  status <- c("event","end of study","drop-out")
  tot_len <- cbind(dat$ot,Dur,dat$ot_drop)[cbind(seq_along(min_ind),min_ind)]
  dat$t_val <- tot_len-t0
  dat$cnsr_desc <- status[min_ind]
  dat$cnsr <- 1-(min_ind==1) #cnsr=1 indicates censoring

  final <- with (dat,data.frame(
    id=id
    ,group=trt
    ,entry.time=ent
    ,event.time=time
    ,drop.time=drop_time
    ,obs.time=tot_len
    ,aval=t_val
    ,cnsr=cnsr
    ,cnsr.desc=cnsr_desc
     ,event=1-cnsr
  )
  )
#table(final$group,final$cnsr.desc)

  return(final)
}
