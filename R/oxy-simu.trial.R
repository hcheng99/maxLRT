
#' @title Simulate Survival Trial Data
#' @description \code{simu.trial} simulates survival data allowing flexible input
#' of design parameters. It supports both event-driven design and fixed study duration
#' design.
#' @param type indicates whether event-driven trial ("event) or fixed study duration
#' trial ("time"), Option: c("event", "time")
#' @param trial_param a vector of length 3 with components for required subject size, enrollment
#' time and required number of events ("event" type trial)/follow-up time ("time" type trial)
#' @param bsl_dist indicates the survival distribution for control group, option: c("weibull", "loglogistic")
#' @param bsl_param a vector of length 2 with the shape and scale parameter for the
#'  survival distribution of control group
#' @param drop_param a vector of length 2 with shape and scale parameter for the
#' weibull distribution of drop-out time
#' @param entry_pdf a function describing the pdf of the entry time. The default
#' is uniform enrollment, Default: function(x) {
#'    (1/trial_param[2]) * (x >= 0 & x <= trial_param[2])
#'}
#' @param HR_fun a function describing the hazard ratio function between treatment
#' and control group
#' @param allocation_ratio allocation ratio between treatment and control group.
#' For example, \code{allocation_ratio}=2 if 2:1 allocation is used.
#' @param seed a single value, specifying seed for random generation
#' @param upInt a value indicating the upper bound usedin the uniroot function.
#' See details.  Default: 100
#' @param print a logical indicating whether basic information summary is printed
#' to the console or not, Default: TRUE
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
#' @rdname simu.trial
#' @export
#' @importFrom stats rweibull integrate rlogis uniroot runif
simu.trial <- function(type=c("event","time")
                       ,trial_param # include the total sample size,entry time,
                       # target event (event type)/fup time (time )
                       ,bsl_dist=c("weibull","loglogistic")
                       ,bsl_param   # alpha=1 corresponds to exponential
                       ,drop_param
                       ,entry_pdf=function(x){(1/trial_param[2])*(x>=0&x<=trial_param[2])}
                       ,HR_fun #the non proportion hazard function
                       ,allocation_ratio # # of trt/# of placebo
                       ,seed
                       ,upInt=100
                       ,print=TRUE
){

  if (length(trial_param) !=3) {stop("The trial parameters must include
                                     total sample size, entry time, targeted events/
                                     follow-up time ")
  }else {
    t_p1 <- trial_param[1]
    t_p2 <- trial_param[2]
    t_p3 <- trial_param[3]
  }
  if (!type %in% c("event","time")){stop("The type must be in ('event','time')")}

  # to note: n1 is the # of subjects in treatment group, which is assumed to
  # have fewer events
  prop <- ratio/(ratio+1)
  n1 <- ceiling(t_p1*prop)
  n0 <- t_p1-n1
  if (ceiling(n1)!=n1|ceiling(n0)!=n0){
    stop("The number of subjects in each group must be interger. Check!")
  }
  trt <- c(rep(0,n0),rep(1,n1))
  a=bsl_param[1]
  b=bsl_param[2]
  #****************************
  #* Simulate Event Time
  #****************************
  if (bsl_dist=="weibull"){
    set.seed(seed)
    T_0 <- stats::rweibull(n0,a,1/b)
    #-- get the cummulative hazard and survival function----#
    Hf <- function(t){exp(-1* a*b*stats::integrate( function(x){(x*b)^(a-1)*HR_fun(x)},0,t)$value)}
  }else if (bsl_dist=="loglogistic"){
    set.seed(seed)
    T_0 <- exp(stats::rlogis(n0,log(b),1/a))
    Hf <- function(t){exp(-1* a/b*integrate( function(x){(x/b)^{a-1}/(1+(x/b)^a)*HR_fun(x)},0,t)$value)}
  }
  gen_t <- function(y){stats::uniroot(function(x){Hf(x)-y},interval = c(0,upInt),extendInt="yes")$root}
  set.seed(seed*10)
  U1 <- stats::runif(n1)
  T_1 <-lapply(U1, gen_t) %>%unlist%>%as.vector()
  Tm<-c(T_0,T_1)
  ## in extreme case, time is 0, add 0.1;
  #****************************
  #* Simulate drop-out Time
  #****************************
  Tm <- ifelse(Tm==0,0.1,Tm)
  if (!missing(drop_param)){
    a_drop <- drop_param[1]
    b_drop <- drop_param[2]
    set.seed(seed*99)
    drop_time <- stats::rweibull(t_p1,a_drop,1/b_drop)
  }else{
    drop_time <- rep(1:length(t_p1),Inf)
  }
  #****************************
  #* Simulate entry Time
  #****************************
  #- create the CDF;
  ent_cdf <- function(t){stats::integrate(entry_pdf,lower=0,upper=t)$value}
  gen_ent <- function(y){stats::uniroot(function(x){ent_cdf(x)-y},interval = c(0,upInt),extendInt="yes")$root}

  set.seed(seed+1)
  tu0 <- runif(t_p1)
  t0 <-lapply(tu0, gen_ent) %>%unlist%>%as.vector()
  ot <- t0+Tm

  dat <- data.frame(id=1:t_p1,ent=t0,time=Tm,trt=trt,ot=ot,
                    drop_time=drop_time,ot_drop=t0+drop_time)
  if (type=="event"){
    # find the smallest time between drop-out and event time
    min_ind0 <- apply(cbind(dat$time,dat$drop_time),1,which.min)
    dat$i1 <- min_ind0==1
    dat <- dat[order(dat$ot),]
    dat$c0 <- cumsum(dat$i1)
    Dur <- min(dat[dat$c0==t_p3,]$ot )
    }else{
    Dur <- t_p2+fup # the length of study
  }
  min_ind <- apply(cbind(dat$ot,Dur,dat$ot_drop),1,which.min)
  status <- c("event","end of study","drop-out")
  tot_len <- cbind(dat$ot,Dur,dat$ot_drop)[cbind(seq_along(min_ind),min_ind)]
  dat$t_val <- tot_len-dat$ent
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
  if (print==TRUE){
    s0 <- paste("The specified trial type is",type,"\n")
    s1 <- paste("The specified entry time is",t_p2,"\n")
    s2 <- paste("The max observed time is",Dur,"\n")
    s3 <- paste("Number of events is",sum(final$event),"\n")
    cat(s0)
    cat(s1)
    cat(s2)
    cat(s3)
  }

  list <- list(data=final,
               type=type,
               entrytime=t_p2,
               maxobs=Dur)
  class(list) <- 'SimuTrial'
  return(list)
}
