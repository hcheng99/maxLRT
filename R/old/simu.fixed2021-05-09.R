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
simu.fixed <- function(bsl_dist=c("weibull","loglogistic")
                            ,param
                            # alpha=1 corresponds to exponential
                            ,totN
                            ,fHR #the non proportion hazard function
                            ,prop #proportion of treatment group
                            ,entry =18 #entry time
                            ,fup=24
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
    T_0 <- stats::rweibull(n0,a,1/b)
    #-- get the cummulative hazard function----#
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

  ## for entry time
  set.seed(seed+1)
  t0 <- runif(totN,0,entry)
  ot <- t0+Tm

  dat <- data.frame(id=1:totN,ent=t0,time=Tm,trt=trt,ot=ot)
  dat$t_cnsr <- with(dat,ifelse(ot>entry+fup,0,1))
  dat$t_val <- with(dat,ifelse(ot>entry+fup,entry+fup-ent,Tm))


  return(dat)
}
