
#*********************************************

###input parameters
pwr2n.maxLR<- function(entry   = 1
                     ,fup      = 1
                     ,k        = 100
                     ,trans.prob
                     ,hazR
                     ,Wlist
                     ,ratio    = 1
                     ,alpha    = 0.05
                     ,beta     = 0.1
                     ,two.side = TRUE
                     ,plot     = TRUE
                     ,nocensor = FALSE
                     ,criteria = 100
){
  tot_time <- entry+fup
  num <- k*tot_time

  x <- seq(0.01,tot_time,length=k)
  haz_val <- hazR(x)*trans.prob[5]
  haz_point <- x*k
  hazard <- stats::stepfun(haz_point,c(haz_val,haz_val[length(haz_val)]),
                    right=TRUE)
  ##initial states for treatment
  D1 <- c(0,0,1,0,0)
  ##initial states for control
  D0 <- c(0,0,0,1,0)
  temp <- diag(rep(1,4))
  pdat <- matrix(NA,nrow=num,ncol=21)
  ti <-0
  l <- 1
  L_trans <- list()
  for (i in 1:num){
    j <- ceiling(i/k)
    trans.prob[2] <- hazard(i)
    # before a, there is no censor
    if (nocensor==FALSE){
      if (j <=fup){
        L_trans[[i]] <- cbind(rbind(CrtTM(Plist = trans.prob,K=k),rep(0,4)),
                              c(0,0,0,0,1))

      }else {
        L_trans [[i]] <- CrtTM_C(Plist = trans.prob,K=k,a=entry,b=fup,i=i)

      }
    }else {
      L_trans[[i]] <- cbind(rbind(CrtTM(Plist = trans.prob,K=k),rep(0,4)),
                            c(0,0,0,0,1))
    }

    ti <- ti+1/k
    if (i==1){temp1 <- L_trans[[i]]%*%D1 ; temp0 <- L_trans[[i]]%*%D0
    } else {temp1 <- L_trans[[i]]%*%temp1; temp0 <- L_trans[[i]]%*%temp0;}
    pdat[i,1:11] <- c(ti,t(temp1),t(temp0))
    if (i==1) {
      #phi
      pdat[i,12] <- 1
      ## hazard of dying for trt=1
      pdat[i,13] <- pdat[i,3]
      # for trt=0
      pdat[i,14] <- pdat[i,8]

    }
    else {
      #phi
      pdat[i,12] <- sum(pdat[i-1,5],pdat[i-1,4])/sum(pdat[i-1,10],pdat[i-1,9])
      # f0 <- function(x) {-exp(-pdat[i,1]*x)+exp(-pdat[i-1,1]*x)-pdat[i,3]+pdat[i-1,3]}
      # #f0 <- function(x) {-exp(-pdat[i,1]*x)+exp(-pdat[i-1,1]*x)+exp(-pdat[i,1])-exp(-pdat[i-1,1])}
      # f1 <- function(x) {-exp(-pdat[i,1]*x)+exp(-pdat[i-1,1]*x)-pdat[i,7]+pdat[i-1,7]}
      # pdat[i,11] <- uniroot(f0,interval = c(0.8,1),extendInt = "yes")$root
      # pdat[i,12] <- uniroot(f1,interval = c(0,0.5),extendInt = "yes")$root
      #hazard1
      pdat[i,13] <- (pdat[i,3]-pdat[i-1,3])/(pdat[i-1,5]+pdat[i-1,4])
      #hazard 0
      pdat[i,14] <- (pdat[i,8]-pdat[i-1,8])/(pdat[i-1,9]+pdat[i-1,10])
      #pdat[i,11] <- 1-(pdat[i,4]+pdat[i,5])/(pdat[i-1,5]+pdat[i-1,4])
      #pdat[i,12] <- 1-(pdat[i,8]+pdat[i,9])/(pdat[i-1,9]+pdat[i-1,8])
      # pdat[i,11]=1
      # pdat[i,12]=0.5
    }
    #theta
    pdat[i,15] <- pdat[i,13]/pdat[i,14]
    # pdat[i,13] <- 2
    #gamma
    pdat[i,16] <- pdat[i,12]*pdat[i,15]*ratio/(1+pdat[i,12]*pdat[i,15]*ratio)-
      pdat[i,12]*ratio/(1+pdat[i,12]*ratio)
    #eta
    pdat[i,17] <- pdat[i,12]*ratio/(1+pdat[i,12]*ratio)^2
    # print(pdat[i,])
  }
  #pdat
  rho <- pdat[num,3]+pdat[num,8]

  ## sample size
  for (i in 1:num){
    ## S1
    #pdat[i,20] <- sum(pdat[i,5],pdat[i,4])
    #pdat[i,20] <- exp(-sum(pdat[1:i,13]))
    pdat[i,20] <- prod(1-pdat[1:i,13])
    ## S0
    # pdat[i,21] <- sum(pdat[i,10],pdat[i,9])
    #pdat[i,21] <- exp(-sum(pdat[1:i,14]))
    pdat[i,21] <- prod(1-pdat[1:i,14])
    ## mean survival
    # pdat[i,19] <- mean(sum(pdat[i,5],pdat[i,4]),sum(pdat[i,10],pdat[i,9]))
    pdat[i,19] <- mean(c(pdat[i,20],pdat[i,21]))
    ## rho
    if (i==1){
      pdat[i,18] <- (pdat[i,3]+pdat[i,8])/rho
    }else
    {
      pdat[i,18] <- (pdat[i,3]-pdat[i-1,3]+pdat[i,8]-pdat[i-1,8])/rho

    }
  }
  pdat <- as.data.frame(pdat)
  names(pdat) <-
    c("ti","E_L","E_E","E_Ae","E_Ac","E_C","C_L","C_E","C_Ae",
      "C_Ac","C_C","phi","hazard_E","hazard_C","theta","gamma","eta",
      "rho","S","S1","S0")

  # w_f <- function(x){2*x-1}
  # w_f <- function(x){1}
  # wgt <- wgt_fun(pdat$S)
  ## number of weight functions

  wn <- length(Wlist)
  W <- matrix(NA,nrow=nrow(pdat),ncol=wn)
  ## calculate the variance-covariance matrix
  Vmat <- matrix(NA,nrow=wn,ncol=wn)
  event <- c()
  if (two.side==TRUE){
    part2=(qnorm(1-alpha/2)+qnorm(1-beta))^2
  }
  for (j in 1:wn){
    W[,j] <- Wlist[[j]](pdat$S)
    wgt <- W[,j]
    part1=t(pdat$rho)%*%(pdat$eta*wgt^2)
    part3=t(pdat$rho)%*%(pdat$gamma*wgt)
    part3=part3^2
    event[j] <- part1*part2/part3

  }
  #print(c(part1,part2,part3))
  ## get the min and max sample size
  E_min=min(event)
  E_max=max(event)
  #print(event)
  eta=-1
  dnum <- E_min
  count <- 0
  while(eta <0){
    count <- count+1
    for (k1 in 1:wn){
      for (k2 in 1:wn){
        Vmat[k1,k2] <-  dnum*t(W[,k1]*W[,k2]) %*%(pdat$rho*pdat$eta)
      }
    }
    rho_est <- cov2cor(Vmat)
    mu <- as.vector(dnum*t(W)%*%(pdat$rho*pdat$gamma))/sqrt(diag(Vmat))
    crit <- mvtnorm::qmvnorm(1-alpha,tail="both.tails",
                    mean=rep(0,wn),sigma = rho_est)$quantile
    #print(crit)
    # power <- 1-pmvnorm(-Inf,crit,mean=mu,sigma = rho_est)+
    #   pmvnorm(-Inf,-crit,mean=mu,sigma = rho_est)
    power <- 1-mvtnorm::pmvnorm(-crit,crit,mean=mu,sigma = rho_est)
    #print(crit)
    #print(mu)
    #print(power)
    if (power<1-beta&count<criteria) {dnum<-dnum+1}
    else if (count>=criteria){
      warning(paste0("the algoritm doesn't converge within ",criteria," iterations;
                     the current power is ",power,"; event size: ",dnum))
      ;break}
    else {break}
  }




  Nsize <- dnum/mean(c(pdat$C_E[num],pdat$E_E[num]))
  #print(Nsize)
  #print(mean(c(pdat$C_E[num],pdat$E_E[num])))
  #print(num)
  if (plot==TRUE){
    Se1 <- exp(-trans.prob[2]*pdat$ti)
    Se0 <- exp(-trans.prob[5]*pdat$ti)

    graphics::par(mfrow=c(1,1))
    with(pdat,
         plot(ti,theta,cex=0.1,ylab="hazard ratio"))

    with(pdat,{
      plot(ti,S1,cex=0.1,lty=1,ylim=c(0,1),col=1)
      graphics::lines(ti,S0,col=2,lty=2)

    })
    #lines(pdat$ti,Se1,col=3,lty=1)
    #lines(pdat$ti,Se0,col=4,lty=2)
    graphics::par(mar=c(5.1, 4.1, 8, 4.1), xpd=TRUE)
    #legend("top",inset=c(0.1,-0.25),legend = c("Exp-S1", "Exp-S0","S1","S0"), lty=c(1,2,1,2), lwd = 1,
    #       col=c(1,2,3,4), horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
    graphics::legend("top",inset=c(0.1,-0.25),legend = c("S1","S0"), lty=c(1,2,1,2), lwd = 1,
           col=c(1,2,3,4), horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
    ## loss to follow-up
    with(pdat,{
      plot(ti,E_L,cex=0.1,lty=1,col=1,ylim=c(0,max(c(E_L,C_L))))
      graphics::lines(ti,C_L,lty=2,col=2)
    })
  }


  return(list( eventN  = dnum
               ,totalN = Nsize
               ,pwr = as.numeric(power)
               ,prob_event = mean(c(pdat$C_E[num],pdat$E_E[num]))
               ,L_trans = L_trans
               ,pdat = pdat

  )
  )

}

#*********************************************
