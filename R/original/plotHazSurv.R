
plotHazSurv <- function(
  bsl_dist=c("weibull","loglogistic")
  ,param=c(1.2,0.03)
  ,fun_list=fun_list
  ,end
  ,tit= c("Hazard Function","Survival Function")
  ,pos=c(1,2)
  ,f_n=1

){
  t <- seq(0.5,end,length=300)
  fn <- length(tit)

  a <- param[1]
  b <- param[2]
  graphics::par(mar=c(4,3,3,1),mfrow=pos)
  #par(mfrow=pos)
  new_fun <- fun_list[f_n]
  print(new_fun)
  for (i in 1:length(new_fun)){
    ##-----------------draw the hazard function------------------##
    if (bsl_dist=="weibull"){
      hr1 <- a*b*(b*t)^(a-1)

    }else if (bsl_dist=="loglogistic"){
      hr1 <- a/b*(t/b)^(a-1)/(1+(t/b)^a)
    }
    Hft <- function(t){
      if (t==0) {1
      }else {Hf(t)}
    }
    hr <- new_fun[[i]](t)*hr1
    plot(t,hr,cex=0.3,ylab="h(t)",cex.main=1,
         ylim = c(0,max(hr,hr1)*1.1),type="l",lty=2,col="blue",xlab="Time")
    graphics::lines(t,hr1, lty = 1,,col="red")
    dis <-tit[1]
    graphics::mtext(dis, side=3,line=1,adj=0,,cex=0.8)
  graphics::par(xpd=TRUE)
    graphics::legend("topright",legend = c("Treatment Group", "Control Group"), lty=c(1,2), lwd = 1,
           col=c("red","blue"), xpd = TRUE, horiz = FALSE, cex = 0.5)


  }


  for (i in 1:length(new_fun)){
    ##-----------------draw the survival function------------------##
    if (bsl_dist=="weibull"){
      S1 <- exp(-(b*t)^a)
      Hf <- function(t){exp(-1*a*b*stats::integrate( function(x){(x*b)^(a-1)*new_fun[[i]](x)},0,t)$value)}
    }else if (bsl_dist=="loglogistic"){
      S1 <- b^a/(b^a+t^a)
      Hf <- function(t){exp(-1* a/b*stats::integrate( function(x){(x/b)^{a-1}/(1+(x/b)^a)*new_fun[[i]](x)},0,t)$value)}
    }
    Hft <- function(t){
      if (t==0) {1
      }else {Hf(t)}
    }
    S0 <- as.vector(unlist(lapply(t,Hft)))
    plot(t,S0,cex=0.3,cex.main=1,ylab = "S(t)",type="l",
         ylim = c(0,1),lty=2,col="blue",xlab="Time")
    graphics::lines(t,S1,lty=1,col="red")
    dis <-tit[2]
    graphics::mtext(dis, side=3,line=1,adj=0,,cex=0.8)
    graphics::legend("topright",legend = c("Treatment Group", "Control Group"), lty=c(1,2), lwd = 1,
           col=c("red","blue"), xpd = TRUE, horiz = FALSE, cex = 0.5)


  }
  graphics::par(mfrow=c(1,1))

}
