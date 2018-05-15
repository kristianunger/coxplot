cp.plot <- function(time, status, strat, col = c("lightseagreen","darkred","blue","purple"), ep, main, baseline = 1, pos.hr = "bottomleft", pos.cols = "bottomright", pos.bas = "topright", lndist = 300){
  
  #######function for comprehensive plotting of cox proportional hazard analysis###
  c.df <- data.frame(time=time, status=status, strat = as.factor(strat))
  
  bl <- levels(c.df$strat)[1]
  
  lv.strat <- levels(c.df$strat)
  
  
  m1 <- coxph(Surv(time, status) ~ strat)
  s.m1 <- summary(m1)
  p.val <- round(s.m1$sctest[3],5)
  
  if(length(lv.strat) == 2) conf.int <- paste(c(round(s.m1$conf.int[3], 2), round(s.m1$conf.int[4], 2)), collapse = "-") else {
  conf.int <- c()
  for(c in 1:(length(lv.strat)-1))
  {
  conf.int.c <- paste(c(round(s.m1$conf.int[c,3], 2), round(s.m1$conf.int[c,4], 2)), collapse = "-")
  conf.int <- c(conf.int,conf.int.c)
  }
  }
  
  if(length(lv.strat)>2) hazard_r <- abs(round(s.m1$coefficients[,2],2)) else hazard_r <- abs(round(s.m1$coefficients[2],2))
  
  max.time <- max(time)
  
  nr.hr <- data.frame()
  for(f in 1:length(lv.strat))
  {
  nr.hr.j <- c()
  for(j in c(seq(0,max.time,500),max.time))
  {
    time.j <- time[strat == lv.strat[f]]
    status.j <- status[strat == lv.strat[f]]
    nr.a <- summary(survfit(Surv(time.j, status.j) ~ 1), j)$n.risk
    if(length(nr.a) > 0) nr.def <- nr.a
    nr.hr.j <- c(nr.hr.j, nr.def)
  }
  nr.hr <- rbind(nr.hr, nr.hr.j)
  }
  
  
  for(k in 1:nrow(nr.hr))
  {
    nr.hr.k <- nr.hr[k,]
    min.nr.hr.k <- which.min(nr.hr.k)
    nr.hr.k[(min.nr.hr.k+1):length(nr.hr.k)] <- ""
    nr.hr[k,] <- nr.hr.k
  }
  
  

  par(oma=c(5,3,2,2))
  par(mar=c(5,3,2,2))
  
  plot(survfit(Surv(time, status) ~ c.df$strat), mark.time = T, lwd = 6, mark="|", cex=1.5, col=col,ylab=ep, xlab="", cex.axis=2, cex.lab=2,xmax=max.time, main=main, cex.main=2, xlim=c(0,3100))
  legend(pos.hr,c(paste("HR:",hazard_r," (95% CI ",conf.int,")",sep=""),paste("p=",p.val,sep="")),bty = "n",cex=2)
  legend(pos.cols,c(lv.strat),col=col,pch="-", bty="n",cex = 2)
  legend(pos.bas,paste("baseline: ",lv.strat[baseline], sep=""), bty="n",cex = 2)
  for(m in 1:nrow(nr.hr))
  {
  mtext(c(lv.strat[m],nr.hr[m,],""),1,at=c(-lndist,seq(0,max.time,500)[1:ncol(nr.hr)], max.time),padj=1+(2*m),cex=2, outer=F)
  }
  
  
}

