cp.plot <- function(time, status, strat, col = c("lightseagreen","darkred","blue","purple"), ep = "endpoint", main ="", baseline = 2, pos.hr = "bottomleft", pos.cols = "bottomright", pos.bas = "topright", lndist = 300, intv = 500, roundfac = 5, max.time.add = 200, med.surv = T, med.surv.pos = 0.1, med.surv.x = 0, cex.med.surv = 1.5, med.surv.round = 2, cex.legend = 1.5, out.file = T, out.file.name = "coxph.out.txt"){
  
  #######function for comprehensive plotting of cox proportional hazard analysis###
  c.df <- data.frame(time=time, status=status, strat = as.factor(strat))
  
  strat.t <- c.df$strat
  
  bl <- levels(strat.t)[baseline]
  
  
  strat.t <- factor(strat.t, levels = unique(c(as.character(strat.t[strat.t==bl]),as.character(strat.t[!strat.t==bl]))))
  
  c.df$strat <- strat.t
  
  lv.strat <- levels(c.df$strat)
  
  m1 <- coxph(Surv(time, status) ~ c.df$strat)
  s.m1 <- summary(m1)
  p.val <- round(s.m1$sctest[3],roundfac)
  
  
  fit <- summary(survfit(Surv(time, status) ~ c.df$strat))
  
  strats <- lv.strat[!lv.strat==bl]
  
  if(length(lv.strat) == 2) conf.int <- paste(c(round(s.m1$conf.int[3], 2), round(s.m1$conf.int[4], 2)), collapse = "-") else {
  conf.int <- c()
    for(c in 1:(length(lv.strat)-1))
  {
  conf.int.c <- paste(c(round(s.m1$conf.int[c,3], 2), round(s.m1$conf.int[c,4], 2)), collapse = "-")
  conf.int <- c(conf.int,conf.int.c)
  
  }
  }
  
  pvs <- c()
  for(s in 1:length(strats))
  {
    strat.s <- strats[s]
    c.s <- coxph(Surv(time[c.df$strat%in%c(bl, strat.s)], status[c.df$strat%in%c(bl, strat.s)]) ~ c.df$strat[c.df$strat%in%c(bl, strat.s)])
    s.c.s <- summary(c.s)
    pvs <- c(pvs, round(s.c.s$sctest[3], roundfac))
  }
  
  ######calculate all possible log-ranks
  cbs <- combn(lv.strat, m = 2)
  
  conf.ints <- c()
  hzds <- c()
  pvds <- c()
  for(d in 1:ncol(cbs))
  {
    sts <- factor(as.character(c.df$strat[c.df$strat%in%cbs[,d]]), levels = cbs[,d])
    cp.d <- coxph(Surv(time[c.df$strat%in%cbs[,d]], status[c.df$strat%in%cbs[,d]]) ~ sts)
    s.dp.d <- summary(cp.d)
    conf.int.d <- paste(c(round(s.dp.d$conf.int[3], 2), round(s.dp.d$conf.int[4], 2)), collapse = "-")
    conf.ints <- c(conf.ints, conf.int.d)
    hzds <- c(hzds, round(s.dp.d$conf.int[1], 2))
    pvds <- c(pvds, round(s.dp.d$sctest[3], roundfac))
  }
  
  out.tab <- data.frame(comparison = paste(cbs[1,],cbs[2,], sep ="-"), p.value = pvds, hazard.ratio = hzds, CI95 = conf.ints, direction = rep("flip", length(pvds)))
  
  
  cbs.rev <- rbind(cbs[2,],cbs[1,])
  conf.ints.rev <- c()
  hzds.rev <- c()
  pvds.rev <- c()
  for(d in 1:ncol(cbs.rev))
  {
    sts.rev <- factor(as.character(c.df$strat[c.df$strat%in%cbs.rev[,d]]), levels = cbs.rev[,d])
    cp.d.rev <- coxph(Surv(time[c.df$strat%in%cbs.rev[,d]], status[c.df$strat%in%cbs.rev[,d]]) ~ sts.rev)
    s.dp.d.rev <- summary(cp.d.rev)
    conf.int.d.rev <- paste(c(round(s.dp.d.rev$conf.int[3], 2), round(s.dp.d.rev$conf.int[4], 2)), collapse = "-")
    conf.ints.rev <- c(conf.ints.rev, conf.int.d.rev)
    hzds.rev <- c(hzds.rev, round(s.dp.d.rev$conf.int[1], 2))
    pvds.rev <- c(pvds.rev, round(s.dp.d.rev$sctest[3], roundfac))
  }
  
  out.tab.rev <- data.frame(comparison = paste(cbs.rev[1,],cbs.rev[2,], sep ="-"), p.value = pvds.rev, hazard.ratio = hzds.rev, CI95 = conf.ints.rev, direction = rep("flop", length(pvds)))
  
  out.tab <- rbind(out.tab, out.tab.rev)
  
  
  if(out.file) write.table(out.tab, file = out.file.name, sep = "\t", row.names = F, quote = F)
  
  if(length(lv.strat)>2) hazard_r <- abs(round(s.m1$coefficients[,2],2)) else hazard_r <- abs(round(s.m1$coefficients[2],2))
  
  max.time <- max(time) + max.time.add
  
  nr.hr <- matrix("", ncol = length(c(seq(0,max.time,intv),max.time)), nrow = length(lv.strat))
  for(f in 1:length(lv.strat))
  {
  nr.hr.j <- rep("", length(c(seq(0,max.time,intv),max.time)) )
  it = 0
  for(j in c(seq(0,max.time,intv),max.time))
  {
    it = it +1
    time.j <- time[strat == lv.strat[f]]
    status.j <- status[strat == lv.strat[f]]
    nr.a <- summary(survfit(Surv(time.j, status.j) ~ 1), j)$n.risk
    if(length(nr.a) > 0) nr.def <- nr.a else nr.def <- ""
    nr.hr.j[it] <- nr.def
    nr.hr.j <- as.character(nr.hr.j)
  }
  nr.hr[f,] <- nr.hr.j
  }
  
  
  
  plot(survfit(Surv(time, status) ~ c.df$strat), mark.time = T, lwd = 6, mark="|", cex=1.5, col=col,ylab=ep, xlab="", cex.axis=2, cex.lab=2,xmax=max.time, main=main, cex.main=2, xlim=c(0,max.time))
  legend(pos.hr,c(paste(strats,"-",bl," HR:",hazard_r," (95% CI ",conf.int,"), p: ",pvs,sep=""),paste("p=",p.val,sep="")),bty = "n",cex=cex.legend)
  legend(pos.cols,c(lv.strat),col=col,pch="-", bty="n",cex = cex.legend)
  legend(pos.bas,paste("baseline: ",lv.strat[1], sep=""), bty="n",cex = cex.legend)
  for(m in 1:nrow(nr.hr))
  {
  mtext(c(lv.strat[m],nr.hr[m,],""),1,at=c(-lndist,seq(0,max.time,intv), max.time),padj=1+(2*m),cex=2, outer=F)
  }
  
  if(med.surv){
  
    abline(h = 0.5, lty = 2)
    abline(v = fit$table[,7], lty = 2)
    text(fit$table[,7]+med.surv.x, y = med.surv.pos, round(fit$table[,7],med.surv.round), cex = cex.med.surv)
    lines(survfit(Surv(time, status) ~ c.df$strat), mark.time = T, lwd = 6, mark="|", cex=1.5, col=col,ylab=ep, xlab="", cex.axis=2, cex.lab=2,xmax=max.time, main=main, cex.main=2, xlim=c(0,max.time))
  }
}

