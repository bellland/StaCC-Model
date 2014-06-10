##import BLA

import.bla <- function(){  
  ff <- list.files()
    ff <- ff[grep('BLA',ff)]
	ff <- ff[grep('csv',ff)]

  BLA <- list()	
  spBLA <- list()
  for(yy in 1:length(ff)){
    BLA[[yy]]   <- read.csv(ff[yy],header=T)
	spBLA[[yy]] <- substr(BLA[[yy]][,1],1,4)
	}
  
  xl <- expression(italic(G[ref]))
  yl <- expression(lambda)
jpeg('blaplot.jpg',width=6,height=6,unit='in',res=1000)
  par(mfrow=c(2,2))
  for(j in 1:4){
    plot(BLA[[j]][,2],-BLA[[j]][,3],pch=PT[match(spBLA[[j]],spall)],
	  main = yrall[j],xlab=xl,ylab=yl)
    points(tapply(BLA[[j]][,2],spBLA[[j]],mean),-tapply(BLA[[j]][,3],spBLA[[j]],mean),
	  pch=PT[match(names(tapply(BLA[[j]][,2],spBLA[[j]],mean)),spall)],lwd=2)
    abline(0,.6)
  }  
  legend('topleft',legend=spname,pch=PT,cex=.8,bty='n')
dev.off()
list(BLA = BLA, spBLA = spBLA)
}

##data figure file names
data.figures <- function(){
  ff <- list.files()
    f.vpd <- ff[grep('VPD',ff)]
    f.sm  <- ff[grep('sm',ff)]
    f.par <- ff[grep('PAR',ff)]
    f.lai <- ff[grep('LAI',ff)]

    yr.vpd <- as.integer(substr(f.vpd,4,7))
    yr.sm  <- as.integer(substr(f.sm,3,6))
    yr.par <- as.integer(substr(f.par,4,7))
    yr.lai <- as.integer(substr(f.lai,4,7))

## import data
  for(j in 1:length(yr.vpd)){
    tmp <- read.csv(f.vpd[j])
	if(j==1) vpd <- tmp
	if(j>1)  vpd <- rbind(vpd,tmp)
    }

  for(j in 1:length(yr.lai)){
    tmp <- read.csv(f.lai[j])
	if(j==1) lai <- tmp
	if(j>1)  lai <- rbind(lai,tmp)
    }

  for(j in 1:length(yr.sm)){
    tmp <- read.csv(f.sm[j])
	if(j==1) sm <- tmp
	if(j>1)  sm <- rbind(sm,tmp)
    }

  for(j in 1:length(yr.par)){
    tmp <- read.csv(f.par[j])
	if(j==1) par <- tmp
	if(j>1)  par <- rbind(par,tmp)
    }

#rescale time
t.vpd <- vpd[,'Time']
t.vpd[t.vpd>30] <- trunc(t.vpd[t.vpd>30]/100) + 
  (t.vpd[t.vpd>30] - trunc(t.vpd[t.vpd>30]/100)*100)/30
t.vpd <- vpd[,'year'] + 
         (vpd[,'DOY']-1)/(365+as.integer(vpd[,'year']==2004)) +
		 t.vpd/24/(365+as.integer(vpd[,'year']==2004))
  
t.lai <- lai[,'Time']
t.lai[t.lai>30] <- trunc(t.lai[t.lai>30]/100) + 
  (t.lai[t.lai>30] - trunc(t.lai[t.lai>30]/100)*100)/30
t.lai <- lai[,'year'] + 
         (lai[,'DOY']-1)/(365+as.integer(lai[,'year']==2004)) +
		 t.lai/24/(365+as.integer(lai[,'year']==2004))

t.par <- par[,'Time']
t.par[t.par>30] <- trunc(t.par[t.par>30]/100) + 
  (t.par[t.par>30] - trunc(t.par[t.par>30]/100)*100)/30
t.par <- par[,'year'] + 
         (par[,'DOY']-1)/(365+as.integer(par[,'year']==2004)) +
		 t.par/24/(365+as.integer(par[,'year']==2004))

t.sm <- sm[,'Time']
t.sm[t.sm>30] <- trunc(t.sm[t.sm>30]/100) + 
  (t.sm[t.sm>30] - trunc(t.sm[t.sm>30]/100)*100)/30
t.sm   <- sm[,'year'] + 
         (sm[,'DOY']-1)/(365+as.integer(sm[,'year']==2004)) +
		 t.sm/24/(365+as.integer(sm[,'year']==2004))

#plot data
jpeg('data.plot.jpg',width=2.8,height=5,units='in',res=1000)
sp.combo <- F
  par(mfrow=c(3,1),mar=c(3,5,1,1))
    plot(t.vpd,vpd[,5],type='l',xlab="",ylab=expression(italic(D[t])~~(kPa)),cex.lab=1.2)
	  rect(2001,-1,2002.3,10,density=10,col='black')
	  rect(2002.8,-1,2003.3,10,density=10,col='black')
	  rect(2003.8,-1,2004.3,10,density=10,col='black')
	  rect(2004.8,-1,2005.3,10,density=10,col='black')
	  rect(2005.8,-1,2006.3,10,density=10,col='black')
	  rect(2005.8,-1,2006.3,10,density=10,col='black')

    plot(t.sm,sm[,5],type='l',xlab="",ylab=expression(italic(M[t])~~(mm^3~~mm^-3)),cex.lab=1.2)
	  rect(2001,-1,2002.3,10,density=10,col='black')
	  rect(2002.8,-1,2003.3,10,density=10,col='black')
	  rect(2003.8,-1,2004.3,10,density=10,col='black')
	  rect(2004.8,-1,2005.3,10,density=10,col='black')
	  rect(2005.8,-1,2006.3,10,density=10,col='black')

    if(sp.combo) plot(t.lai,apply(lai[,5:10],1,sum),type='l',
	  xlab="",ylab=expression(summed~~italic(A[L])),cex.lab=1.2)
    if(!sp.combo) {
	  plot(t.lai,lai[,5],type='l',xlab="",ylab=expression(summed~~italic(A[L])),cex.lab=1.2,ylim=c(0,3))
      for(j in 1:6) {
        pkeep <- trunc(seq(0,length(t.lai),length=24))
		points(t.lai[pkeep],lai[pkeep,4+j],pch=j)
	    lines(t.lai,lai[,4+j])
        }
	  }
	  rect(2001,-1,2002.3,10,density=10,col='black')
	  rect(2002.8,-1,2003.3,10,density=10,col='black')
	  rect(2003.8,-1,2004.3,10,density=10,col='black')
	  rect(2004.8,-1,2005.3,10,density=10,col='black')
	  rect(2005.8,-1,2006.3,10,density=10,col='black')
	  
	  legend(2001.9,3,legend=spabbr,pch=1:6,cex=.70,horiz=TRUE,x.intersp = .5,y.intersp=.2)
	  
#    plot(t.par,par[,5],type='l')
dev.off()

}

#plot convergences by species, year, and kappa

converge <- function(thin){
  for(j in 1:dim(spyr)[1]){
    postscript(paste('converge',spyr[j,1],spyr[j,2],'ps',sep='.'))
	  tmp <- which(splist == spyr[j,1] & yrlist == spyr[j,2])
	  if(length(tmp)==0) next
	  for(i in 1:length(tmp)){
      gtrim <- seq(burnin,g.all[tmp[i]],by=thin)
	  par(mfrow=c(3,3))
        plot(gtrim,ggibbs.all[[tmp[i]]][gtrim,1],type='l',
		  ylab=expression(italic(G[ref])),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,agibbs.all[[tmp[i]]][gtrim,1],type='l',
		  ylab=expression(lambda),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,lgibbs.all[[tmp[i]]][gtrim,1],type='l',
		  ylab=expression(italic(L)[1]),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,lgibbs.all[[tmp[i]]][gtrim,2],type='l',
		  ylab=expression(italic(L)[2]),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,mgibbs.all[[tmp[i]]][gtrim,1],type='l',
		  ylab=expression(italic(M[1])),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,mgibbs.all[[tmp[i]]][gtrim,2],type='l',
		  ylab=expression(italic(M[2])),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,vgibbs.all[[tmp[i]]][gtrim,1],type='l',
		  ylab=expression(sigma^2),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,vgibbs.all[[tmp[i]]][gtrim,2],type='l',
		  ylab=expression(italic(V[J])),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
        plot(gtrim,vgibbs.all[[tmp[i]]][gtrim,3],type='l',
		  ylab=expression(italic(V[a])),
		  main=paste(splist[tmp[i]],yrlist[tmp[i]],kplist[tmp[i]]))
	  }
    dev.off()  
  }
}  

##
hist.exp <- function(){
  jpeg('hist.exp.jpg',width=169/25.4,height=3,units='in',res=1000)
  par(mfrow=c(2,3),mar=c(5,5,1,1))
  for(stmp in spall){
	    hist(as.numeric(r2[r2[,'species']==stmp,5]),main=spname[which(spall==stmp)],
		  xlab="",ylab="",cex.lab=1.4,cex.axis=1.2)
    }
    
  mtext(expression(italic(r)^2),side=1,outer=TRUE,line=-1.5)
  mtext("No. of Probe-Sets",side=2,outer=TRUE,line=-1.5)
  
  dev.off()
  }

  ## pred vs. observed -- TC
pred.vs.obs.TC <- function(){

r2 <- matrix(NA,0,length(kpall)+3)
  colnames(r2) <- c('species','year','probe',kpall)

pdf('pred.obs.pdf')
par(mfrow=c(3,3),mar=c(5,5,2,1))
for(stmp in spall){
  for(ytmp in yrall){
    tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
	if(length(tmp)==0) next
    for(n in 1:dim(Jdata.all[[tmp[1]]])[2]){
	  tp <- rep(NA,length(kpall))
	  for(j in 1:length(kpall)){
        tvec <- sample(dim(Jdata.all[[tmp[j]]])[1],1000)
	    td <- cbind(Jdata.all[[tmp[j]]][tvec,n],Jpred.all[[tmp[j]]][tvec,n])
        plot(td,xlab='observed',ylab='predicted',main=paste(stmp,ytmp,n,kpall[j]))
          abline(0,1,lty=2)
        text(min(td[,1],na.rm=T)+(max(td[,1],na.rm=T)-min(td[,1],na.rm=T))*.05,
          max(td[,2],na.rm=T),expression(italic(r)^2))
        tp[j] <- cor(td,use='complete.obs')[1,2]^2
		text(min(td[,1],na.rm=T)+(max(td[,1],na.rm=T)-min(td[,1],na.rm=T))*.25,
          max(td[,2],na.rm=T),labels=paste("=",round(as.numeric(tp[j]),2)))
        }
        r2 <- rbind(r2,c(stmp,ytmp,n,tp))
	  }
    }
  }
dev.off()

r2
}

##another pred vs obs

exp.var.sp <- function(){
pdf('explained.variation.sp.year.pdf')
for(stmp in spall){
  par(mfrow=c(2,2))
  for(ytmp in yrall){
    np <- r2[r2[,'species']==stmp & r2[,'year']==ytmp,3]
    l<-1
	for(n in 1:length(np)){
	  keep <- which(r2[,'species']==stmp & r2[,'year']==ytmp & r2[,'probe']==np[n])
	  if(length(keep)==0) next
	  if(l==1) plot(trunc(kpall*30),as.numeric(r2[keep,3+1:length(kpall)]),
	             type='l',ylim=c(.25,1),main=paste(stmp,ytmp),
				 ylab=expression(italic(r)^2),xlab='time constant')
	  l<-0
	  lines(trunc(30*kpall),r2[keep,3+1:length(kpall)])
	    maxk <- which(r2[keep,3+1:length(kpall)]==max(r2[keep,3+1:length(kpall)]))
	  points(trunc(30*kpall)[maxk],r2[keep,3+1:length(kpall)][maxk],pch=1,cex=2)
	    maxk <- which(round(as.numeric(r2[keep,3+1:length(kpall)]),2)==round(max(as.numeric(r2[keep,3+1:length(kpall)])),2))
	  points(trunc(30*kpall)[maxk],r2[keep,3+1:length(kpall)][maxk],pch=20)
	  
	  }
    }
  }
  dev.off()
  }
  
## pred vs. observed -- species

#explained variation
exp.var.ind <- function(){
postscript('explained.variation.eps',width=6,height=4,onefile=FALSE,horizontal=FALSE)

let <- c('(a)','(b)','(c)','(d)','(e)','(f)')

par(mfrow=c(2,3),mar=c(5,5,.5,.5))
  for(ss in 1:6){
	stmp <- spall[ss]
	l<-1
	for(yy in 1:4){
	  ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
	    plot(1-exp(-1/kplist[tmp]),R2[tmp],ylab="",xlab='',
	         pch=20,ylim=c(.25,1))
	  points(1-exp(-1/kplist[tmp]),R2[tmp],col=yy,pch=20)
	  lines(1-exp(-1/kplist[tmp]),R2[tmp],col=yy)
	  
	  l <- 0
	}
	legend('topleft',legend=let[ss],bty="n")
	if(ss==2)   legend('bottomleft',legend=yrall,text.col=1:4,cex=1.1,bty='n')

  }
  
  mtext(expression(Squared~~Pearson~~Correlation~~italic(r)^2),side=2,outer=TRUE,line=-1.7)
  mtext(expression(Capacitance~~Parameter~~italic(beta)),side=1,outer=TRUE,line=-1.7)
  dev.off()
}
	

	#error biplot
err.biplot <- function(){
	postscript('err.biplot.eps',width=6,height=4,onefile=FALSE,horizontal=FALSE)
  par(mfrow=c(2,3),mar=c(5,5,.5,.5))
  for(ss in 1:6){  
    stmp <- spall[ss]
	l <- 1 
    xl <- ""
	if(stmp=='qumi') xl <- expression(sqrt(italic(V[J])))
	if(stmp %in% c('cato')) par(mar=c(4,5,2,0))
	if(stmp %in% c('list')) par(mar=c(4,3,2,2))
	if(stmp %in% c('litu')) par(mar=c(4,1,2,4))
	if(stmp %in% c('qual')) par(mar=c(4,5,2,0))
	if(stmp %in% c('qumi')) par(mar=c(4,3,2,2))
	if(stmp %in% c('quph')) par(mar=c(4,1,2,4))
	
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),xlab="",main=spname[ss],
	      ylab="",pch=20,xlim=range(sqrt(vjvec)),
		  ylim=range(sqrt(sigvec)),log='xy',cex.lab=1.4)
      
	    points(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),pch=14+yy,col=yy,cex=1.5)
	    #text(vjvec[tmp],sigvec[tmp],labels=round(kplist[tmp]*30,0),pos=1)
	    arrows(sqrt(vjvec[tmp[-length(tmp)]]),sqrt(sigvec[tmp[-length(tmp)]]),
	           sqrt(vjvec[tmp[-1]]),sqrt(sigvec[tmp[-1]]),length=.075,col=yy)
	  l<-0
		  }
    
		if(ss==6){
		  mtext(expression(sigma),side=2,outer=TRUE,line=-1.7)
		  mtext(expression(sqrt(italic(S))),side=1,outer=TRUE,line=-1.7)
		}
		if(stmp=='litu')
	  legend('topright',legend=yrall,col=1:4,pch=15:18,text.col=1:4,bty='n',cex=1)
	}
	dev.off()


	postscript('G.biplot.eps',width=6,height=4,onefile=FALSE,horizontal=FALSE)
  par(mfrow=c(2,3),mar=c(5,5,.5,.5))
  for(ss in 1:6){  
    stmp <- spall[ss]
	l <- 1 
    xl <- ""
	if(stmp=='qumi') xl <- expression(sqrt(italic(V[J])))
	if(stmp %in% c('cato')) par(mar=c(4,5,2,0))
	if(stmp %in% c('list')) par(mar=c(4,3,2,2))
	if(stmp %in% c('litu')) par(mar=c(4,1,2,4))
	if(stmp %in% c('qual')) par(mar=c(4,5,2,0))
	if(stmp %in% c('qumi')) par(mar=c(4,3,2,2))
	if(stmp %in% c('quph')) par(mar=c(4,1,2,4))
	
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(gvec[tmp],agvec[tmp],main=spname[ss],
	      ylab="",xlab="",pch=20,xlim=range(gvec),
		  ylim=range(agvec),log='xy',cex.lab=1.4)
      
	    points(gvec[tmp],agvec[tmp],pch=14+yy,col=yy,cex=1.5)
	    #text(vjvec[tmp],sigvec[tmp],labels=round(kplist[tmp]*30,0),pos=1)
	    arrows(gvec[tmp[-length(tmp)]],agvec[tmp[-length(tmp)]],
	           gvec[tmp[-1]],agvec[tmp[-1]],length=.075,col=yy)
	  l<-0
		  }
    
		if(ss==6){
		  mtext(expression(italic(lambda/G[ref])),side=2,outer=TRUE,line=-1.7)
		  mtext(expression(italic(G[ref])),side=1,outer=TRUE,line=-1.7)
		}
		if(stmp=='litu')
	  legend('bottomright',legend=yrall,col=1:4,pch=15:18,text.col=1:4,bty='n',cex=1)
	}
	dev.off()

	postscript('L.biplot.eps',width=6,height=4,onefile=FALSE,horizontal=FALSE)
  par(mfrow=c(2,3),mar=c(5,5,.5,.5))
  for(ss in 1:6){  
    stmp <- spall[ss]
	l <- 1 
    xl <- ""
	if(stmp=='qumi') xl <- expression(sqrt(italic(V[J])))
	if(stmp %in% c('cato')) par(mar=c(4,5,2,0))
	if(stmp %in% c('list')) par(mar=c(4,3,2,2))
	if(stmp %in% c('litu')) par(mar=c(4,1,2,4))
	if(stmp %in% c('qual')) par(mar=c(4,5,2,0))
	if(stmp %in% c('qumi')) par(mar=c(4,3,2,2))
	if(stmp %in% c('quph')) par(mar=c(4,1,2,4))
	
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(l1vec[tmp],l2vec[tmp],main=spname[ss],
	      ylab="",xlab="",pch=20,xlim=range(l1vec),
		  ylim=range(l2vec),log='xy',cex.lab=1.4)
      
	    points(l1vec[tmp],l2vec[tmp],pch=14+yy,col=yy,cex=1.5)
	    #text(vjvec[tmp],sigvec[tmp],labels=round(kplist[tmp]*30,0),pos=1)
	    arrows(l1vec[tmp[-length(tmp)]],l2vec[tmp[-length(tmp)]],
	           l1vec[tmp[-1]],l2vec[tmp[-1]],length=.075,col=yy)
	  l<-0
		  }
    
		if(ss==6){
		  mtext(expression(italic(alpha[2])),side=2,outer=TRUE,line=-1.7)
		  mtext(expression(italic(alpha[1])),side=1,outer=TRUE,line=-1.7)
		}
		if(stmp=='litu')
	  legend('bottomright',legend=yrall,col=1:4,pch=15:18,text.col=1:4,bty='n',cex=1)
	}
	dev.off()

	postscript('M.biplot.eps',width=6,height=4,onefile=FALSE,horizontal=FALSE)
  par(mfrow=c(2,3),mar=c(5,5,.5,.5))
  for(ss in 1:6){  
    stmp <- spall[ss]
	l <- 1 
    xl <- ""
	if(stmp=='qumi') xl <- expression(sqrt(italic(V[J])))
	if(stmp %in% c('cato')) par(mar=c(4,5,2,0))
	if(stmp %in% c('list')) par(mar=c(4,3,2,2))
	if(stmp %in% c('litu')) par(mar=c(4,1,2,4))
	if(stmp %in% c('qual')) par(mar=c(4,5,2,0))
	if(stmp %in% c('qumi')) par(mar=c(4,3,2,2))
	if(stmp %in% c('quph')) par(mar=c(4,1,2,4))
	
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(m1vec[tmp],m2vec[tmp],main=spname[ss],
	      ylab="",xlab="",pch=20,xlim=range(m1vec),
		  ylim=range(m2vec),log='xy',cex.lab=1.4)
      
	    points(m1vec[tmp],m2vec[tmp],pch=14+yy,col=yy,cex=1.5)
	    #text(vjvec[tmp],sigvec[tmp],labels=round(kplist[tmp]*30,0),pos=1)
	    arrows(m1vec[tmp[-length(tmp)]],m2vec[tmp[-length(tmp)]],
	           m1vec[tmp[-1]],m2vec[tmp[-1]],length=.075,col=yy)
	  l<-0
		  }
    
		if(ss==6){
		  mtext(expression(italic(alpha[4])),side=2,outer=TRUE,line=-1.7)
		  mtext(expression(italic(alpha[3])),side=1,outer=TRUE,line=-1.7)
		}
		if(stmp=='litu')
	  legend('bottomright',legend=yrall,col=1:4,pch=15:18,text.col=1:4,bty='n',cex=1)
	}
	dev.off()
}

#species biplot
sp.biplot <- function(stmp){
	jpeg(paste("sp",stmp,'biplot.jpg',sep="."),width=6,height=5,units='in',res=1000)
  par(mfrow=c(2,2),mar=c(5,5,1,1))
   l <- 1 

    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(gvec[tmp],agvec[tmp],xlab=expression(italic(G[ref])),#main=spname[which(spall==stmp)],
	      ylab=expression(lambda/italic(G[ref])),pch=20,xlim=range(gvec[splist==stmp]),col='white',
		  ylim=range(agvec[splist==stmp]),log='xy',cex.lab=1.4)
      
	    #points(gvec[tmp],agvec[tmp],pch=14+yy,col=yy,cex=1.5)
	    text(gvec[tmp],agvec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
	    #arrows(gvec[tmp[-length(tmp)]],agvec[tmp[-length(tmp)]],
	    #       gvec[tmp[-1]],agvec[tmp[-1]],length=.075,col=yy)
	  l<-0
	}

    l <- 1
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(l1vec[tmp],l2vec[tmp],xlab=expression(alpha[1]),#main=spname[which(spall==stmp)],
	      ylab=expression(alpha[2]),pch=20,xlim=range(l1vec[splist==stmp]),
		  ylim=range(l2vec[splist==stmp]),log='xy',cex.lab=1.4,col='white')
      
	    #points(l1vec[tmp],l2vec[tmp],pch=14+yy,col=yy,cex=1.5)
	    text(l1vec[tmp],l2vec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
	    #arrows(l1vec[tmp[-length(tmp)]],l2vec[tmp[-length(tmp)]],
	    #       l1vec[tmp[-1]],l2vec[tmp[-1]],length=.075,col=yy)
	  l<-0
    }
	
    l <- 1
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(m1vec[tmp],m2vec[tmp],xlab=expression(alpha[4]),#main=spname[which(spall==stmp)],
	      ylab=expression(alpha[3]),pch=20,xlim=range(m1vec[splist==stmp]),col='white',
		  ylim=range(m2vec[splist==stmp]),log='xy',cex.lab=1.4)
      
	    #points(m1vec[tmp],m2vec[tmp],pch=14+yy,col=yy,cex=1.5)
	    text(m1vec[tmp],m2vec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
	    #arrows(m1vec[tmp[-length(tmp)]],m2vec[tmp[-length(tmp)]],
	    #       m1vec[tmp[-1]],m2vec[tmp[-1]],length=.075,col=yy)
	  l<-0
    }
	
    l <- 1
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
	  if(length(tmp)==0) next
      if(l==1)
        plot(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),#main=spname[which(spall==stmp)],
		  xlab=expression(sqrt(italic(V[J]))),col='white',
	      ylab=expression(sigma),pch=20,xlim=range(sqrt(vjvec[splist==stmp])),
		  ylim=range(sqrt(sigvec[splist==stmp])),log='xy',cex.lab=1.4)
      
	    #points(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),pch=14+yy,col=yy,cex=1.5)
	    text(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),labels=round(kplist[tmp]*30,0),col=yy)
	    #arrows(sqrt(vjvec[tmp[-length(tmp)]]),sqrt(sigvec[tmp[-length(tmp)]]),
	    #       sqrt(vjvec[tmp[-1]]),sqrt(sigvec[tmp[-1]]),length=.075,col=yy)
	  l<-0
	}
	legend('topright',legend=yrall,text.col=1:4,bty='n')

	dev.off()
}

#example
sp.biplot <- function(stmp){
	jpeg(paste("sp",stmp,'biplot.ex.jpg',sep="."),width=4,height=3.5,units='in',res=1000)
  par(mfrow=c(2,2),mar=c(4,5,1,1))
   l <- 1 

    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
		tmp <- tmp[kplist[tmp] %in% kpall]
	  if(length(tmp)==0) next
      if(l==1)
        plot(gvec[tmp],agvec[tmp],xlab=expression(italic(G[ref])),#main=spname[which(spall==stmp)],
	      ylab=expression(lambda/italic(G[ref])),pch=20,xlim=range(gvec[splist==stmp]),col='white',
		  ylim=range(agvec[splist==stmp]),log='xy',cex.lab=1.4)
      
	    points(gvec[tmp],agvec[tmp],pch=14+which(spall==stmp),col=yy,cex=.8)
	    #text(gvec[tmp],agvec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
	    arrows(gvec[tmp[-length(tmp)]],agvec[tmp[-length(tmp)]],
	           gvec[tmp[-1]],agvec[tmp[-1]],length=.075,col=yy)
	  l<-0
	}

    l <- 1
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
		tmp <- tmp[kplist[tmp] %in% kpall]
	  if(length(tmp)==0) next
      if(l==1)
        plot(l1vec[tmp],l2vec[tmp],xlab=expression(alpha[1]),#main=spname[which(spall==stmp)],
	      ylab=expression(alpha[2]),pch=20,xlim=range(l1vec[splist==stmp]),
		  ylim=range(l2vec[splist==stmp]),log='xy',cex.lab=1.4,col='white')
      
	    points(l1vec[tmp],l2vec[tmp],pch=14+which(spall==stmp),col=yy,cex=.8)
	    #text(l1vec[tmp],l2vec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
	    arrows(l1vec[tmp[-length(tmp)]],l2vec[tmp[-length(tmp)]],
	           l1vec[tmp[-1]],l2vec[tmp[-1]],length=.075,col=yy)
	  l<-0
    }
	
    l <- 1
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
		tmp <- tmp[kplist[tmp] %in% kpall]
	  if(length(tmp)==0) next
      if(l==1)
        plot(m2vec[tmp],m1vec[tmp],xlab=expression(alpha[3]),#main=spname[which(spall==stmp)],
	      ylab=expression(alpha[4]),pch=20,xlim=range(m2vec[splist==stmp]),col='white',
		  ylim=range(m1vec[splist==stmp]),log='xy',cex.lab=1.4)
      
	    points(m2vec[tmp],m1vec[tmp],pch=14+which(spall==stmp),col=yy,cex=.8)
	    #text(m1vec[tmp],m2vec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
	    arrows(m2vec[tmp[-length(tmp)]],m1vec[tmp[-length(tmp)]],
	           m2vec[tmp[-1]],m1vec[tmp[-1]],length=.075,col=yy)
	  l<-0
    }
	
    l <- 1
    for(yy in 1:4){
      ytmp <- sort(unique(yrlist))[yy]
      tmp  <- which(splist==stmp & yrlist==ytmp) 
        tmp <- tmp[order(kplist[tmp])]
		tmp <- tmp[kplist[tmp] %in% kpall]
	  if(length(tmp)==0) next
      if(l==1)
        plot(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),#main=spname[which(spall==stmp)],
		  xlab=expression(sqrt(italic(V[J]))),col='white',
	      ylab=expression(sigma),pch=20,xlim=range(sqrt(vjvec[splist==stmp])),
		  ylim=range(sqrt(sigvec[splist==stmp])),log='xy',cex.lab=1.4)
      
	    points(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),pch=14+which(spall==stmp),col=yy,cex=.8)
	    #text(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),labels=round(kplist[tmp]*30,0),col=yy)
	    arrows(sqrt(vjvec[tmp[-length(tmp)]]),sqrt(sigvec[tmp[-length(tmp)]]),
	           sqrt(vjvec[tmp[-1]]),sqrt(sigvec[tmp[-1]]),length=.075,col=yy)
	  l<-0
	}
	legend('topright',legend=yrall,text.col=1:4,bty='n',cex=.7)

	dev.off()
}


#example
all.kappa.effects <- function(stmp){
  jpeg(paste("sp",stmp,'biplot.ex.jpg',sep="."),width=4,height=3.5,units='in',res=1000)
  par(mfrow=c(2,2),mar=c(4,5,1,1))
  l <- 1 
  
  for(yy in 1:4){
    ytmp <- sort(unique(yrlist))[yy]
    tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
    tmp <- tmp[kplist[tmp] %in% kpall[c(1,3,5,7,9)]]
    if(length(tmp)==0) next
    if(l==1)
      plot(gvec[tmp],agvec[tmp],xlab=expression(italic(G[ref])),#main=spname[which(spall==stmp)],
           ylab=expression(lambda/italic(G[ref])),pch=20,xlim=range(gvec[splist==stmp]),col='white',
           ylim=range(agvec[splist==stmp]),log='xy',cex.lab=1.4)
    
    points(gvec[tmp],agvec[tmp],pch=14+which(spall==stmp),col=yy,cex=.8)
    #text(gvec[tmp],agvec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
    arrows(gvec[tmp[-length(tmp)]],agvec[tmp[-length(tmp)]],
           gvec[tmp[-1]],agvec[tmp[-1]],length=.075,col=yy)
    l<-0
  }
  
  l <- 1
  for(yy in 1:4){
    ytmp <- sort(unique(yrlist))[yy]
    tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
    tmp <- tmp[kplist[tmp] %in% kpall[c(1,3,5,7,9)]]
    if(length(tmp)==0) next
    if(l==1)
      plot(l1vec[tmp],l2vec[tmp],xlab=expression(alpha[1]),#main=spname[which(spall==stmp)],
           ylab=expression(alpha[2]),pch=20,xlim=range(l1vec[splist==stmp]),
           ylim=range(l2vec[splist==stmp]),log='xy',cex.lab=1.4,col='white')
    
    points(l1vec[tmp],l2vec[tmp],pch=14+which(spall==stmp),col=yy,cex=.8)
    #text(l1vec[tmp],l2vec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
    arrows(l1vec[tmp[-length(tmp)]],l2vec[tmp[-length(tmp)]],
           l1vec[tmp[-1]],l2vec[tmp[-1]],length=.075,col=yy)
    l<-0
  }
  
  l <- 1
  for(yy in 1:4){
    ytmp <- sort(unique(yrlist))[yy]
    tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
    tmp <- tmp[kplist[tmp] %in% kpall[c(1,3,5,7,9)]]
    if(length(tmp)==0) next
    if(l==1)
      plot(m2vec[tmp],m1vec[tmp],xlab=expression(alpha[3]),#main=spname[which(spall==stmp)],
           ylab=expression(alpha[4]),pch=20,xlim=range(m2vec[splist==stmp]),col='white',
           ylim=range(m1vec[splist==stmp]),log='xy',cex.lab=1.4)
    
    points(m2vec[tmp],m1vec[tmp],pch=14+which(spall==stmp),col=yy,cex=.8)
    #text(m1vec[tmp],m2vec[tmp],labels=round(kplist[tmp]*30,0),col=yy)
    arrows(m2vec[tmp[-length(tmp)]],m1vec[tmp[-length(tmp)]],
           m2vec[tmp[-1]],m1vec[tmp[-1]],length=.075,col=yy)
    l<-0
  }
  
  l <- 1
  for(yy in 1:4){
    ytmp <- sort(unique(yrlist))[yy]
    tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
    tmp <- tmp[kplist[tmp] %in% kpall[c(1,3,5,7,9)]]
    if(length(tmp)==0) next
    if(l==1)
      plot(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),#main=spname[which(spall==stmp)],
           xlab=expression(sqrt(italic(V[J]))),col='white',
           ylab=expression(sigma),pch=20,xlim=range(sqrt(vjvec[splist==stmp])),
           ylim=range(sqrt(sigvec[splist==stmp])),log='xy',cex.lab=1.4)
    
    points(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),pch=14+which(spall==stmp),col=yy,cex=.8)
    #text(sqrt(vjvec[tmp]),sqrt(sigvec[tmp]),labels=round(kplist[tmp]*30,0),col=yy)
    arrows(sqrt(vjvec[tmp[-length(tmp)]]),sqrt(sigvec[tmp[-length(tmp)]]),
           sqrt(vjvec[tmp[-1]]),sqrt(sigvec[tmp[-1]]),length=.075,col=yy)
    l<-0
  }
  legend('topright',legend=yrall,text.col=1:4,bty='n')
  
  dev.off()
}

#
par.by.kappa <- function(){
pdf('par.by.kappa.pdf',width=9,height=9)	
for(stmp in c('cato','list','litu','qual','qumi','quph')){

par(mfrow=c(3,3),mar=c(4,4,1,1))
l <- 1

move <- seq(-1,1,length=length(yrall))

for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],gvec[tmp],type='l',col=j,ylim=range(gvec),
             xlab='time constant',ylab=expression(italic(G[ref])))
    points(kplist[tmp]*30+move[j],gvec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],gvec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],glo[tmp],kplist[tmp]*30+move[j],ghi[tmp],col=j)
	l <- 0
  }
l <- 1
  
for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],agvec[tmp],type='l',col=j,ylim=c(.4,.85),
             xlab='time constant',ylab=expression(lambda/italic(G[ref])))
    points(kplist[tmp]*30+move[j],agvec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],agvec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],aglo[tmp],kplist[tmp]*30+move[j],aghi[tmp],col=j)
	l <- 0  
	}
l <- 1
  
for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],l1vec[tmp],type='l',col=j,ylim=range(l1vec),
             xlab='time constant',ylab=expression(italic(L[1])))
    points(kplist[tmp]*30+move[j],l1vec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],l1vec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],l1lo[tmp],kplist[tmp]*30+move[j],l1hi[tmp],col=j)
	l <- 0  
	}
l <- 1
  
for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],l2vec[tmp],type='l',col=j,ylim=range(l2vec),
             xlab='time constant',ylab=expression(italic(L[2])))
    points(kplist[tmp]*30+move[j],l2vec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],l2vec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],l2lo[tmp],kplist[tmp]*30+move[j],l2hi[tmp],col=j)
	l <- 0  
	}
l <- 1

for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],m1vec[tmp],type='l',col=j,ylim=range(m1vec),
             xlab='time constant',ylab=expression(italic(M[1])))
    points(kplist[tmp]*30+move[j],m1vec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],m1vec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],m1lo[tmp],kplist[tmp]*30+move[j],m1hi[tmp],col=j)
	l <- 0  
	}
l <- 1

for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],m2vec[tmp],type='l',col=j,ylim=c(min(unlist(M.all)),max(unlist(M.all))),
             xlab='time constant',ylab=expression(italic(M[2])))
    points(kplist[tmp]*30+move[j],m2vec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],m2vec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],m2lo[tmp],kplist[tmp]*30+move[j],m2hi[tmp],col=j)
	l <- 0  
	}
legend('topright',legend=yrall,text.col=1:4)
l <- 1

for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],sigvec[tmp],type='l',col=j,ylim=range(sigvec),
             xlab='time constant',ylab=expression(sigma^2))
    points(kplist[tmp]*30+move[j],sigvec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],sigvec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],siglo[tmp],kplist[tmp]*30+move[j],sighi[tmp],col=j)
	l <- 0  
	}
l <- 1

for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],vjvec[tmp],type='l',col=j,ylim=range(vjvec),
             xlab='time constant',ylab=expression(italic(V[j])^2))
    points(kplist[tmp]*30+move[j],vjvec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],vjvec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],vjlo[tmp],kplist[tmp]*30+move[j],vjhi[tmp],col=j)
	l <- 0  
	}
l <- 1

for(j in 1:4) {
  ytmp <- unique(yrlist)[j]
  tmp  <- which(splist==stmp & yrlist==ytmp) 
    tmp <- tmp[order(kplist[tmp])]
  if(length(tmp)==0) next
  if(l==1) plot(kplist[tmp]*30+move[j],vavec[tmp],type='l',col=j,ylim=c(min(valo),max(vahi)),
             xlab='time constant',ylab=expression(italic(V[a])^2))
    points(kplist[tmp]*30+move[j],vavec[tmp],col=j)
    lines(kplist[tmp]*30+move[j],vavec[tmp],type='l',col=j)
	segments(kplist[tmp]*30+move[j],valo[tmp],kplist[tmp]*30+move[j],vahi[tmp],col=j)
	l <- 0  
	}
    }  
  dev.off()
}

	  
##biplot Gref and lambda
gref.lam.rat.jpg <- function(){
jpeg('gref.lambda.ratio.1.1.line.jpg',width=6.5,height=2,units='in',res=1000)
par(mar=c(5,5,1,1),mfrow=c(1,3))

gtmp <- atmp <- agtmp <- matrix(NA,dim(spyr)[1],3)

for(j in 1:dim(spyr)[1]){
  stmp <- spyr[j,1]
  ytmp <- spyr[j,2]

      tmp  <- which(splist==stmp & yrlist==ytmp) 
      tmp <- tmp[which(kplist[tmp]==tc)]
    
	if(length(tmp)==0) next
    
	atmp[j,]  <- quantile(agibbs.all[[tmp]][burnin:ng,],c(.5,.025,.975))
	gtmp[j,]  <- quantile(ggibbs.all[[tmp]][burnin:ng,],c(.5,.025,.975))
	agtmp[j,] <- quantile(agibbs.all[[tmp]][burnin:ng,]/ggibbs.all[[tmp]][burnin:ng,],c(.5,.025,.975))
  }


pbla <- c(2,3,4)
ptmp <- c('gtmp','atmp','agtmp')
lo   <- c(5,5,.2)
hi   <- c(240,150,1)
pname <- c(expression(italic(G[ref])),expression(lambda),expression(lambda/italic(G[ref])))
for(i in 1:3){
  mt <- 1
  if(i==2) mt <- -1
  for(j in 1:4) {
	  tmpsp <- spyr[as.numeric(spyr[,2])==yrall[j],1]
	  blsp  <- names(tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean))
	  bmatch <- which(blsp %in% tmpsp)
	if(j==1)
	  plot(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
	       get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],
		   xlim=c(lo[i],hi[i]),ylim=c(lo[i],hi[i]),pch=PT[match(blsp[bmatch],spall)],log='xy',
		   xlab='boundary line estimates',ylab='state-space estimates',col=j,main=pname[i],cex=1.2,cex.lab=1.4)
      points(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
	       get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],pch=PT[match(blsp[bmatch],spall)],col=j,cex=1.2)

      segments(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],quantile,.025,na.rm=T)[bmatch],
	          get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],
	          mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],quantile,.975,na.rm=T)[bmatch],
			  get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],col=j)
      segments(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
	          get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],2],
	          mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
			  get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],3],col=j)
	  }
	  abline(0,1,lty=2)

   if(i==1) legend("topleft", legend=yrall, text.col=1:4,bty='n',cex=.75)
   if(i==3) legend("bottomright", legend=spname, pch=PT, bty='n',cex=.75)
 
  }
	  
#	splm <- lm(atmp[,1]~gtmp[,1])$coef
#	segments(min(gtmp[,1]),splm[1]+splm[2]*min(gtmp[,1]),max(gtmp[,1]),splm[1]+splm[2]*max(gtmp[,1]))
  
 
  dev.off()
## 
  list(atmp = atmp, gtmp = gtmp, agtmp = agtmp)
}


##biplot Gref and lambda
gref.lam.rat.jpg2 <- function(){
jpeg('gref.lambda.ratio2.1.1.line.jpg',width=6.5,height=3.25,units='in',res=1000)
par(mar=c(5,5,1,1),mfrow=c(1,2))

gtmp <- atmp <- agtmp <- matrix(NA,dim(spyr)[1],3)

for(j in 1:dim(spyr)[1]){
  stmp <- spyr[j,1]
  ytmp <- spyr[j,2]

      tmp  <- which(splist==stmp & yrlist==ytmp) 
      tmp <- tmp[which(kplist[tmp]==tc)]
    
	if(length(tmp)==0) next
    
	atmp[j,]  <- quantile(agibbs.all[[tmp]][burnin:ng,],c(.5,.025,.975))
	gtmp[j,]  <- quantile(ggibbs.all[[tmp]][burnin:ng,],c(.5,.025,.975))
	agtmp[j,] <- quantile(agibbs.all[[tmp]][burnin:ng,]/ggibbs.all[[tmp]][burnin:ng,],c(.5,.025,.975))
  }


pbla <- c(2,3)
ptmp <- c('gtmp','atmp')
lo   <- c(10,5)
hi   <- c(250,180)
pname <- c(expression(italic(G[ref])),expression(lambda),expression(lambda/italic(G[ref])))
for(i in 1:2){
  mt <- 1
  if(i==2) mt <- -1
  for(j in 1:4) {
	  tmpsp <- spyr[as.numeric(spyr[,2])==yrall[j],1]
	  blsp  <- names(tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean))
	  bmatch <- which(blsp %in% tmpsp)
	if(j==1)
	  plot(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
	       get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],
		   xlim=c(lo[i],hi[i]),ylim=c(lo[i],hi[i]),pch=PT[match(blsp[bmatch],spall)],log='xy',
		   xlab='boundary line estimates',ylab='state-space estimates',col=j,main=pname[i],cex=1.2,cex.lab=1.4)
      points(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
	       get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],pch=PT[match(blsp[bmatch],spall)],col=j,cex=1.2)

      nseg  <- tapply(BLA[[j]][,pbla[i]][is.finite(BLA[[j]][,pbla[i]])],
	             spBLA[[j]][is.finite(BLA[[j]][,pbla[i]])],length)[bmatch] 
	  loseg <- mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch] - 
	            1.96*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],sd,na.rm=T)[bmatch]/sqrt(nseg)
	  hiseg <- mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch] + 
	            1.96*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],sd,na.rm=T)[bmatch]/sqrt(nseg)
				
	    loseg[loseg<1] <- 1
		
	  
	  segments(loseg,
	          get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],
	          hiseg,
			  get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],1],col=j)
      segments(mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
	          get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],2],
	          mt*tapply(BLA[[j]][,pbla[i]],spBLA[[j]],mean,na.rm=T)[bmatch],
			  get(ptmp[i])[as.numeric(spyr[,2])==yrall[j],3],col=j)
	  }
	  abline(0,1,lty=2)

   if(i==2) legend(lo[i],hi[i]*1.1, legend=yrall, text.col=1:4, title='Year',bty='n',cex=.8)
   if(i==1) legend(lo[i],hi[i]*1.1, legend=spname, pch=PT, bty='n',cex=.8)
 
  }
	  
#	splm <- lm(atmp[,1]~gtmp[,1])$coef
#	segments(min(gtmp[,1]),splm[1]+splm[2]*min(gtmp[,1]),max(gtmp[,1]),splm[1]+splm[2]*max(gtmp[,1]))
  
 
  dev.off()
## 
  list(atmp = atmp, gtmp = gtmp, agtmp = agtmp)
}

##get tree sizes
get.dbh <- function(){
   wd <-  strsplit(getwd(),'/')
    wd <- unlist(wd)[length(unlist(wd))]

  setwd('../')
  setwd('datafiles')

  f   <- list.files()
  dm  <- f[grep('id.csv',f)]
  
  b2 <- read.csv(dm,header=T)
  dbh <- b2[,c('colname',colnames(b2)[grep('DBH',colnames(b2))])]
    
   
  setwd('../')
  setwd(wd)
  
  dbh
  }
##

gref.lambda <- function(){ 
jpeg('gref.lambda.jpg',width=3.5,height=3.5,units='in',res=1000)
par(mar=c(5,5,1,1),mfrow=c(1,1))
  plot(gtmp[,1],atmp[,1],xlab=xl,ylab=yl,
    pch=PT[match(spyr[,1],spall)],col=match(spyr[,2],yrall),log='xy')
  segments(gtmp[,1],atmp[,2],gtmp[,1],atmp[,3],col=match(spyr[,2],yrall))
  segments(gtmp[,2],atmp[,1],gtmp[,3],atmp[,1],col=match(spyr[,2],yrall))
  x <- seq(1,140,by=1)
  lines(x,x*.6)

  legend(80,30, legend=yrall, text.col=1:4,bty='n',cex=.8)
  legend(15,80, legend=spname, title='Species',pch=PT, bty='n',cex=.8)
dev.off()
}
##

##response functions

resp <- function(){
samp <- sample(burnin:ng,2000,replace=T)

Qvec <- array(NA,dim=c(dim(spyr)[1],dim(Qm)[1],length(samp))) 
Mvec <- Dvec <- Qvec
tc <- 1

for(j in 1:length(samp)) {
  for(stmp in spall){
    for(ytmp in yrall) {
      tmp  <- which(splist==stmp & yrlist==ytmp & kplist==tc)
	  keep <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	  if(length(keep)==0) next
	  
	  blitespec <- matrix(lgibbs.all[[tmp]][samp[j],],2,1)
	  bmoistspec <- matrix(mgibbs.all[[tmp]][samp[j],],2,1)
	  gspec <- ggibbs.all[[tmp]][samp[j],]
	  aspec <- agibbs.all[[tmp]][samp[j],]
	  sigma <- vgibbs.all[[tmp]][samp[j],1]
	  
	  Dvec[keep,,j] <- gspec - aspec*log(Dm[,keep])
	  
	    Mtmp <- exp(-0.5 * ((Mm[,keep] - bmoistspec[2])/bmoistspec[1])^2)
	    Mtmp[Mm[,keep] > bmoistspec[2]] <- 1
      Mvec[keep,,j] <- Mtmp # rnorm(dim(Mm)[1],gspec*Mtmp,sqrt(sigma))
	  
	  Qvec[keep,,j] <- 1 - blitespec[1]*exp(-Qm[,keep]/blitespec[2])
	  
	  }
	}
  }
list(Qvec = Qvec, Mvec = Mvec, Dvec = Dvec)
  }
  
##plotting relations
resp.by.year.plot <- function(){
jpeg('Dresp.by.year.jpg',width=5,height=5,unit='in',res=1000)
par(mfrow=c(2,2),mar=c(5,5,2,1))
for(ytmp in yrall){
  l <- 1
  for(stmp in spall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	if(length(tmp)==0) next
	
	ciD <- apply(Dvec[tmp,,],1,quantile,c(.5,.025,.975))
	
	if(l==1) plot(Dm[,tmp],ciD[1,],col=which(spall==stmp),type='l',ylim=c(min(Dvec),max(Dvec)),
	              xlab=expression(italic(D[t])),ylab='Relative Response',
				  main=ytmp,xlim=c(.6,max(Dm)))
	  lines(Dm[,tmp],ciD[1,],col=which(spall==stmp))
	  lines(Dm[,tmp],ciD[2,],col=which(spall==stmp),lty=2)
	  lines(Dm[,tmp],ciD[3,],col=which(spall==stmp),lty=2)
    l <- 0
	}
  }	  
  legend('topright',legend=spall,text.col=1:6,bty='n')
  dev.off()

jpeg('Qresp.by.year.jpg',width=5,height=5,unit='in',res=1000)
par(mfrow=c(2,2),mar=c(5,5,2,1))
for(ytmp in yrall){
  l <- 1
  for(stmp in spall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	if(length(tmp)==0) next
	
	ciQ <- apply(Qvec[tmp,,],1,quantile,c(.5,.025,.975))
	
	if(l==1) plot(Qm[,tmp],ciQ[1,],col=which(spall==stmp),type='l',ylim=c(min(Qvec),max(Qvec)),
	              xlab=expression(italic(Q[t])),ylab='Relative Response',log='xy',
				  main=ytmp)
	  lines(Qm[,tmp],ciQ[1,],col=which(spall==stmp))
	  lines(Qm[,tmp],ciQ[2,],col=which(spall==stmp),lty=2)
	  lines(Qm[,tmp],ciQ[3,],col=which(spall==stmp),lty=2)
    l <- 0
	}
  }	  
  legend('bottomright',legend=spall,text.col=1:6,bty='n')
  dev.off()
  
  
##plotting relations
jpeg('Mresp.by.year.jpg',width=5,height=5,unit='in',res=1000)
par(mfrow=c(2,2),mar=c(5,5,2,1))
for(ytmp in yrall){
  l <- 1
  for(stmp in spall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	if(length(tmp)==0) next
	
	ciM <- apply(Mvec[tmp,,],1,quantile,c(.5,.025,.975))
	
	if(l==1) plot(Mm[,tmp],ciM[1,],col=which(spall==stmp),type='l',ylim=c(min(Mvec),max(Mvec)),
	              xlab=expression(italic(M[t])),ylab='Relative Response',
				  main=ytmp,xlim=c(0.1,.4))
	  lines(Mm[,tmp],ciM[1,],col=which(spall==stmp))
	  lines(Mm[,tmp],ciM[2,],col=which(spall==stmp),lty=2)
	  lines(Mm[,tmp],ciM[3,],col=which(spall==stmp),lty=2)
    l <- 0
	}
  }	  
  legend('bottomright',legend=spall,text.col=1:6,bty='n')
  dev.off()
  }
  

##plotting relations
resp.by.species.plot <- function(){
jpeg('Dresp.by.sp.jpg',width=6.5,height=4,unit='in',res=1000)
par(mfrow=c(2,3),mar=c(5,5,2,1))
  for(stmp in spall){
  l <- 1
    for(ytmp in yrall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	if(length(tmp)==0) next
	
	ciD <- apply(Dvec[tmp,,],1,quantile,c(.5,.025,.975))
	
	if(l==1) plot(Dm[,tmp],ciD[1,],col=which(yrall==ytmp),type='l',ylim=c(min(Dvec),max(Dvec)),
	              main=spname[which(spall==stmp)],ylab='Relative Response',log='x',
				  xlab=expression(Vapor~~Pressure~~Deficit~~italic(D[t])~~(kPa)),xlim=c(0.6,max(Dm)))
	  lines(Dm[,tmp],ciD[1,],col=which(yrall==ytmp))
	  lines(Dm[,tmp],ciD[2,],col=which(yrall==ytmp),lty=2)
	  lines(Dm[,tmp],ciD[3,],col=which(yrall==ytmp),lty=2)
    l <- 0
	}
  }	  
  legend('topright',legend=yrall,text.col=1:4,bty='n')
  dev.off()
  
jpeg('Qresp.by.sp.jpg',width=6.5,height=4,unit='in',res=1000)
par(mfrow=c(2,3),mar=c(5,5,2,1))
  for(stmp in spall){
  l <- 1
    for(ytmp in yrall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	if(length(tmp)==0) next
	
	ciQ <- apply(Qvec[tmp,,],1,quantile,c(.5,.025,.975))
	
	if(l==1) plot(Qm[,tmp],ciQ[1,],col=which(yrall==ytmp),type='l',ylim=c(min(Qvec),max(Qvec)),
	              xlab="",ylab='Relative Response',log='x',
				  main=spname[spall==stmp],xaxt='n')
	  axis(1,at=c(2,20,200,2000))
	  lines(Qm[,tmp],ciQ[1,],col=which(yrall==ytmp))
	  lines(Qm[,tmp],ciQ[2,],col=which(yrall==ytmp),lty=2)
	  lines(Qm[,tmp],ciQ[3,],col=which(yrall==ytmp),lty=2)
    l <- 0
	}
  }	  
  legend('bottomright',legend=yrall,text.col=1:4,bty='n')
  dev.off()
  
  
##plotting relations
jpeg('Mresp.by.sp.jpg',width=6.5,height=4,unit='in',res=1000)
par(mfrow=c(2,3),mar=c(5,5,2,1))
  for(stmp in spall){
  l <- 1
    for(ytmp in yrall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
	if(length(tmp)==0) next
	
	ciM <- apply(Mvec[tmp,,],1,quantile,c(.5,.025,.975))
	
	if(l==1) plot(Mm[,tmp],ciM[1,],col=which(yrall==ytmp),type='l',ylim=c(min(Mvec),max(Mvec)),
	              xlab=expression(italic(M[t])),ylab='Relative Response',
				  main=spname[which(spall==stmp)],xlim=c(0.1,.35))
	  lines(Mm[,tmp],ciM[1,],col=which(yrall==ytmp))
	  lines(Mm[,tmp],ciM[2,],col=which(yrall==ytmp),lty=2)
	  lines(Mm[,tmp],ciM[3,],col=which(yrall==ytmp),lty=2)
    l <- 0
	}
  }	  
  legend('bottomright',legend=yrall,text.col=1:4,bty='n')
  dev.off()
  }
    





##plotting relations combined
resp.by.species.plot <- function(w,h,cc=1){

  let <- c('(a)','(b)','(c)','(d)','(e)','(f)',   #D response
           '(g)','(h)','(i)','(j)','(k)','(l)',   #Q response
           '(m)','(n)','(o)','(p)','(q)','(r)')  #M response
  
  postscript('combined.species.responses.eps',width=w,height=h,onefile=FALSE,horizontal=FALSE)

  plot.new()
  
  
    for(stmp in spall){
      
      par(new = "TRUE",plt = plt.mat[(which(spall == stmp)-1)*3 +1,],las = 1,cex.axis = 1)

      l <- 1
      for(ytmp in yrall){
        tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
        if(length(tmp)==0) next
        
        ciD <- apply(Dvec[tmp,,],1,quantile,c(.5,.025,.975))
        
        if(l==1) plot(Dm[,tmp],ciD[1,],col=which(yrall==ytmp),type='l',ylim=c(min(Dvec),max(Dvec)),
                      main="",ylab='',axes=FALSE,frame.plot=TRUE,
                      xlab="",xlim=c(0.6,max(Dm)))
        lines(Dm[,tmp],ciD[1,],col=which(yrall==ytmp))
        lines(Dm[,tmp],ciD[2,],col=which(yrall==ytmp),lty=2)
        lines(Dm[,tmp],ciD[3,],col=which(yrall==ytmp),lty=2)
        
        l <- 0  
      }

      legend('topright',legend=let[which(spall == stmp)],bty='n',cex=cc)
      
      axis(2,at=c(0,50,100),labels=rep("",3),tck=.05,cex=cc)
      mtext(c(0,50,100),side=2,at=c(0,50,100),cex=cc)
      
      if(spall[1] == stmp){
        legend('topleft',legend=yrall,text.col=1:4,bty='n',cex=cc)
        mtext("Vapor Pres.",side=3,line=.75,cex=cc)
        mtext("Deficit",side=3,cex=cc)
        
      }
      
      if(spall[6] == stmp){
        axis(1,at=1:4,labels=rep("",4),tck=.05,cex=cc)
        mtext(1:4,side=1,at=1:4,cex=cc)
        }
    }
  mtext(expression(italic(D[t])~~(kPa)),side=1,line=1,cex=cc)
  
  mtext(expression(Steady-State ~~ Canopy~~Conductance~~italic(G[list(S,t)])~~(mmol~~ m^-2 ~~ s^-1)),
        side=2,outer=TRUE,line=-1.1,las=3,cex=cc)

  for(stmp in spall){
    
  par(new = "TRUE",plt = plt.mat[(which(spall == stmp)-1)*3 +2,],las = 1,cex.axis = 1)
  
  l <- 1
  for(ytmp in yrall){
    tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
    if(length(tmp)==0) next
    
    ciQ <- apply(Qvec[tmp,,],1,quantile,c(.5,.025,.975))
    
    if(l==1) plot(Qm[,tmp],ciQ[1,],col=which(yrall==ytmp),type='l',ylim=c(min(Qvec),max(Qvec)),
                  xlab="",ylab='',axes=FALSE,frame.plot=TRUE,
                  main="")
    lines(Qm[,tmp],ciQ[1,],col=which(yrall==ytmp))
    lines(Qm[,tmp],ciQ[2,],col=which(yrall==ytmp),lty=2)
    lines(Qm[,tmp],ciQ[3,],col=which(yrall==ytmp),lty=2)
    l <- 0
  }

  axis(2,at=c(0.0,0.25,.5,.75,1),labels=rep("",5),tck=.05,cex=cc)
  mtext(c(0.25,0.75),side=2,at=c(0.25,0.75),cex=cc)
 
  if(spall[1] == stmp)   mtext("Light",side=3,line=.75,cex=cc)

  legend('topright',legend=let[6 + which(spall == stmp)],bty='n',cex=cc)
  
  if(spall[6] == stmp){
    axis(1,at=seq(0,2000,by=500),labels=rep("",5),tck=.05,cex=cc)
    mtext(c(500,1500),side=1,at=c(500,1500),cex=cc)
  }
  
  
  } 

  mtext(expression(italic(Q[t])),side=1,line=1,cex=cc)
  mtext(expression((mu~mol~~m^-2~~sec^-1)),side=1,line=2,cex=cc)
  
  mtext("Relative Response",side=2,outer=TRUE,line=-6.8,las=3,cex=cc)
  
  for(stmp in spall){
    par(new = "TRUE",plt = plt.mat[(which(spall == stmp)-1)*3 +3,],las = 1,cex.axis = 1)

    l <- 1
    for(ytmp in yrall){
      tmp <- which(spyr[,1]==stmp & spyr[,2]==ytmp)
      if(length(tmp)==0) next
      
      ciM <- apply(Mvec[tmp,,],1,quantile,c(.5,.025,.975))
      
      if(l==1) plot(Mm[,tmp],ciM[1,],col=which(yrall==ytmp),type='l',ylim=c(min(Mvec),max(Mvec)),
                    xlab="",ylab="",axes=FALSE,frame.plot=TRUE,
                    main="",xlim=c(0.1,.35))
      lines(Mm[,tmp],ciM[1,],col=which(yrall==ytmp))
      lines(Mm[,tmp],ciM[2,],col=which(yrall==ytmp),lty=2)
      lines(Mm[,tmp],ciM[3,],col=which(yrall==ytmp),lty=2)
      l <- 0
    }
    
    legend('topright',legend=let[12 + which(spall == stmp)],bty='n',cex=cc)
    
    if(spall[1] == stmp)   mtext("Moisture",side=3,line=.75,cex=cc)

    if(spall[6] == stmp){
      axis(1,at=c(.1,.2,.3),labels=rep("",3),tck=.05,cex=cc)
      mtext(c(.1,.2,.3),side=1,at=c(.1,.2,.3),cex=cc)
    }
    
    mtext(spname[which(spall == stmp)],side=4,las=3,cex=cc)
    
  }    
  
  mtext(expression(italic(M[t])~~(m~~m^-1)),side=1,line=1,cex=cc)
  
   dev.off()
}


Dpar.by.time.bar <- function(){
  jpeg('Dpar.by.time.bar.jpg',width=3.5,height=4,units='in',res=1000) 
par(mfrow=c(3,1),mar=c(3,4,1,0))
#  plot(match(spyr[,1],spall)+spoff[match(spyr[,2],yrall)],gtmp[,1],pch=PT[match(spyr[,1],spall)],xaxt='n',
#    ylab=expression(italic(G[ref])),xlab="")#,xlim=c(min(yrall)-.5,max(yrall)+.5))

count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- atmp[,1]
  bx <- as.vector(barplot(count,beside=T,col=1:4,
    ylab=expression(italic(G[ref])),plot=FALSE))
    
    bx.seq <- c(1,6,11,17,2,7,12,18,3,8,13,19,4,9,14,20,15,21,5,10,16,22)
    
  colbx <- rep(c('black','red','green','blue'),times=6)
  colbx <- colbx[-c(17,18)]
  bx <- bx[-c(17,18)]
  
plot(bx,gtmp[bx.seq,1],pch=20,frame.plot=TRUE,axes=FALSE,ylab="",xlab="",col=colbx,log='y',
		ylim=c(10,120))

  arrows(bx,gtmp[bx.seq,2],bx,gtmp[bx.seq,3],angle=90,length=.05,lwd=2,col=colbx)
  
  axis(2,c(10,20,40,80,120))
  
count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- atmp[,1]
  bx <- barplot(t(count),beside=T,col=1:4,
    ylab=expression(lambda))
  
counthi <- matrix(NA,nsp,4)
  counthi[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- atmp[,3]
  arrows(bx,t(count),bx,t(counthi),angle=90,length=.05,lwd=2)

  
count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- agtmp[,1]
  bx <- barplot(t(count),beside=T,col=1:4,
    ylab=expression(lambda/italic(G[ref])),ylim=c(.4,.8),xpd=F)
  
counthi <- matrix(NA,nsp,4)
  counthi[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- agtmp[,3]
  arrows(bx,t(count),bx,t(counthi),angle=90,length=.05,lwd=2)

	mtext(spname,side=1,at=seq(.2,.9,length=6),outer=TRUE,line=-2,cex=1)

  dev.off()
  
  }


Dpar.by.time.old <- function(){
  jpeg('Dpar.by.time.old.jpg',width=3.5,height=4,units='in',res=1000) 
par(mfrow=c(3,1),mar=c(2,4,1,0))
  plot(spyr[,2],gtmp[,1],pch=PT[match(spyr[,1],spall)],xaxt='n',
    ylab=expression(italic(G[ref])),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(spyr[spyr[,1]==spall[j],2],gtmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2]),gtmp[,2],as.integer(spyr[,2]),gtmp[,3])
  
  plot(spyr[,2],atmp[,1],pch=PT[match(spyr[,1],spall)],xaxt='n',
    ylab=expression(lambda),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(spyr[spyr[,1]==spall[j],2],atmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2]),atmp[,2],as.integer(spyr[,2]),atmp[,3])
  
  plot(spyr[,2],agtmp[,1],pch=PT[match(spyr[,1],spall)],xaxt='n',
    ylab=expression(lambda/italic(G[ref])),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(spyr[spyr[,1]==spall[j],2],agtmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2]),agtmp[,2],as.integer(spyr[,2]),agtmp[,3])
  legend('bottomright',legend=spall,pch=PT,cex=.9,bty='n')
  dev.off()
  
  }
VS.par <- function(){

Stmp <- Vtmp <- vatmp <- matrix(NA,dim(spyr)[1],3)

for(j in 1:dim(spyr)[1]){
  stmp <- spyr[j,1]
  ytmp <- spyr[j,2]

      tmp  <- which(splist==stmp & yrlist==ytmp) 
      tmp <- tmp[which(kplist[tmp]==tc)]
    
	if(length(tmp)==0) next
    
	Vtmp[j,]  <- quantile(vgibbs.all[[tmp]][burnin:ng,1],c(.5,.025,.975))
	Stmp[j,]  <- quantile(vgibbs.all[[tmp]][burnin:ng,2],c(.5,.025,.975))
	vatmp[j,]  <- quantile(vgibbs.all[[tmp]][burnin:ng,3],c(.5,.025,.975))
 }
  list(Vtmp = Vtmp, Stmp = Stmp, vatmp = vatmp)
  }  
  
  

QM.par <- function(){

q1tmp <- q2tmp <- m1tmp <- m2tmp <- matrix(NA,dim(spyr)[1],3)

for(j in 1:dim(spyr)[1]){
  stmp <- spyr[j,1]
  ytmp <- spyr[j,2]

      tmp  <- which(splist==stmp & yrlist==ytmp) 
      tmp <- tmp[which(kplist[tmp]==tc)]
    
	if(length(tmp)==0) next
    
	q1tmp[j,]  <- quantile(lgibbs.all[[tmp]][burnin:ng,1],c(.5,.025,.975))
	q2tmp[j,]  <- quantile(lgibbs.all[[tmp]][burnin:ng,2],c(.5,.025,.975))
	m1tmp[j,]  <- quantile(mgibbs.all[[tmp]][burnin:ng,1],c(.5,.025,.975))
	m2tmp[j,]  <- quantile(mgibbs.all[[tmp]][burnin:ng,2],c(.5,.025,.975))
  }
  list(q1tmp = q1tmp, q2tmp = q2tmp, m1tmp = m1tmp, m2tmp = m2tmp)
  }  
  

QMpar.by.time.bar <- function(){
  jpeg('QMpar.by.time.bar.jpg',width=3.5,height=5.33,units='in',res=1000) 
par(mfrow=c(4,1),mar=c(2,4,1,0))
#  plot(match(spyr[,1],spall)+spoff[match(spyr[,2],yrall)],gtmp[,1],pch=PT[match(spyr[,1],spall)],xaxt='n',
#    ylab=expression(italic(G[ref])),xlab="")#,xlim=c(min(yrall)-.5,max(yrall)+.5))

count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- q1tmp[,1]
  bx <- barplot(t(count),beside=T,col=1:4,
    ylab=expression(alpha[1]),xpd=F,ylim=c(.75,.95))
  
counthi <- matrix(NA,nsp,4)
  counthi[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- q1tmp[,3]
  arrows(bx,t(count),bx,t(counthi),angle=90,length=.05,lwd=2)

count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- q2tmp[,1]
  bx <- barplot(t(count),beside=T,col=1:4,
    ylab=expression(alpha[2]),xpd=F,ylim=c(50,550))
  
counthi <- matrix(NA,nsp,4)
  counthi[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- q2tmp[,3]
  arrows(bx,t(count),bx,t(counthi),angle=90,length=.05,lwd=2)

count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- m2tmp[,1]
  bx <- barplot(t(count),beside=T,col=1:4,
    ylab=expression(alpha[3]),xpd=F,ylim=c(.15,.45))
  
counthi <- matrix(NA,nsp,4)
  counthi[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- m2tmp[,3]
  arrows(bx,t(count),bx,t(counthi),angle=90,length=.05,lwd=2)

count <- matrix(NA,nsp,4)
  count[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- m1tmp[,1]
  bx <- barplot(t(count),beside=T,col=1:4,
    ylab=expression(alpha[4]),xpd=F,ylim=c(.02,.4))
  
counthi <- matrix(NA,nsp,4)
  counthi[cbind(match(spyr[,1],spall),match(spyr[,2],yrall))] <- m1tmp[,3]
  arrows(bx,t(count),bx,t(counthi),angle=90,length=.05,lwd=2)

  dev.off()
  
  }

  
QMpar.by.time.old <- function(){
  jpeg('QMpar.by.time.old.jpg',width=3.5,height=5.33,units='in',res=1000) 
par(mfrow=c(4,1),mar=c(3,4,1,0))
  plot(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],q1tmp[,1],
    pch=PT[match(spyr[,1],spall)],xaxt='n',ylim=c(min(q1tmp[,2]),max(q1tmp[,2])),
    ylab=expression(italic(L)[1]),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(as.integer(spyr[spyr[,1]==spall[j],2])+spoff[j],
                        q1tmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],q1tmp[,2],
           as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],q1tmp[,3])
  
  legend('topleft',legend=spall,pch=PT,cex=.9,bty='n')

  plot(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],q2tmp[,1],
    pch=PT[match(spyr[,1],spall)],xaxt='n',ylim=c(min(q2tmp[,2]),max(q2tmp[,3])),
    ylab=expression(italic(L)[2]),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(as.integer(spyr[spyr[,1]==spall[j],2])+spoff[j],
                        q2tmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],q2tmp[,2],
           as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],q2tmp[,3])
  
  plot(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],m1tmp[,1],
    pch=PT[match(spyr[,1],spall)],xaxt='n',ylim=c(min(m1lo),max(m1hi)),
    ylab=expression(italic(M)[1]),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(as.integer(spyr[spyr[,1]==spall[j],2])+spoff[j],
                        m1tmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],m1tmp[,2],
           as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],m1tmp[,3])
  
  plot(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],m2tmp[,1],
    pch=PT[match(spyr[,1],spall)],xaxt='n',ylim=c(min(m2lo),max(m2hi)),
    ylab=expression(italic(M)[2]),xlab="",xlim=c(min(yrall)-.5,max(yrall)+.5))
  axis(1,at=2002:2005,las=2)
  for(j in 1:nsp) lines(as.integer(spyr[spyr[,1]==spall[j],2])+spoff[j],
                        m2tmp[spyr[,1]==spall[j],1])
  segments(as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],m2tmp[,2],
           as.integer(spyr[,2])+spoff[match(spyr[,1],spall)],m2tmp[,3])
    
  dev.off()
  
  }
  
month.vec <- function(){
  month <- as.numeric()

  for(ytmp in yrall){
    ml <- as.integer(ytmp==2004)
    mvec  <- c(31,28+ml,31,30,31,30,31,31,30,31,30,31)
    #mvec <- (1:52)*7
	month <- c(month,(ytmp + (cumsum(mvec)-mvec[1])/sum(mvec)))
    }
  month <- c(month,ytmp+1)
  month
  }

#provide a vector of days
day.vec <- function(){
  day <- as.numeric()

  for(ytmp in yrall){
    ml <- as.integer(ytmp==2004)
    mvec  <- rep(1,length=365+ml)
    #mvec <- (1:52)*7
	day <- c(day,(ytmp + (cumsum(mvec)-mvec[1])/sum(mvec)))
    }
  day <- c(day,ytmp+1)
  day
  }

#provide a vector of days
week.vec <- function(){
  week <- as.numeric()

  for(ytmp in yrall){
    ml <- as.integer(ytmp==2004)
    #mvec  <- rep(1,length=365+ml)
    mvec <- seq(1,365+ml,by=7)
	week <- c(week,(ytmp + (cumsum(mvec)-mvec[1])/sum(mvec)))
    }
  week <- c(week,ytmp+1)
  week
  }
    
#provide a start and end to each day
day.mat <- function(){
  day <- matrix(NA,nrow=0,ncol=2)

  for(ytmp in yrall){
    ml <- as.integer(ytmp==2004)
    mvec  <- rep(1,length=365+ml)
    #mvec <- (1:52)*7
	day <- rbind(day,cbind((ytmp + (cumsum(mvec)-mvec[1])/sum(mvec)) +  8/24/(365+ml),
	                       (ytmp + (cumsum(mvec)-mvec[1])/sum(mvec)) + 17/24/(365+ml)))
    }
  #day <- rbind(day,ytmp+1)
  
  day
  }
 
  
E.G.plot <- function(){
  jpeg('trans.and.cond.jpg',width=3.5,height=7,units='in',res=1000)
  par(mfrow=c(2,1),mar=c(3,5,1,1))
  for(j in 1:length(spall)){
    if(j==1) plot(month,Et[,j],pch=PT[j],xlab="",ylab=expression(Transpiration~~(mm~~d^-1)),
	              ylim=c(min(Et,na.rm=T),max(Et,na.rm=T)+20))
	points(month,Et[,j],pch=PT[j])
	lines(month,Et[,j])
    }
	legend('topleft',legend=spname,pch=PT,bty='n',cex=.7)
  for(j in 1:length(spall)){
    if(j==1) plot(month,Gt[,j],pch=PT[j],xlab="",ylab='Canopy Conductance',ylim=c(range(Gt,na.rm=T)))
	points(month,Gt[,j],pch=PT[j])
	lines(month,Gt[,j])
    }
	dev.off()	
  }
  

E.G.plot <- function(){
  jpeg('trans.and.cond.jpg',width=6.5,height=7,units='in',res=1000)
  par(mfrow=c(3,1),mar=c(3,5,1,1))
  for(j in 1:length(spall)){
    if(j==1) plot(month,Et[,j],pch=PT[j],xlab="",ylab=expression(italic(E[C])~~(mm~~month^-1)),
	              ylim=c(min(Et,na.rm=T),max(Et,na.rm=T)+20),cex.lab=1.4,cex.axis=1.2)
	points(month,Et[,j],pch=PT[j])
	lines(month,Et[,j])
    }
	legend('topleft',legend=spname,pch=PT,bty='n',cex=1)

  for(j in 1:length(spall)){
    if(j==1) plot(month,EL[,j],pch=PT[j],xlab="",ylab=expression(italic(E[L])~~(mm~~month^-1)),
	              ylim=c(min(EL,na.rm=T),max(EL,na.rm=T)),cex.lab=1.4,cex.axis=1.2)
	points(month,EL[,j],pch=PT[j])
	lines(month,EL[,j])
    }

  for(j in 1:length(spall)){
    if(j==1) plot(month,Gt[,j],pch=PT[j],xlab="",ylab=expression(italic(G[t])~~(mmol~~ m^-2 ~~ s^-1)),
	  ylim=c(range(Gt,na.rm=T)),cex.lab=1.4,cex.axis=1.2)
	points(month,Gt[,j],pch=PT[j])
	lines(month,Gt[,j])
    }
	dev.off()	
  }
  
  
  
E.G.plot <- function(){
  jpeg('trans.and.cond.plot.jpg',width=3.5,height=7.5,units='in',res=1000)
  par(mfrow=c(3,1),mar=c(3,5,1,1))
  for(j in 1:length(spall)){
    if(j==1) plot(month,Et[,j],pch=PT[j],xlab="",cex=.8,ylab=expression(italic(E[C])~~(mm~~month^-1)),
	              ylim=c(min(Et,na.rm=T),max(Et,na.rm=T)*1.5),cex.lab=1.4,cex.axis=1.2)
	points(month,Et[,j],pch=PT[j])
	lines(month,Et[,j])
    }
	legend('topleft',legend=spname[1:3],pch=PT[1:3],bty='n',cex=1)
	legend('topright',legend=spname[4:6],pch=PT[4:6],bty='n',cex=1)

	monthkeep <- month*0+1
    monthkeep[(month-trunc(month))<.2 | (month-trunc(month))>.85] <- NA

	for(j in 1:length(spall)){
    if(j==1) plot(month,EL[,j]*monthkeep,pch=PT[j],xlab="",cex=.8,
			   ylab=expression(italic(E[L])~~(mm~~month^-1)),log='y',
	           ylim=c(min(EL,na.rm=T),max(EL,na.rm=T)),cex.lab=1.4,cex.axis=1.2)
	points(month,EL[,j]*monthkeep,pch=PT[j])
	lines(month,EL[,j]*monthkeep)
    }
	
for(j in 1:length(spall)){
    if(j==1) plot(month,Gt[,j]*monthkeep,pch=PT[j],xlab="",cex=.8,
	           ylab=expression(italic(G[t])~~(mmol ~~ m^-2 ~~ s^-1)),
			   ylim=c(range(Gt,na.rm=T)),cex.lab=1.4,cex.axis=1.2,log='y')
	points(month,Gt[,j]*monthkeep,pch=PT[j])
	lines(month,Gt[,j]*monthkeep)
    }
	dev.off()	
  }
  

E.G.pred <- function(){
  jpeg('trans.and.cond,pred.jpg',width=6,height=7.5,units='in',res=1000)
  cvec <- c(1,2,3,4,5,6)
  
  bord.col <- c(rgb(.5,.5,.5),rgb(0.9,0,0),rgb(0,0,0.9),
                rgb(.1,.1,.1),rgb(0,0.9,0),rgb(0.5,0.5,0))
  fill.col <- c(rgb(.5,.5,.5,alpha=.5),rgb(0.9,0,0,alpha=.5),rgb(0,0,0.9,alpha=.5),
                rgb(.1,.1,.1,alpha=.2),rgb(0,0.9,0,alpha=.5),rgb(0.5,0.5,0,alpha=.5))
  
  par(mfrow=c(3,2),mar=c(3,5,1,1))
  
  Etlo[Etlo==0] <- .000001

  monthkeep <- month*0+1
  monthkeep[(month-trunc(month))<.2 | (month-trunc(month))>.85] <- NA
  start.plot <- c(4,16,28,40)
  stops.plot <- c(11,23,35,47)
  noplot <- c(1:3,12:15,24:27,36:39,48:50)
  
  for(j in 1:length(spall)){
    
    if(j %in% c(1,4)) 
      plot(month,Et[,j],pch=PT[j],xlab="",cex=.8,ylab=expression(italic(E[C])~~(mm~~month^-1)),
                  ylim=c(.05,max(Ethi[-noplot,],na.rm=T)*1.5),
                  cex.lab=1.4,cex.axis=1.2,col='white',log='y')
    
    #polygon
    for(i in 1:length(start.plot))
    polygon(c(month[start.plot[i]:stops.plot[i]],rev(month[start.plot[i]:stops.plot[i]])),
            c(Etlo[start.plot[i]:stops.plot[i],j],rev(Ethi[start.plot[i]:stops.plot[i],j])),
            border=bord.col[cvec[j]],
            col=fill.col[cvec[j]])
    
    #points(month,Et[,j],pch=PT[j])
    #lines(month,Et[,j])
  
    if(j == 1) {
      legend('bottomleft',legend=spname[1:3],fill=fill.col[1:3],bty='n',cex=1)
      legend('topright',legend='(a)',bty='n')
    }
    if(j == 4) {
      legend('bottomleft',legend=spname[4:6],fill=fill.col[4:6],bty='n',cex=1)
      legend('topright',legend='(b)',bty='n')
    }
    
  }
  
  
  for(j in 1:length(spall)){
    if(j %in% c(1,4)) plot(month,EL[,j]*monthkeep,pch=PT[j],xlab="",cex=.8,
                  ylab=expression(italic(E[L])~~(mm~~month^-1)),log='y',
                  ylim=c(min(EL[-noplot,],na.rm=T),max(EL[-noplot,],na.rm=T)),cex.lab=1.4,cex.axis=1.2,col='white')
    for(i in 1:length(start.plot))
      polygon(c(month[start.plot[i]:stops.plot[i]],rev(month[start.plot[i]:stops.plot[i]])),
              c(ELlo[start.plot[i]:stops.plot[i],j],rev(ELhi[start.plot[i]:stops.plot[i],j])),
              border=bord.col[cvec[j]],
              col=fill.col[cvec[j]])
 
    if(j == 1) legend('topright',legend='(c)',bty='n')
    
    if(j == 4) legend('topright',legend='(d)',bty='n')
    
  }
  
  for(j in 1:length(spall)){
    if(j %in% c(1,4)) plot(month,Gt[,j]*monthkeep,pch=PT[j],xlab="",cex=.8,
                  ylab=expression(italic(G[t])~~(mmol ~~ m^-2 ~~ s^-1)),
                  ylim=c(range(Gt,na.rm=T)),cex.lab=1.4,cex.axis=1.2,log='y',col='white')
    for(i in 1:length(start.plot))
      polygon(c(month[start.plot[i]:stops.plot[i]],rev(month[start.plot[i]:stops.plot[i]])),
              c(Gtlo[start.plot[i]:stops.plot[i],j],rev(Gthi[start.plot[i]:stops.plot[i],j])),
              border=bord.col[cvec[j]],
              col=fill.col[cvec[j]])
 
    if(j == 1) legend('topright',legend='(e)',bty='n')
 
    if(j == 4) legend('topright',legend='(f)',bty='n')
    
  }
  dev.off()	
}


E.G.err.plot <- function(){
  jpeg('trans.and.cond.err.jpg',width=6.5,height=7,units='in',res=1000)
  par(mfrow=c(3,1),mar=c(3,5,1,1))
  for(j in 1:length(spall)){
    if(j==1) plot(month,Etlo[,j],pch=PT[j],xlab="",ylab=expression(italic(E[C])~~(mm~~month^-1)),
	              ylim=c(min(Etlo,na.rm=T),max(Ethi,na.rm=T)+20),cex.lab=1.4,cex.axis=1.2)
	points(month,Etlo[,j],pch=PT[j])
	lines(month,Etlo[,j])
	points(month,Ethi[,j],pch=PT[j])
	lines(month,Ethi[,j])
    }
	legend('topleft',legend=spname[1:3],pch=PT[1:3],bty='n',cex=1)
	legend('topright',legend=spname[4:6],pch=PT[4:6],bty='n',cex=1)

  for(j in 1:length(spall)){
    if(j==1) plot(month,ELlo[,j],pch=PT[j],xlab="",ylab=expression(italic(E[L])~~(mm~~month^-1)),
	              ylim=c(min(ELlo,na.rm=T),max(ELhi,na.rm=T)),cex.lab=1.4,cex.axis=1.2)
	points(month,ELlo[,j],pch=PT[j])
	lines(month,ELlo[,j])
	points(month,ELhi[,j],pch=PT[j])
	lines(month,ELhi[,j])
    }

  for(j in 1:length(spall)){
    if(j==1) plot(month,Gt[,j],pch=PT[j],xlab="",ylab=expression(italic(G[t])~~(mmol~~ m^-2 ~~ s^-1)),
	  ylim=c(range(Gt,na.rm=T)),cex.lab=1.4,cex.axis=1.2)
	points(month,Gt[,j],pch=PT[j])
	lines(month,Gt[,j])
    }
	dev.off()	
  }
  
getdata <- function(tt){
  wd <-  strsplit(getwd(),'/')
    wd <- unlist(wd)[length(unlist(wd))]

  setwd('../')
  setwd('datafiles')

  f   <- list.files()
  dm  <- f[grep('VPD',f)]
  mm  <- f[grep('sm',f)]
  qm  <- f[grep('PAR',f)]
  lai <- f[grep('LAI',f)]
  sap <- f[grep('SA',f)]
  tm  <- f[grep('Ta',f)]
  
	  d.vec <- read.csv(dm[tt],header=T)
	  m.vec <- read.csv(mm[tt],header=T)
	  q.vec <- read.csv(qm[tt],header=T)
	  l.vec <- read.csv(lai[tt],header=T)
	  s.vec <- read.csv(sap[tt],header=T)
	  t.vec <- read.csv(tm[tt],header=T)

	leap    <- as.numeric(d.vec[,'year']==2004)
	year    <- d.vec[,'year']
	TT      <- year + (d.vec[,'DOY']-1)/(365+leap) + (d.vec[,'Time'])/24/(365+leap)
	TTshort <- unique(year + (d.vec[,'DOY']-1)/(365+leap))
	
	d.vec <- matrix(unlist(d.vec[,-(1:4)]),length(TT))
	q.vec <- matrix(unlist(q.vec[,-(1:4)]),length(TT))
	m.vec <- matrix(unlist(m.vec[,-(1:4)]),length(TT))
	t.vec <- matrix(unlist(t.vec[,-(1:4)]),length(TT))
	s.vec <- matrix(unlist(s.vec[,-(1:4)]),length(TT))
	l.vec <- matrix(unlist(l.vec[findInterval(TT,TTshort),-(1:4)]),length(TT))
    
   
  setwd('../')
  setwd(wd)
  
  list(TT = TT, d.vec = d.vec, q.vec = q.vec, m.vec = m.vec, t.vec = t.vec, 
      s.vec = s.vec, l.vec = l.vec)
  }  
  