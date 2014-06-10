
  #wd <- '/Users/davidbell/Desktop/SF.figures'
  wd <- "/Volumes/dbell9$/Sap Flux"
  
  setwd(wd)

  source('SF.post.functions.r')


##import BLA
 yrall <- 2002:2005	
spname <- c(expression(italic(C.~~tomentosa)),
	        expression(italic(L.~~styraciflua)),
	        expression(italic(L.~~tulipifera)),
	        expression(italic(Q.~~alba)),
	        expression(italic(Q.~~michauxii)),
	        expression(italic(Q.~~phellos)))

spabbr <- c('CATO','LIST','LITU','QUAL','QUMI','QUPH')
  
#list of species
	spall <- c('cato','list','litu','qual','qumi','quph')			

#point type?
PT <- c(0,1,2,4,5,6)

tmp <- import.bla()
  BLA <- tmp$BLA     #yearly boundary layer based estimates of Gref and lambda
  spBLA <- tmp$spBLA #species identifiers for each year

##data figure file names
  wd <- 'datafiles'
  setwd(wd)

  data.figures()

  setwd('../')
  
##import runs
wd <- 'output'
setwd(wd)

ww <- list.files()
  ww<- ww[grep('min',ww)]
  
TCww <- as.integer(substr(ww,1,2)) #vector of time constants

for(W in 1:length(TCww)) {

  setwd(ww[W])
  
  ff <- list.files()
  r2 <- rep(0,length(ff))

  fyear <- matrix(unlist(strsplit(ff,'_')),length(ff),byrow=T)
    fsp   <- as.character(fyear[,7])
	fkp   <- as.numeric(fyear[,5])
    fyear <- as.numeric(fyear[,3])
  
  if(W == 1) {
    splist <- fsp   #list of species for output matrices
	yrlist <- fyear #list of years for output matrices
	kplist <- fkp
	}
  
  if(W > 1) {
    splist <- c(splist,fsp)
	yrlist <- c(yrlist,fyear)
	kplist <- c(kplist,fkp)
	}

  
  for(J in 1:length(ff)) {

    load(ff[J])
	JD <- Jdata[D>.6,]
	JP <- Jpred[D>.6,]
	r2[J] <- cor(JD[is.finite(JD)],JP[is.finite(JD)])^2

    if(J==1 & W==1){
      agibbs.all  <- list(matrix(agibbs[1:(g-1),],g-1))
      ggibbs.all  <- list(matrix(ggibbs[1:(g-1),],g-1))
      dgibbs.all  <- list(matrix(dgibbs[1:(g-1),],g-1))
      Ggibbs.all  <- list(Ggibbs)
      Ggibbs2.all <- list(Ggibbs2)
      Jgibbs.all  <- list(Jgibbs)
      Jgibbs2.all <- list(Jgibbs2)
      lgibbs.all  <- list(matrix(lgibbs[1:(g-1),],g-1))
      mgibbs.all  <- list(matrix(mgibbs[1:(g-1),],g-1))
      rgibbs.all  <- list(matrix(rgibbs[1:(g-1),],g-1))
      vgibbs.all  <- list(matrix(vgibbs[1:(g-1),],g-1))
	  Jdata.all   <- list(Jdata)
	  Jpred.all   <- list(Jpred)
	  timeall.all <- list(timeall)
      probe.all   <- list(probe)
	  gtlo.all    <- list(priorgt.lo)
      gthi.all    <- list(priorgt.hi)
      D.all     <- list(D)
      M.all     <- list(M)
      Q.all     <- list(Q)
      Temp.all  <- list(Temp)
	  time.all  <- list(timeall)
	  LAI.all   <- list(LAIspecies)
	  SAP.all   <- list(SAPspecies)
	  if('SPECVARS' %in% ls()){
	    sv.all    <- list(specvars)
	    SV.all    <- list(SPECVARS)
        }
	  }
	  
    if(J>1 | W>1){
      agibbs.all[[length(agibbs.all)+1]]  <- matrix(agibbs[1:(g-1),],g-1)
      ggibbs.all[[length(ggibbs.all)+1]]  <- matrix(ggibbs[1:(g-1),],g-1)
      dgibbs.all[[length(dgibbs.all)+1]]  <- matrix(dgibbs[1:(g-1),],g-1)
      Ggibbs.all[[length(Ggibbs.all)+1]]  <- Ggibbs
      Ggibbs2.all[[length(Ggibbs.all)+1]] <- Ggibbs2
      Jgibbs.all[[length(Jgibbs.all)+1]]  <- Jgibbs
      Jgibbs2.all[[length(Jgibbs.all)+1]] <- Jgibbs2
      lgibbs.all[[length(lgibbs.all)+1]]  <- matrix(lgibbs[1:(g-1),],g-1)
      mgibbs.all[[length(mgibbs.all)+1]]  <- matrix(mgibbs[1:(g-1),],g-1)
      rgibbs.all[[length(rgibbs.all)+1]]  <- matrix(rgibbs[1:(g-1),],g-1)
      vgibbs.all[[length(vgibbs.all)+1]]  <- matrix(vgibbs[1:(g-1),],g-1)
	  Jdata.all[[length(Jdata.all)+1]]    <- Jdata
	  Jpred.all[[length(Jpred.all)+1]]    <- Jpred
      timeall.all[[length(timeall.all)+1]]<- timeall
      probe.all[[length(probe.all)+1]]    <- probe
	  gtlo.all[[length(timeall.all)+1]]   <- priorgt.lo
      gthi.all[[length(timeall.all)+1]]   <- priorgt.hi
	  D.all[[length(D.all)+1]]       <- D
      M.all[[length(M.all)+1]]       <- M
      Q.all[[length(Q.all)+1]]       <- Q
      Temp.all[[length(Temp.all)+1]]    <- Temp
      time.all[[length(time.all)+1]]    <- timeall
      LAI.all[[length(LAI.all)+1]]     <- LAIspecies
      SAP.all[[length(SAP.all)+1]]     <- SAPspecies
	  if('SPECVARS' %in% ls()){
	    sv.all[[length(sv.all)+1]]       <- specvars
	    SV.all[[length(SV.all)+1]]       <- SPECVARS
	    }
	  }
    }
  if(W==1) R2 <- r2
  if(W>1)  R2 <- c(R2,r2)
  print(W)
  setwd('../')
  }

#species and year 
  g.all <- matrix(unlist(lapply(ggibbs.all,dim)),ncol=2,byrow=T)[,1]

#set burnin and ng for comparisons across runs
  burnin <- 2000
  ng <- min(g.all)-1

#plot convergences by species, year, and kappa
  spyr <- unique(cbind(splist,yrlist))
  converge(20)
  
## parameter response to capacitance
avec <- alo <- ahi<- gvec <- glo <- ghi <- agvec <- aglo <- aghi <- 
l1vec <- l1lo <- l1hi <- l2vec <- l2lo <- l2hi <- m1vec <- m1lo <- m1hi <- 
  m2vec <- m2lo <- m2hi <- sigvec <- siglo <- sighi <- vjvec <- vjlo <- vjhi <-
  vavec <- valo <- vahi <- yrlist*0

for(j in 1:length(avec)){  
    avec[j]    <- mean(agibbs.all[[j]][burnin:ng,]) 
    alo[j]     <- quantile(agibbs.all[[j]][burnin:ng,],.025) 
    ahi[j]     <- quantile(agibbs.all[[j]][burnin:ng,],.975) 
    gvec[j]    <- mean(ggibbs.all[[j]][burnin:ng,])
    glo[j]     <- quantile(ggibbs.all[[j]][burnin:ng,],.025) 
    ghi[j]     <- quantile(ggibbs.all[[j]][burnin:ng,],.975) 
    agvec[j]   <- mean(agibbs.all[[j]][burnin:ng,]/ggibbs.all[[j]][burnin:ng,])
    aglo[j]    <- quantile(agibbs.all[[j]][burnin:ng,]/ggibbs.all[[j]][burnin:ng,],.025) 
    aghi[j]    <- quantile(agibbs.all[[j]][burnin:ng,]/ggibbs.all[[j]][burnin:ng,],.975) 
    l1vec[j]   <- mean(lgibbs.all[[j]][burnin:ng,1])
    l1lo[j]    <- quantile(lgibbs.all[[j]][burnin:ng,1],.025) 
    l1hi[j]    <- quantile(lgibbs.all[[j]][burnin:ng,1],.975) 
    l2vec[j]  <- mean(lgibbs.all[[j]][burnin:ng,2])
    l2lo[j]    <- quantile(lgibbs.all[[j]][burnin:ng,2],.025) 
    l2hi[j]    <- quantile(lgibbs.all[[j]][burnin:ng,2],.975) 
    m1vec[j]  <- mean(mgibbs.all[[j]][burnin:ng,1])
    m1lo[j]    <- quantile(mgibbs.all[[j]][burnin:ng,1],.025) 
    m1hi[j]    <- quantile(mgibbs.all[[j]][burnin:ng,1],.975) 
    m2vec[j]  <- mean(mgibbs.all[[j]][burnin:ng,2])
    m2lo[j]    <- quantile(mgibbs.all[[j]][burnin:ng,2],.025) 
    m2hi[j]    <- quantile(mgibbs.all[[j]][burnin:ng,2],.975) 
    sigvec[j] <- mean(vgibbs.all[[j]][burnin:ng,1])
    siglo[j]    <- quantile(vgibbs.all[[j]][burnin:ng,1],.025) 
    sighi[j]    <- quantile(vgibbs.all[[j]][burnin:ng,1],.975) 
    vjvec[j]  <- mean(vgibbs.all[[j]][burnin:ng,2])
    vjlo[j]    <- quantile(vgibbs.all[[j]][burnin:ng,2],.025) 
    vjhi[j]    <- quantile(vgibbs.all[[j]][burnin:ng,2],.975) 
    vavec[j]  <- mean(vgibbs.all[[j]][burnin:ng,3])
    valo[j]    <- quantile(vgibbs.all[[j]][burnin:ng,3],.025) 
    vahi[j]    <- quantile(vgibbs.all[[j]][burnin:ng,3],.975) 
	}
	
spname <- c(expression(italic(C.~~tomentosa)),
	        expression(italic(L.~~styraciflua)),
	        expression(italic(L.~~tulipifera)),
	        expression(italic(Q.~~alba)),
	        expression(italic(Q.~~michauxii)),
	        expression(italic(Q.~~phellos)))
spall <- c('cato','list','litu','qual','qumi','quph')			

kpall <- sort(unique(kplist))
 yrall <- 2002:2005	

##list of iterations to sample
  sample <- seq(1,nrow(agibbs.all[[1]]),by=100)
  
## pred vs. observed -- TC
r2 <- pred.vs.obs.TC()

##another pred vs obs
exp.var.sp()
  
## pred vs. observed -- examples
     sp <- c('qual','list','cato')
     yr <- c(2004,2005,2005)
  probe <- c(3, 2, 3)
     tc <- 1 
	 nsamp <- 1000
     panname <- c('(a)','(b)','(c)')
	 
	 lim <- c(0,70)
	 
  postscript('obs.pred.ex.eps',width=6.5,height=2,onefile=FALSE,horizontal=FALSE)	
  par(mfrow=c(1,3),mar=c(5,5,1,1))
  for(j in 1:length(sp)){
    tmp <- which(splist==sp[j] & yrlist==yr[j] & kplist==tc)
    samp <- sample(1:dim(Jdata.all[[tmp]])[1],nsamp)
	#lim <- range(c(Jdata.all[[tmp]][samp,probe[j]],Jpred.all[[tmp]][samp,probe[j]]),na.rm=T)
    plot(Jdata.all[[tmp]][samp,probe[j]],Jpred.all[[tmp]][samp,probe[j]],
	  xlab=expression(Observed~~italic(J[it])),cex.lab=1.2,
	  ylab=expression(Predicted~~italic(J[it])),ylim=lim,xlim=lim,cex=.75)
	legend('topleft',legend=panname[j],bty='n')
	abline(0,1,col='darkgray',lty=2,lwd=2)
    }
	dev.off()
	
#explained variation
exp.var.ind()

#error biplot
err.biplot()

#sp biplots
for(stmp in spall)
  sp.biplot(stmp)

#pararmeter by kappa
par.by.kappa()

##dbh by random effects
dbh <- get.dbh()

##
hist.exp()

##
jpeg('rand.by.diam.jpg',width=7,height=5,units='in',res=1000)
tc <- 1
par(mfrow=c(2,3),mar=c(5,5,1,1),family='serif')
for(stmp in spall){
   l <- 1
  for(ytmp in yrall){
    tmp <- which(splist==stmp & yrlist==ytmp & kplist==tc)
	if(length(tmp)==0) next
		
	db <- dbh[dbh[,1] %in% probe.all[[tmp]],]
	  db <- db[match(db[,1],probe.all[[tmp]]),1+which(yrall==ytmp)]
	rb <- apply(rgibbs.all[[tmp]],2,mean)
	
	if(l==1) plot(db,rb,xlab='diameter (cm)',ylab=expression(Random~~Effect~~italic(a[ij])),
	  main=spname[spall==stmp],pch=PT[spall==stmp],col=which(yrall==ytmp),cex=1.2,
	  xlim=range(dbh[grep(stmp,dbh[,1]),-1]),ylim=c(.5,1.8),cex.lab=1.4,cex.axis=1.2)
	  points(db,rb,pch=PT[spall==stmp],col=which(yrall==ytmp),cex=1.2)
	l <- 0
    }
  }
dev.off()  

##biplot Gref and lambda
yl <- expression(lambda)
xl <- expression(italic(G[ref]))
tc <- 1

tmp <- gref.lam.rat.jpg()
atmp  <- tmp$atmp
gtmp  <- tmp$gtmp
agtmp <- tmp$agtmp

## 
gref.lambda() 
##
##response functions

Qmax <- tapply(unlist(lapply(Q.all,max)),paste(yrlist,splist),max)
Qmin <- tapply(unlist(lapply(Q.all,min)),paste(yrlist,splist),min)
Mmax <- tapply(unlist(lapply(M.all,max)),paste(yrlist,splist),max)
Mmin <- tapply(unlist(lapply(M.all,min)),paste(yrlist,splist),min)
Dmax <- tapply(unlist(lapply(D.all,max)),paste(yrlist,splist),max)
Dmin <- tapply(unlist(lapply(D.all,min)),paste(yrlist,splist),min)*0+.6

  rlength <- 100
  Qm <- Mm <- Dm <- matrix(NA,rlength,length(Qmax))

  for(j in 1:length(Qmax)){
    Qm[,j] <- exp(seq(log(max(Qmin[j],1)),log(Qmax[j]),length=rlength))
    Mm[,j] <- exp(seq(log(Mmin[j]),log(Mmax[j]),length=rlength))
    Dm[,j] <- exp(seq(log(Dmin[j]),log(Dmax[j]),length=rlength))
  }
  
tmp <- resp()
  Qvec <- tmp$Qvec
  Mvec <- tmp$Mvec
  Dvec <- tmp$Dvec
  
##plotting relations
##four panels, one for each year
  resp.by.year.plot()

##six panels, one for each species
 # resp.by.species.plot()
  
##combined responses with covariate as columns and species as rows and four years plotted within each
plt.mat <- rbind(c(.10,.35,.81,.95),
                 c(.45,.70,.81,.95),
                 c(.70,.95,.81,.95),
                 c(.10,.35,.67,.81),
                 c(.45,.70,.67,.81),
                 c(.70,.95,.67,.81),
                 c(.10,.35,.53,.67),
                 c(.45,.70,.53,.67),
                 c(.70,.95,.53,.67),
                 c(.10,.35,.39,.53),
                 c(.45,.70,.39,.53),
                 c(.70,.95,.39,.53),
                 c(.10,.35,.25,.39),
                 c(.45,.70,.25,.39),
                 c(.70,.95,.25,.39),
                 c(.10,.35,.11,.25),
                 c(.45,.70,.11,.25),
                 c(.70,.95,.11,.25))
  resp.by.species.plot(w=3.5,h=6,cc = .65)
  
##one panel, use symbols and colors from previous panels
nsp <- length(spall) 
spoff <- seq(-.4,.4,length=nsp)
  Dpar.by.time.bar()

  tc <- 1	
  tmp <- QM.par()
    q1tmp <- tmp$q1tmp
	q2tmp <- tmp$q2tmp
	m1tmp <- tmp$m1tmp
	m2tmp <- tmp$m2tmp

  QMpar.by.time.bar()

 tmp <- VS.par()
 	Vtmp <- tmp$Vtmp
 	Stmp <- tmp$Stmp
 	vatmp <- tmp$vatmp
 	

##
tc <- 1

	  
  samp  <- sample(burnin:ng,100,replace=T)
  l <- 1
  
  err <- 0 #zero to not include, one to include
  
	MCUT <- matrix(NA,length(samp),length(yrlist[kplist==tc]))

 for(ytmp in yrall){
    for(stmp in spall){
      print(paste(stmp,ytmp))
	  
	  tmp  <- which(yrlist==ytmp & splist==stmp & kplist == tc)
	  tmp2 <- which(yrlist[kplist==tc]==ytmp & splist[kplist==tc]==stmp)
	  
	  if(length(tmp)==0) next
      
	  tmpdata <- getdata(which(yrall==ytmp))
        timeall     <- tmpdata$TT
        D           <- tmpdata$d.vec
        M           <- tmpdata$m.vec
        Q           <- tmpdata$q.vec
        LAIspecies  <- tmpdata$l.vec
        SAPspecies  <- tmpdata$s.vec
        Temp        <- tmpdata$t.vec

	    nt <- length(timeall)
	    dt <- diff(timeall)

		gkeep <- gmdkeep <- matrix(NA,nt,length(samp))
	    
      for(j in 1:length(samp)){
		gspec      <- ggibbs.all[[tmp]][samp[j],1]
		aspec      <- agibbs.all[[tmp]][j,1]
		blitespec  <- matrix(lgibbs.all[[tmp]][samp[j],],2)
		bmoistspec <- matrix(mgibbs.all[[tmp]][samp[j],],2)
		sigma      <- vgibbs.all[[tmp]][samp[j],1]
		
		MCUT[j,tmp2] <- bmoistspec[2]
		
		Gsmat      <- gssMat(gspec,aspec,blitespec,bmoistspec)
		Gtmat      <- Gsmat
		for(t in 2:dim(Gtmat)[1]) 
			Gtmat[t,] <- glagMat(Gtmat[t-1,],Gsmat[t,],dt[t-1],tau) +
			             rnorm(1,0,sqrt(sigma))*err

		Q          <- Q * 0 + 2000 
		GMDmat     <- gssMat(gspec,aspec,blitespec,bmoistspec)
		Q           <- tmpdata$q.vec
		
		gkeep[,j] <- Gtmat
		gmdkeep[,j] <- GMDmat
        }
		
      if(l==1) {
      	gtemp <- list(gkeep)
      	gmdtemp <- list(gmdkeep)
      	}
	  if(l==0) {
	  	gtemp[[length(gtemp)+1]] <- gkeep
	  	gmdtemp[[length(gmdtemp)+1]] <- gmdkeep
      	}
	  l<-0
      }
	
    }

	month <- month.vec()

	Et <- EL <- Gt <- Gmd <- Lt <- St <- matrix(NA,length(month),nsp)
	Etlo <- ELlo <- Gtlo <- Gmdlo  <- matrix(NA,length(month),nsp)
	Ethi <- ELhi <- Gthi <- Gmdhi  <- matrix(NA,length(month),nsp)
	M.Mon <- Q.Mon <- D.Mon <- matrix(NA,length(month),ncol=1)
  
   
  for(ytmp in yrall){
    for(stmp in spall){
	  tmp  <- which(yrlist == ytmp & splist == stmp & kplist==tc)
      if(length(tmp)==0) next
	  
	  print(paste(stmp,ytmp))
	  
	  tmpdata <- getdata(which(yrall==ytmp))
        timeall     <- tmpdata$TT
        D           <- tmpdata$d.vec
        M           <- tmpdata$m.vec
        Q           <- tmpdata$q.vec
        LAIspecies  <- tmpdata$l.vec
        SAPspecies  <- tmpdata$s.vec
        Temp        <- tmpdata$t.vec
	    if('SPECVARS' %in% ls()){
		  specvars    <- sv.all[[tmp]]
		  SPECVARS    <- SV.all[[tmp]]
          }
	    nt <- length(timeall)
	    dt <- diff(timeall)
	  
	  tind <- findInterval(timeall,month)
	  
    ghold <- gtemp[[which(spyr[,1]==stmp & spyr[,2]==ytmp)]]
    ghold[ghold<0] <- 0
    gmdhold <- gmdtemp[[which(spyr[,1]==stmp & spyr[,2]==ytmp)]]
    gmdhold[ghold<0] <- 0
    
    	gt <- gmd <- eC <- eL <- matrix(NA,length(unique(tind)),ncol(ghold))

	  	e.qt <- matrix(qt(LAIspecies[,which(spall %in% stmp)],
	            		  SAPspecies[,which(spall %in% stmp)],D[,site.all[,1]],
						  Temp[,site.all[,1]]),nt)

    	for(j in 1:ncol(gt)) {
    	  gtmp <- ghold
    	  gtmp[Q==0,] <- NA
    	  gtmp   <- tapply(gtmp[,j],tind,mean,na.rm=TRUE) 
    	  gt[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  
    	  gtmp <- gmdhold
    	  gtmp[Q==0,] <- NA
    	  gtmp   <- tapply(gtmp[,j],tind,mean,na.rm=TRUE) 
    	  gmd[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  
    	  
    	  #canopy transpiration EC = G * e.qt ( S*ratio)
    	  gtmp   <- tapply(ghold[,j]*
    	                  e.qt*SAPspecies[,which(spall %in% stmp)]*
    	                  specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           				  SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']],
           				  tind,mean,na.rm=TRUE)
    	  eC[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp 
    	  
    	  #leaf transpiration EL = G * e.qt ( S*ratio) / L
    	  gtmp   <- tapply(ghold[,j]*
    	                  e.qt*SAPspecies[,which(spall %in% stmp)]*
    	                  specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           				  SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
    	                  LAIspecies[,which(spall %in% stmp)],
    	                  tind,mean,na.rm=TRUE) 
    	  eL[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  }
      
      	rm(list=c('gtmp'))	
	  
	  
      Gt[unique(tind),which(spall==stmp)]     <- apply(gt,1,quantile,.5,na.rm=T)
	    Gtlo[unique(tind),which(spall==stmp)] <- apply(gt,1,quantile,.025,na.rm=T)
	    Gthi[unique(tind),which(spall==stmp)] <- apply(gt,1,quantile,.975,na.rm=T)

      Gmd[unique(tind),which(spall==stmp)]     <- apply(gmd,1,quantile,.5,na.rm=T)
	    Gmdlo[unique(tind),which(spall==stmp)] <- apply(gmd,1,quantile,.025,na.rm=T)
	    Gmdhi[unique(tind),which(spall==stmp)] <- apply(gmd,1,quantile,.975,na.rm=T)
	  
	  Et[unique(tind),which(spall==stmp)]   <- apply(eC,1,quantile,.5,na.rm=T)*60*60*24*30/1000
	  	Etlo[unique(tind),which(spall==stmp)] <- apply(eC,1,quantile,.025,na.rm=T)*60*60*24*30/1000
	  	Ethi[unique(tind),which(spall==stmp)] <- apply(eC,1,quantile,.975,na.rm=T)*60*60*24*30/1000
	  
	  EL[unique(tind),which(spall==stmp)]   <- apply(eL,1,quantile,.5,na.rm=T)*60*60*24*30/1000
	  	ELlo[unique(tind),which(spall==stmp)] <- apply(eL,1,quantile,.025,na.rm=T)*60*60*24*30/1000
	  	ELhi[unique(tind),which(spall==stmp)] <- apply(eL,1,quantile,.975,na.rm=T)*60*60*24*30/1000
	  
	  Lt[unique(tind),which(spall==stmp)]   <- tapply(LAIspecies[,spall==stmp],tind,quantile,.5,na.rm=T)
	  St[unique(tind),which(spall==stmp)]   <- tapply(SAPspecies[,spall==stmp],tind,quantile,.5,na.rm=T)

	  M.Mon[unique(tind),] <- tapply(M,tind,quantile,.5,na.rm=T)
	  Q.Mon[unique(tind),] <- tapply(Q,tind,quantile,.5,na.rm=T)
	  D.Mon[unique(tind),] <- tapply(D,tind,quantile,.5,na.rm=T)

	  }
	}


E.G.plot()


E.G.pred()
  
  
#synthetic drought response

	day <- day.vec()	#day times
	
	Et.daily <- EL.daily <- Gt.daily <- Gmd.daily<- Lt.daily <- St.daily <- matrix(NA,length(day),nsp)
	Etlo.daily <- ELlo.daily <- Gtlo.daily <- Gmdlo.daily <- matrix(NA,length(day),nsp)
	Ethi.daily <- ELhi.daily <- Gthi.daily <- Gmdhi.daily <- matrix(NA,length(day),nsp)
	M.daily <- Q.daily <- D.daily <- matrix(NA,length(day),ncol=1)


 for(ytmp in yrall){
    for(stmp in spall){
	  tmp  <- which(yrlist == ytmp & splist == stmp & kplist==tc)
      if(length(tmp)==0) next
	  
	  print(paste(stmp,ytmp))
	  
	  tmpdata <- getdata(which(yrall==ytmp))
        timeall     <- tmpdata$TT
        D           <- tmpdata$d.vec
        M           <- tmpdata$m.vec
        Q           <- tmpdata$q.vec
        LAIspecies  <- tmpdata$l.vec
        SAPspecies  <- tmpdata$s.vec
        Temp        <- tmpdata$t.vec
	    if('SPECVARS' %in% ls()){
		  specvars    <- sv.all[[tmp]]
		  SPECVARS    <- SV.all[[tmp]]
          }
	    nt <- length(timeall)
	    dt <- diff(timeall)
	  
	  tind <- findInterval(timeall,day)
	   
    ghold <- gtemp[[which(spyr[,1]==stmp & spyr[,2]==ytmp)]]
    ghold[ghold<0] <- 0
    gmdhold <- gmdtemp[[which(spyr[,1]==stmp & spyr[,2]==ytmp)]]
    gmdhold[ghold<0] <- 0
    
    #get daily predictions for each set of parameters
    	gt <- gmd <- eC <- eL <- matrix(NA,length(unique(tind)),ncol(ghold))

	  	e.qt <- matrix(qt(LAIspecies[,which(spall %in% stmp)],
	            		  SAPspecies[,which(spall %in% stmp)],D[,site.all[,1]],
						  Temp[,site.all[,1]]),nt)

    	for(j in 1:ncol(gt)) {
    	  gtmp <- ghold
    	  gtmp[Q==0,] <- NA
    	  gtmp   <- tapply(gtmp[,j],tind,mean,na.rm=TRUE) 
    	  gt[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp

    	  gtmp <- gmdhold
    	  gtmp[Q==0,] <- NA
    	  gtmp   <- tapply(gtmp[,j],tind,mean,na.rm=TRUE) 
    	  gmd[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  
    	  gtmp   <- tapply(ghold[,j]*
    	                  e.qt*SAPspecies[,which(spall %in% stmp)]*
    	                  specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           				  SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']],
           				  tind,mean,na.rm=TRUE) 
    	  eC[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  
    	  gtmp   <- tapply(ghold[,j]*
    	                  e.qt*SAPspecies[,which(spall %in% stmp)]*
    	                  specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           				  SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
    	                  LAIspecies[,which(spall %in% stmp)],tind,mean,na.rm=TRUE) 
    	  eL[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  }
      
      	rm(list=c('gtmp'))	
	  
	  
      Gt.daily[unique(tind),which(spall==stmp)]     <- apply(gt,1,quantile,.5,na.rm=T)
	    Gtlo.daily[unique(tind),which(spall==stmp)] <- apply(gt,1,quantile,.025,na.rm=T)
	    Gthi.daily[unique(tind),which(spall==stmp)] <- apply(gt,1,quantile,.975,na.rm=T)
	  
      Gmd.daily[unique(tind),which(spall==stmp)]     <- apply(gmd,1,quantile,.5,na.rm=T)
	    Gmdlo.daily[unique(tind),which(spall==stmp)] <- apply(gmd,1,quantile,.025,na.rm=T)
	    Gmdhi.daily[unique(tind),which(spall==stmp)] <- apply(gmd,1,quantile,.975,na.rm=T)

	  Et.daily[unique(tind),which(spall==stmp)]   <- apply(eC,1,quantile,.5,na.rm=T)*60*60*24*30/1000
	  Etlo.daily[unique(tind),which(spall==stmp)] <- apply(eC,1,quantile,.025,na.rm=T)*60*60*24*30/1000
	  Ethi.daily[unique(tind),which(spall==stmp)] <- apply(eC,1,quantile,.975,na.rm=T)*60*60*24*30/1000
	  
	  EL.daily[unique(tind),which(spall==stmp)]   <- apply(eL,1,quantile,.5,na.rm=T)*60*60*24*30/1000
	  ELlo.daily[unique(tind),which(spall==stmp)] <- apply(eL,1,quantile,.025,na.rm=T)*60*60*24*30/1000
	  ELhi.daily[unique(tind),which(spall==stmp)] <- apply(eL,1,quantile,.975,na.rm=T)*60*60*24*30/1000
	  
	  Lt.daily[unique(tind),which(spall==stmp)]   <- tapply(LAIspecies[,spall==stmp],tind,quantile,.5,na.rm=T)
	  St.daily[unique(tind),which(spall==stmp)]   <- tapply(SAPspecies[,spall==stmp],tind,quantile,.5,na.rm=T)
	  
	  M.daily[unique(tind),] <- tapply(M,tind,quantile,.5,na.rm=T)
	  Q.daily[unique(tind),] <- tapply(Q,tind,quantile,.5,na.rm=T)
	  D.daily[unique(tind),] <- tapply(D,tind,quantile,.5,na.rm=T)
	  }
	}

#weekly numbers

	week <- week.vec()	#week times
	yr.tmp <- match(trunc(week),yrall)
	
	Et.week <- EL.week <- Gt.week <- Gmd.week<- Lt.week <- St.week <- matrix(NA,length(week),nsp)
	Etlo.week <- ELlo.week <- Gtlo.week <- Gmdlo.week <- matrix(NA,length(week),nsp)
	Ethi.week <- ELhi.week <- Gthi.week <- Gmdhi.week <- matrix(NA,length(week),nsp)
	M.week <- Q.week <- D.week <- matrix(NA,length(week),ncol=1)

	B.DROUGHT <- array(NA,dim=c(2,length(yrlist),4))
		dimnames(B.DROUGHT) <- list(c('B0','B1'),paste(splist,yrlist),c('mean','sd','2.5','9.5'))
		BMD.DROUGHT <- B.DROUGHT
		
	samp  <- sample(burnin:ng,100,replace=T)
	l <- 1

 for(ytmp in yrall){
    for(stmp in spall){

	  tmp  <- which(yrlist == ytmp & splist == stmp & kplist==tc)

 	  tmp2 <- which(yrlist[kplist == tc]==ytmp & splist[kplist == tc]==stmp)

      if(length(tmp)==0) next
	  
	  print(paste(stmp,ytmp))
	  
	  tmpdata <- getdata(which(yrall==ytmp))
        timeall     <- tmpdata$TT
        D           <- tmpdata$d.vec
        M           <- tmpdata$m.vec
        Q           <- tmpdata$q.vec
        LAIspecies  <- tmpdata$l.vec
        SAPspecies  <- tmpdata$s.vec
        Temp        <- tmpdata$t.vec
	    if('SPECVARS' %in% ls()){
		  specvars    <- sv.all[[tmp]]
		  SPECVARS    <- SV.all[[tmp]]
          }
	    nt <- length(timeall)
	    dt <- diff(timeall)
	  
	  tind <- findInterval(timeall,week)
	   
    ghold <- gtemp[[which(spyr[,1]==stmp & spyr[,2]==ytmp)]]
    ghold[ghold<0] <- 0
    gmdhold <- gmdtemp[[which(spyr[,1]==stmp & spyr[,2]==ytmp)]]
    gmdhold[ghold<0] <- 0
    
    #get daily predictions for each set of parameters
    	gt <- gmd <- eC <- eL <- matrix(NA,length(unique(tind)),ncol(ghold))

	  	e.qt <- matrix(qt(LAIspecies[,which(spall %in% stmp)],
	            		  SAPspecies[,which(spall %in% stmp)],D[,site.all[,1]],
						  Temp[,site.all[,1]]),nt)

    	for(j in 1:ncol(gt)) {
    	  gtmp <- ghold
    	  gtmp[Q==0,] <- NA
    	  gtmp   <- tapply(gtmp[,j],tind,mean,na.rm=TRUE) 
    	  gt[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp

    	  gtmp <- gmdhold
    	  gtmp[Q==0,] <- NA
    	  gtmp   <- tapply(gtmp[,j],tind,mean,na.rm=TRUE) 
    	  gmd[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  
    	  
    	  gtmp   <- tapply(ghold[,j]*
    	                  e.qt*SAPspecies[,which(spall %in% stmp)]*
    	                  specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           				  SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']],
           				  tind,mean,na.rm=TRUE) 
    	  eC[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  
    	  gtmp   <- tapply(ghold[,j]*
    	                  e.qt*SAPspecies[,which(spall %in% stmp)]*
    	                  specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           				  SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
    	                  LAIspecies[,which(spall %in% stmp)],tind,mean,na.rm=TRUE) 
    	  eL[match(as.integer(names(gtmp)),sort(unique(tind))),j] <- gtmp
    	  }
      
      	rm(list=c('gtmp'))	
      	
	 
	 #save conductance 
     Gt.week[unique(tind),which(spall==stmp)]     <- apply(gt,1,quantile,.5,na.rm=T)
	    Gtlo.week[unique(tind),which(spall==stmp)] <- apply(gt,1,quantile,.025,na.rm=T)
	    Gthi.week[unique(tind),which(spall==stmp)] <- apply(gt,1,quantile,.975,na.rm=T)
	  
      Gmd.week[unique(tind),which(spall==stmp)]     <- apply(gmd,1,quantile,.5,na.rm=T)
	    Gmdlo.week[unique(tind),which(spall==stmp)] <- apply(gmd,1,quantile,.025,na.rm=T)
	    Gmdhi.week[unique(tind),which(spall==stmp)] <- apply(gmd,1,quantile,.975,na.rm=T)

	  Et.week[unique(tind),which(spall==stmp)]   <- apply(eC,1,quantile,.5,na.rm=T)*60*60*24*30/1000
	  Etlo.week[unique(tind),which(spall==stmp)] <- apply(eC,1,quantile,.025,na.rm=T)*60*60*24*30/1000
	  Ethi.week[unique(tind),which(spall==stmp)] <- apply(eC,1,quantile,.975,na.rm=T)*60*60*24*30/1000
	  
	  EL.week[unique(tind),which(spall==stmp)]   <- apply(eL,1,quantile,.5,na.rm=T)*60*60*24*30/1000
	  ELlo.week[unique(tind),which(spall==stmp)] <- apply(eL,1,quantile,.025,na.rm=T)*60*60*24*30/1000
	  ELhi.week[unique(tind),which(spall==stmp)] <- apply(eL,1,quantile,.975,na.rm=T)*60*60*24*30/1000
	  
	  Lt.week[unique(tind),which(spall==stmp)]   <- tapply(LAIspecies[,spall==stmp],tind,quantile,.5,na.rm=T)
	  St.week[unique(tind),which(spall==stmp)]   <- tapply(SAPspecies[,spall==stmp],tind,quantile,.5,na.rm=T)
	  
	  M.week[unique(tind),] <- tapply(M,tind,quantile,.5,na.rm=T)
	  Q.week[unique(tind),] <- tapply(Q,tind,quantile,.5,na.rm=T)
	  D.week[unique(tind),] <- tapply(D,tind,quantile,.5,na.rm=T)
	  

     #regression on drought conductance
     
     		kk <- yrall[yr.tmp]==ytmp & is.finite(yrall[yr.tmp]==ytmp)
     		xm <- cbind(1,M.week[kk][M.week[kk] <= max(apply(MCUT,2,mean)[splist[kplist==tc] == stmp])])
     		ym <- gt[M.week[kk] <= max(apply(MCUT,2,mean)[splist[kplist==tc] == stmp]),]
     		
     		print(nrow(ym))
     		
     		if(nrow(ym)>=4) {
     		
     		B.DROUGHT[,tmp,1] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,mean)
     		B.DROUGHT[,tmp,2] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,sd)
     		B.DROUGHT[,tmp,3] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,quantile,.025)
     		B.DROUGHT[,tmp,4] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,quantile,.975)
     		
     		ym <- gmd[M.week[kk] <= max(apply(MCUT,2,mean)[splist[kplist==tc] == stmp]),]
     		
     		BMD.DROUGHT[,tmp,1] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,mean)
     		BMD.DROUGHT[,tmp,2] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,sd)
     		BMD.DROUGHT[,tmp,3] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,quantile,.025)
     		BMD.DROUGHT[,tmp,4] <- apply(solve(crossprod(xm)) %*% crossprod(xm,ym),1,quantile,.975)
     		}
   	
	  }
	}

#plot conductance response

monthkeep <- which((month - trunc(month))>.2 & (month - trunc(month))<.85)

par(mfrow=c(6,2),mar=c(4,4,1,1))

for(j in 1:length(spall)){
	
	plot(D.Mon[monthkeep],Gt[monthkeep,j],xlab='VPD',ylab='f(D)g(M)')
	plot(M.Mon[monthkeep],Gt[monthkeep,j],xlab='Moisture',ylab='f(D)g(M)')
	
	
	}

daykeep <- which((day - trunc(day))>.2 & (day - trunc(day))<.85)
weekkeep <- which((week - trunc(week))>.2 & (week - trunc(week))<.85)

plt.mat <- rbind(c(.20,.55,.81,.95),
                 c(.57,.92,.81,.95),
                 c(.20,.55,.67,.81),
                 c(.57,.92,.67,.81),
                 c(.20,.55,.53,.67),
                 c(.57,.92,.53,.67),
                 c(.20,.55,.39,.53),
                 c(.57,.92,.39,.53),
                 c(.20,.55,.25,.39),
                 c(.57,.92,.25,.39),
                 c(.20,.55,.11,.25),
                 c(.57,.92,.11,.25))


D.time <- D.week
M.time <- M.week
G.time <- Gmd.week
keep <- weekkeep


postscript('Drought.week.eps',width=3,height=7,onefile=FALSE,horizontal=FALSE)

plot.new()

let <- c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)')

for(j in 1:length(spall)){
	
    par(new = "TRUE",plt = plt.mat[(j-1)*2 +1,],las = 1,cex.axis = 1)

	plot(D.time[keep],G.time[keep,j],xlab='',ylab='',cex=.6,
	     col=yr.tmp[keep], pch = c(0,15,1,16)[yr.tmp[keep]],axes=FALSE,frame.plot=TRUE)
	     
	     legend(0,max(G.time[keep,j],na.rm=TRUE)*1.05,legend=let[j],cex=.7,bty='n')
	     
	axis(2,cex.axis=.5)
	if(j==6) {
		axis(1,cex.axis=.5)
	    mtext(expression(Vap. ~~ Pres. ~~ Def. ~~ italic(D[t])),side=1,cex=.7,line=2)
		}	
	if(j==1) legend(.77,50,legend=yrall,pch=c(0,15,1,16),bty="n",cex=.55,col=1:4)

    par(new = "TRUE",plt = plt.mat[(j-1)*2 +2,],las = 1,cex.axis = 1)

	plot(M.time[keep],G.time[keep,j],xlab='',ylab='',cex=.6,
	     col=yr.tmp[keep],pch = c(0,15,1,16)[yr.tmp[keep]],axes=FALSE,frame.plot=TRUE)
	
	     legend(.1,max(G.time[keep,j],na.rm=TRUE)*1.05,legend=let[j+6],cex=.7,bty='n')
	
	if(j==6) 	{
	  axis(1,cex.axis=.5)
	  mtext(expression(Soil ~~ Moist. ~~ italic(M[t])),side=1,cex=.7,line=2)
	  }

	#mtext(expression(italic(f)*(italic(D[t]))*italic(h)*(italic(M[t]))),side=2,outer=TRUE,
	#  line=-1.5,las=3,cex=.7)
	 
	if(j==1)
		mtext(expression(italic(f)*(italic(D[t]))*italic(h)*(italic(M[t]))),side=2,outer=TRUE,
	  line=-1.5,las=3,cex=.7)
	  
	mtext(spname[j],side=4,las=3,cex=.65)
	
	
	for(ytmp in 1:length(yrall)){
	
		if(TRUE){
			xm <- cbind(1,M.time[keep][M.time[keep] <= max(apply(MCUT,2,mean)[splist[kplist==tc]==spall[j]]) &
		                  yr.tmp[keep]==ytmp])
	    	ym <- G.time[keep,j][M.time[keep] <= max(apply(MCUT,2,mean)[splist[kplist==tc]==spall[j]]) &
		                  yr.tmp[keep]==ytmp]
		    }
	if(FALSE){
			xm <- cbind(1,M.time[keep][M.time[keep] <= .3 & yr.tmp[keep]==ytmp])
	    	ym <- G.time[keep,j][M.time[keep] <= .3 & yr.tmp[keep]==ytmp]
		    }

		#abline(v=mean(MCUT[,splist==spall[j] & yrlist == yrall[ytmp]]),lty=2,col=ytmp)
	
		if(length(ym)<4) next
	
		B <- solve(crossprod(xm))%*%crossprod(xm,ym)
		
		lines(xm[,2],xm%*%B,col=ytmp,lwd=3)
		
	}
	
	}
	dev.off()
	
yr.tmp <- match(trunc(month),yrall)
	
postscript('Drought.Month.eps',width=3,height=7,onefile=FALSE,horizontal=FALSE)

plot.new()


for(j in 1:length(spall)){
	
    par(new = "TRUE",plt = plt.mat[(j-1)*2 +1,],las = 1,cex.axis = 1)

	plot(D.Mon[monthkeep],Gt[monthkeep,j],xlab='',ylab='',cex=.3,
	     col = yr.tmp[monthkeep],axes=FALSE,frame.plot=TRUE)
	     
	axis(2,cex.axis=.5)
	if(j==6) {
		axis(1,cex.axis=.5)
	    mtext(expression(Vap. ~~ Pres. ~~ Def. ~~ italic(D[t])),side=1,cex=.7,line=2)
		}	
	if(j==1) legend('topright',legend=yrall,text.col=1:4,bty="n",cex=.65)

    par(new = "TRUE",plt = plt.mat[(j-1)*2 +2,],las = 1,cex.axis = 1)

	plot(M.Mon[monthkeep],Gt[monthkeep,j],xlab='',ylab='',cex=.3,
	     col = yr.tmp[monthkeep],axes=FALSE,frame.plot=TRUE)
	
	
	if(j==6) 	{
	  axis(1,cex.axis=.5)
	  mtext(expression(Soil ~~ Moist. ~~ italic(M[t])),side=1,cex=.7,line=2)
	  }

	#mtext(expression(italic(f)*(italic(D[t]))*italic(h)*(italic(M[t]))),side=2,outer=TRUE,
	#  line=-1.5,las=3,cex=.7)
	 
	mtext(expression(Can. ~~ Cond. ~~ G[t]),side=2,outer=TRUE,
	  line=-1.5,las=3,cex=.7)
	  
	mtext(spname[j],side=4,las=3,cex=.65)
	
	
	}
	dev.off()
	
	
postscript('Drought.Interval.eps',width=3,height=7,onefile=FALSE,horizontal=FALSE)

plot.new()

D.vec <- seq(0,2,by=.5)
M.vec <- seq(.1,.6,by=.1)

D.dayInt <- findInterval(D.daily,D.vec)
M.dayInt <- findInterval(M.daily,M.vec)

for(j in 1:length(spall)){
	
    par(new = "TRUE",plt = plt.mat[(j-1)*2 +1,],las = 1,cex.axis = 1)
	
	pl <- 1
	
	for(yy in 1:length(yrall)){
	
		if(length(which(is.finite(Gt.daily[daykeep,j][trunc(day[daykeep])==yrall[yy]])))==0) next
	
		tmp <- tapply(Gt.daily[daykeep,j][trunc(day[daykeep])==yrall[yy]],
		              D.dayInt[daykeep][trunc(day[daykeep])==yrall[yy]],
		              quantile,c(.5,.05,.95))
		tkeep <- as.integer(names(tmp))
		
		tmp <- matrix(unlist(tmp),nrow=3)
	
	if(pl == 1)
		plot(D.daily,Gt.daily[,j],xlab='',ylab='',cex=.1,
	    	 col = "white",axes=FALSE,frame.plot=TRUE)
	
	arrows(D.vec[tkeep]+c(-.09,-.03,.03,.09)[yy]+.25,tmp[2,],
		   D.vec[tkeep]+c(-.09,-.03,.03,.09)[yy]+.25,tmp[3,],
		   col = yy,code=3,angle=90,length=.02)
	
	pl<-0
	
	}
	
	     
	axis(2,cex.axis=.5)
	if(j==6) {
		axis(1,cex.axis=.5)
	    mtext(expression(Vap. ~~ Pres. ~~ Def. ~~ italic(D[t])),side=1,cex=.7,line=2)
		}	
	if(j==1) legend('topright',legend=yrall,text.col=1:4,bty="n",cex=.65)

    par(new = "TRUE",plt = plt.mat[(j-1)*2 +2,],las = 1,cex.axis = 1)

	pl <- 1
	
	for(yy in 1:length(yrall)){
	
		if(length(which(is.finite(Gt.daily[daykeep,j][trunc(day[daykeep])==yrall[yy]])))==0) next
	
		tmp <- tapply(Gt.daily[daykeep,j][trunc(day[daykeep])==yrall[yy]],
		              M.dayInt[daykeep][trunc(day[daykeep])==yrall[yy]],
		              quantile,c(.5,.05,.95))
		tkeep <- as.integer(names(tmp))
		
		tmp <- matrix(unlist(tmp),nrow=3)
	
	if(pl == 1)
		plot(M.daily,Gt.daily[,j],xlab='',ylab='',cex=.1,
	    	 col = "white",axes=FALSE,frame.plot=TRUE)
	
	arrows(M.vec[tkeep]+c(-.06,-.02,.02,.06)[yy]+.05,tmp[2,],
		   M.vec[tkeep]+c(-.06,-.02,.02,.06)[yy]+.05,tmp[3,],
		   col = yy,code=3,angle=90,length=.02)
	
	pl<-0
	
	}
	
	
	if(j==6) 	{
	  axis(1,cex.axis=.5)
	  mtext(expression(Soil ~~ Moist. ~~ italic(M[t])),side=1,cex=.7,line=2)
	  }

	#mtext(expression(italic(f)*(italic(D[t]))*italic(h)*(italic(M[t]))),side=2,outer=TRUE,
	#  line=-1.5,las=3,cex=.7)
	 
	mtext(expression(Can. ~~ Cond. ~~ G[t]),side=2,outer=TRUE,
	  line=-1.5,las=3,cex=.7)
	  
	mtext(spname[j],side=4,las=3,cex=.65)
	
	
	}
	
	dev.off()
	
	