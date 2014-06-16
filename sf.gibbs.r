## Run Gibbs Sampler for Sap Flux Model


for(g in 1:ng){

  if(!SECTION&!EFFSAI){
    tmp <- rescale()
	  specvars   <- tmp$sv2
    }

  ptmp2 <- update_processpars()    #process parameters
  aspec      <- ptmp2$a
  blitespec  <- ptmp2$bl
  bmoistspec <- ptmp2$bm
  gspec      <- ptmp2$gsr
  ul         <- ul + ptmp2$ul
  um         <- um + ptmp2$um
  ug         <- ug + ptmp2$ug

  #priorgt.hi[Q[,1] < 1,] <- .2*gspec[site.all[,'specindex']]  #nighttime conductance can be 20% of maximum
  
  Gsmat <- gssMat(gspec,aspec,blitespec,bmoistspec)  #steady state conductance
  Gsmat[Gsmat>priorgt.hi] <- priorgt.hi[Gsmat>priorgt.hi]
  Gsmat[Gsmat<priorgt.lo] <- priorgt.lo[Gsmat<priorgt.lo]
aig <- update_rint()
va  <- update_va()
 
#  if(CAP) {
#    gtmp  <- update_gt_CAP()
#    Jtmat <- gtmp$Jtmat
#    Wtmat <- gtmp$Wtmat
#	}
#	
#  if(!CAP) gtmp  <- update_gt_NOCAP()
#  
#Gtmat <- gtmp$Gtmat
#
#Gtmat[Gtmat>priorgt.hi] <- priorgt.hi[Gtmat>priorgt.hi]
#Gtmat[Gtmat<priorgt.lo] <- priorgt.lo[Gtmat<priorgt.lo]
#
# 
# 
#  if(!CAP) Jtmat <- Gtmat*e.qt
#  
#  Jpred <- matrix(0,nt,nprobe)
#  pj <- pred_jt(aig,bag,bstat,SECTION)
#    Jpred   <- Jtmat[,site.sensor]*matrix(pj,nt,nprobe,byrow=T)
#    colnames(Jpred) <- probe
#
    
if(CAP){ 
gtmp    <- update_gt_CAP()
Gtmat   <- gtmp$Gtmat
Jtmat   <- gtmp$Jtmat
Wtmat   <- gtmp$Wtmat
Jpred <- matrix(0,nt,nprobe)
pj <- pred_jt(aig,bag,bstat,SECTION)
Jpred   <- Jtmat[,site.sensor]*matrix(pj,nt,nprobe,byrow=T)*lrrat
colnames(Jpred) <- probe
}

if(!CAP){
gtmp    <- update_gt_NOCAP()
Gtmat   <- gtmp$Gtmat
Jtmat <- Gtmat*e.qt
Jpred <- matrix(0,nt,nprobe)
pj <- pred_jt(aig,bag,bstat,SECTION)
Jpred   <- Jtmat[,site.sensor]*matrix(pj,nt,nprobe,byrow=T)*lrrat
colnames(Jpred) <- probe
}
  
  



 if(DPAR){ 
 	pt1   <- update_datapars()
  bag   <- pt1$bag
    ub    <- ub + pt1$ub
  if(BS) {
    bstat <- pt1$bstat
    ubs   <- ubs + pt1$ubs 
	}
  if(!BS) bstat <- rep(1,nspec)}

  
  sigma  <- update_sigma()
  verror <- update_verror()
  #tau    <- update_tau()
  if(CAP&!FIX.CAP){
    tmp  <- update_alpha()
	  alpha <- tmp$alpha
	  uk  <- tmp$uk + uk
    }
  aig <- update_rint()
  va  <- update_va()

#
  dgibbs[g,] <- c(as.vector(bag),bstat)
  ggibbs[g,] <- as.vector(gspec)
  agibbs[g,] <- as.vector(aspec)
  lgibbs[g,] <- as.vector(blitespec)
  mgibbs[g,] <- as.vector(bmoistspec)
  vgibbs[g,] <- c(sigma,verror,va,tau,alpha)
  rgibbs[g,] <- aig
  if(g == gburnin) {
    Ggibbs <- matrix(0,nt,dim(site.all)[1])
      colnames(Ggibbs) <- paste(species[site.all[,'specindex']],
                          SITE[site.all[,'SITE']],sep='.')
    Ggibbs2 <- Ggibbs
    Jgibbs <- matrix(0,nt,dim(site.all)[1])
      colnames(Jgibbs) <- paste(species[site.all[,'specindex']],
                          SITE[site.all[,'SITE']],sep='.')
    Jgibbs2 <- Jgibbs
    Wgibbs <- matrix(0,nt,dim(site.all)[1])
      colnames(Wgibbs) <- paste(species[site.all[,'specindex']],
                          SITE[site.all[,'SITE']],sep='.')
    Wgibbs2 <- Wgibbs
    }
  Ggibbs     <- Ggibbs + Gtmat
  Ggibbs2    <- Ggibbs2 + Gtmat^2

  if(!CAP){
  Jgibbs     <- Jgibbs + apply(Jtmat,1,mean)
  Jgibbs2     <- Jgibbs2 + apply(Jtmat,1,mean)^2
  }
  if(CAP){
  Jgibbs     <- Jgibbs + Jtmat
  Jgibbs2     <- Jgibbs2 + Jtmat^2
  Wgibbs     <- Wgibbs + Wtmat
  Wgibbs2     <- Wgibbs2 + Wtmat^2
  }

  print(g)
  print('data model')
  print(c(bag,bstat))
  print(aig)
  print('process model')
  print(rbind(gspec,aspec))
  print(blitespec)
  print(bmoistspec)
  print('errors and algorithms')
  print(c(sigma,verror,va,tau,alpha))
  print(c("ga  ",paste(ug,sep="  ")),quote=F)
  print(c("light  ",paste(ul,sep="  ")),quote=F)
  print(c("moist  ",paste(um,sep="  ")),quote=F)
  print(c("data  ",paste(ub,sep="  ")),quote=F)
  print(c("alpha  ",paste(uk,sep="  ")),quote=F)
  print(c("pcovga  ",paste(pcovga,sep="  ")),quote=F)
  print(c("pcovQ  ",paste(pcovQ,sep="  ")),quote=F)
  print(c("pcovM  ",paste(pcovM,sep="  ")),quote=F)                   
  
if(g %in% kg){
    if(min(ug) > 12) pcovga <- pcovga*8
    if(min(ug) < 5)  pcovga <- pcovga*.1
    
    if(min(ul) > 12) pcovQ  <- pcovQ*8
    if(min(ul) < 5)  pcovQ  <- pcovQ*.1

    if(min(um) > 12 & pcovM[1,1]<priormat['hiB3',] & pcovM[2,2]<priormat['hiB4',]) pcovM  <- pcovM*8   
    if(min(um) < 5)  pcovM  <- pcovM*.1
											

 if(DPAR){
	  if(min(ub) > 12) pcovba  <- var(dgibbs[(g-40):(g-1),-grep('bstat',colnames(dgibbs))])*Ib
      if(min(ub) < 5)  pcovba  <- pcovba*.2
       ub     <- ub*0
      }
    
    if(CAP & !FIX.CAP){
	  if(min(uk) > 12) pcovK  <- var(vgibbs[(g-40):(g-1),grep('alpha',colnames(vgibbs))])
      if(min(uk) < 5)  pcovK  <- pcovK*.1
      uk     <- uk*0
      }      

    if(BS){
	  if(min(ubs) > 12) pcovbs  <- apply(dgibbs[(g-40):(g-1),grep('bstat',colnames(dgibbs))],2,sd)
      if(min(ubs) < 5)  pcovbs  <- pcovbs*.1
      ubs    <- ubs*0
      }      
    ul     <- ul*0
    um     <- um*0


    ug     <- ug*0

    
#############
#if(!REMOTE) loop.conv()  
#############

  
  }

if (g %in% saveg) save.image(paste(RUN,Treat,year,'-',round(intval[1],2),'.Rdata',sep=''))

}
