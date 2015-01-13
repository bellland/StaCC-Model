
#sap flux model functions

####################################################
tnorm <- function(n,lo,hi,mu,sig){   #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig)

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf] <- hi[z == -Inf]
  z
}


tnorm.mvt <- function(avec,muvec,smat,lo,hi){   
  # truncated multvariate normal
  # avec are current values (e.g., Gibbs sampling)
  # if there are no current values, pass any integer as the first argument
  # muvec is the vector of means
  # smat is the covariance matrix 

  if(length(avec) == 1)avec <- muvec
  avec[avec < lo] <- lo[avec < lo]
  avec[avec > hi] <- hi[avec > hi]
  for(k in 1:length(muvec)){
    piece1 <- smat[-k,k] %*% solve(smat[-k,-k])
    muk <- muvec[k] + piece1 %*% (avec[-k] - muvec[-k])
    sgk <- smat[k,k] - piece1 %*% smat[k,-k]
    if(sgk < .000000001)sgk <- .000000001
    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
  }
  avec
}
####################################################

#Estimate heartwood diameter; s = species; d = diameter (cm)

hd   <- function(s,d) {
 
  ss<- c( 'acru',  'cato',  'list',  'litu',  'pita',  'qual',  'qufa',  'quph')
  a <- c(-3.6293, -5.1282, -2.1939, -0.6551, -9.2071, -3.0033, -3.0033, -3.0033)
  b <- c( 1.7974,  0.8322,  0.4795,  0.5143,  0.7450,  0.9284,  0.9284,  0.9284)
  y <- a[which(ss == s)] + b[which(ss == s)] * d
  if(s == 'acru') y <- exp(a[which(ss == s)] + b[which(ss == s)] * log(d)) 
  y[y<0] <- 0
  y
  }
  
###########################################  

#get effective sapwood area

rescale <- function() {
  sv2 <- specvars
  for(j in 1:nspec){
		sa   <- SA.est(bag[1,j],bag[2,j],rad[specindex==j,'rin'],rad[specindex==j,'rout'])
		for(i in length(SITE)) sv2[species[j],paste('SAP',SITE[i],sep='-')]   <- sum(sa)/10000
  	}
       
  for(j in 1:nspec)
    for(i in 1:nsite)
	  SAPspecies[,site.all[,'SITE']==i & site.all[,'specindex']==j] <- 
	    SAPspecies[,site.all[,'SITE']==i & site.all[,'specindex']==j] * 
		sv2[species[j],paste('SAP',SITE[i],sep='-')] / SPECVARS[species[j],paste('SAP',SITE[i],sep='-')]
  
  SAPtree    <- SAPspecies[,site.sensor]

  list(sv2 = sv2, SAPspecies = SAPspecies, SAPtree = SAPtree)
}

###########################################

#canopy conductance

glagMat <- function(glast,gs,dt,tau) glast + (gs - glast)*(1 - exp(-dt%*%t(tau^-1)))


####################################################

#steady state conductance; gspec = reference conductance, aspec = vpd sensitivity;
#blitespec = light coefficients, bmoistspec = moisture coefficients

gssMat <- function(gspec,aspec,blitespec,bmoistspec){

  out <- func_gQmat(blitespec,Q)*func_gMmat(bmoistspec,M)*func_gDmat(gspec,aspec,D)
    colnames(out) <- paste(species[site.all[,'specindex']],
                         SITE[site.all[,'SITE']],sep='.')
  out  
  }
####################################################

#moisture effects; M = soil moisture; bmoistspec = moisture coefficients

func_gMmat <- function(bmoistspec,M){  #note: using deficit, not M
	
    bm1 <- matrix(bmoistspec[1,site.all[,'specindex']],nt,dim(site.all)[1],byrow=T)
    bm2 <- matrix(bmoistspec[2,site.all[,'specindex']],nt,dim(site.all)[1],byrow=T)
    Mm  <- matrix(M[,site.all[,'SITE']],nt,dim(site.all)[1])
    
      f <- exp(-.5*((Mm - bm2)/bm1)^2)
      f[Mm > bm2] <- 1
	  
    f
}
####################################################

#Light effects; Q = light; blitespec = light parameters
func_gQmat <- function(blitespec,Q){

   bl1 <- matrix(blitespec[1,site.all[,'specindex']],nt,dim(site.all)[1],byrow=T)
   bl2 <- matrix(blitespec[2,site.all[,'specindex']],nt,dim(site.all)[1],byrow=T)
   Qm  <- matrix(Q[,site.all[,'SITE']],nt,dim(site.all)[1])
  
   1 - bl1 * exp(-Qm/(bl2)) 
    
  }
  
####################################################

#VPD effects; D = vpd, gspec = reference conductance, aspec = vpd sensitivity
func_gDmat <- function(gspec,aspec,D){
   
   Dm <- matrix(D[,site.all[,'SITE']],nt,dim(site.all)[1])
   Dm[Dm<exp(-2)] <- exp(-2)
   if(!Dreg){
     out <- matrix(gspec[site.all[,'specindex']],nt,dim(site.all)[1],byrow=T) *             
           (1- matrix(aspec[site.all[,'specindex']],nt,dim(site.all)[1],byrow=T)*log(Dm))  
     } else {
	 out <- matrix(gspec[site.all[,'specindex']],nt,dim(site.all)[1],byrow=T) -               
            matrix(aspec[site.all[,'specindex']],nt,dim(site.all)[1],byrow=T)*log(Dm)       
     }
  out
}

####################################################

#scaling factor for calculating transpiration from conductance
#la = leaf area; sa = sapwood area; dvec = vpd; tvec = temperature

qt <- function(la,sa,dvec,tvec){ #to be multiplied by Gt

 nd <- length(dvec)
 sa.rat <- specvars[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]/
           SPECVARS[,paste('SAP',SITE,sep="-")][site.all[,'specindex']]
 lr <- la/sa/(sa.rat)
 e <- dvec*(tvec + 273)/273/44.6  #mmol m-2 s-1 kPa-1 to mm s-1 
 e <- e/(115.8 + .4236*tvec) #mm s-1 to g m-2 s-1 
 e*lr
}

####################################################

#depth effect times random effect
#deep = relative depth, aig = random effect, ba = depth model parameters
func_ph <- function(deep,aig,ba) {             #mean sensor value
    phi <- exp(-(deep - ba[1,specindex])^2/2/ba[2,specindex]^2)
    phi[deep < ba[1,specindex]] <- 1
  	phi <- phi
  phi * aig
}

####################################################

#(depth times random effect) times canopy status effect
#aig = random effect, ba = depth model parameters, bstat = canopy status effect
pred_jt <- function(aig,ba,bstat){ 
	# sap flux depth, crown status; a vector over trees

  phi <- func_ph(depth,aig,ba)
  bstat <- bstat[specindex]^Status
  phi*bstat
}

###########################################

#capacitance model
#W0 deficit in time t-dt; E = transpiration in time t-dt; DT = time interval, B = proportional lag parameter
WJcalc <- function(W0,E,DT,B) {
	  if(!year %in% seq(1980,2020,by=4)){
        W    <- W0+(E - B*W0/DT)*DT
		W[W<0] <- 0
        J    <- B*W/DT
	    }
  
      #lag years
      if(year %in% seq(1980,2020,by=4)){
        W    <- W0+(E - B*W0/DT)*DT
		W[W<0] <- 0
        J    <- B*W/DT
	    }

		
    list(W = W, J = J)
  }

###########################################
#estimates LAI from sapwood area; s = species, d = diameter (cm), AS = sapwood area

LA.est <- function(s,AS,d) {
  if(s=='acru') dd <- 10^(-0.197 + 0.843 * log10(AS)) 
  if(s=='list') dd <- AS * 0.20
  if(s=='litu') dd <- (AS - 0.157) / 5.028
  if(s=='pita') dd <- AS * 0.20 #from Drake et al. 2010
  if(s=='qual') dd <- (AS - 0.363) / 1.428
  if(s=='qust') dd <- (AS - 0.363) / 1.428
  if(s=='quru') dd <- AS * 0.13 #from Oren et al 
  dd                                                           
  }

###########################################
#estimates effective outer sapwood area of a tree on sapwood profiles and dbh
#a1/a2 = depth parameters, rin = heartwood diameter (cm), rout = sapwood diameter (cm)

SA.est <- function(a1,a2,rin,rout) {
  Adep <- a1 + (sqrt(2*pi)*a2)*(pnorm(1,a1,a2) - .5)           #area under profile
  cj   <- qnorm(((Adep-a1)/2)/(sqrt(2*pi)*a2)+.5,a1,a2)-a1/2   #relative centroid of profile
  dd   <- Adep * (rout-rin) * pi * 2 * (rin + cj * (rout-rin)) #effective area of sapwood
  dd                                                           #units = cm^2
  }
  
#estimate diameter inside bark, s = species, d = diameter (cm)
dib   <- function(s,d) {
  if(s == 'acru') d <- .97 * d
  if(s == 'cato') d <- d - (.105 + 3.070)
  if(s == 'list') d <- d - (3.580 * exp(0.021 * d))/10
  if(s == 'litu') d <- d - (9.687 * exp(0.019 * d))/10
  if(s == 'pita') d <- d - (0.04 + 0.09*d)
  if(s %in% c('qual','qumi','qust')) d <- d - (1.99 * d + 1.322)/10
  if(s %in% c('quco','qufa','quph','quru','quve')) d <- .93 * d
  d
  }
  
#########################################################
#Estimate heartwood diameter, s = species, d = diameter (cm)

hd   <- function(s,d) {
 
  ss<- c( 'acru',  'cato',  'list',  'litu',  'pita',  'qual',  'qufa',  'qumi',  'qust',  'quph')
  
  a <- c(-3.6293, -7.9787, -2.0290, -0.2457, -9.2071, -2.9914, -2.9914, -2.9914, -2.9914, -2.9914)
  b <- c( 1.7974,  1.0623,  0.4858,  0.5307,  0.7450,  0.9462,  0.9462,  0.9462,  0.9462,  0.9462)
  y <- a[which(ss == s)] + b[which(ss == s)] * d
  if(s == 'acru') y <- exp(a[which(ss == s)] + b[which(ss == s)] * log(d)) 
  y[y<0] <- 0
  y
  }
  
###########################################

#update data model parameters
update_datapars <- function(){  #update ba for depth model and bstat for canopy

	jmat <- Jtmat[,site.sensor]
    ub <- rep(0,nspec)

	pba    <- tnorm.mvt(as.vector(bag),as.vector(bag),pcovba,    ##
                          as.vector(loba),as.vector(hiba)) 
	pba    <- matrix(pba,ncol=nspec,byrow=F)
	
    mnow <- jmat*matrix(pred_jt(aig,bag,bstat),nt,nprobe,byrow=T)
    mnew <- jmat*matrix(pred_jt(aig,pba,bstat),nt,nprobe,byrow=T)

    pnow <- rep(0,nspec)
    pnew <- pnow
  
    for(j in 1:nspec){

	  pnow[j] <- sum(dnorm(Jdata[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 mnow[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 sqrt(verror),log=T)) + 
                 sum(dnorm(bag[,j],prior.ba[,j],rbind(sb1[j],sb2[j]),log=T)) 
      pnew[j] <- sum(dnorm(Jdata[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 mnew[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 sqrt(verror),log=T)) + 
                 sum(dnorm(pba[,j],prior.ba[,j],rbind(sb1[j],sb2[j]),log=T)) 
    }

  aa <- exp(pnew - pnow)
  z  <- runif(nspec,0,1)
  wa <- which(z < aa,arr.ind=T)
  bag[,wa]  <- pba[,wa]

  ub[wa] <- ub[wa]+1
  
if(BS){
  ubs <- rep(0,nspec)
  pb   <- tnorm(nspec,priormat['lobstat',],priormat['hibstat',],bstat,pcovbs)

  mnow <- jmat*matrix(pred_jt(aig,bag,bstat),nt,nprobe,byrow=T)
  mnew <- jmat*matrix(pred_jt(aig,bag,pb),nt,nprobe,byrow=T)

  pnow <- rep(0,nspec)
  pnew <- pnow

  for (j in 1:nspec){
    pnow[j] <- sum(dnorm(Jdata[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 mnow[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 sqrt(verror),log=T)) + 
               dnorm(bstat[j],priormat['prior.bstat',j],sqrt(priormat['prior.Vbst',j]),log=T) 
    pnew[j] <- sum(dnorm(Jdata[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 mnew[notmiss[notmiss[,2] %in% which(specindex==j),]],
                 sqrt(verror),log=T)) + 
               dnorm(pb[j],priormat['prior.bstat',j],sqrt(priormat['prior.Vbst',j]),log=T)
  }

  aa <- exp(pnew - pnow)
  z  <- runif(nspec,0,1)
  wa <- which(z < aa,arr.ind=T)
  bstat[wa]  <- pb[wa]
  ubs[wa] <- ubs[wa] + 1
  }

  if(BS) out <- list(bag = bag, bstat = bstat, ub = ub, ubs = ubs)
  if(!BS) out <- list(bag = bag, ub = ub)
  out
}

#update process parameters
update_processpars <- function(){ #update gsr, a, light, sm

  ul   <- rep(0,nspec)
  um   <- ul
  ug   <- ul
  Evec <- 1 - exp(-dt%*%t(tau[site.all[,'specindex']]^-1))

  Gt <- Gtmat
  for(j in 1:nsite) 
    Gt[wdD[wdD[,2]==j,1],site.all[,1]==j] <- NA 
  
  vnow   <- c(gspec,aspec)
  vnew   <- vnow

  ks <- c(1:nspec)
  kn <- ks + nspec

###########################
#VPD Parameters - begin
###########################

 if(!Dreg) 
    vnew <- tnorm.mvt(vnow,vnow,pcovga,c(priormat['loG',],priormat['loA',]),        
            c(priormat['hiG',],priormat['hiA',]))                                 

  #gsref old - begin
  if(Dreg){
    if (dim(site.all)[1]>1)gsprop <- tnorm.mvt(gspec,gspec,pcovga[ks,ks],
                          priormat['loG',],priormat['hiG',])
    if(dim(site.all)[1]==1) gsprop <- tnorm(1,priormat['loG',1],priormat['hiG',1],gspec,
                          sqrt(pcovga[ks,ks]))
  
    vnew   <- c(gsprop,aspec)
    #lamda|gsref
    for(k in kn){
      piece1  <- pcovga[k,-k] %*% solve(pcovga[-k,-k]) 
      muk     <- aspec[k-nspec] + piece1 %*% (vnew[-k] - vnow[-k])
      sgk     <- pcovga[k,k] - piece1 %*% pcovga[k,-k]
      vnew[k] <- tnorm(1,minpropa*gsprop[k-nspec],min(gsprop[k-nspec]/log(max(D)),maxpropa*gsprop[k-nspec]),muk,sqrt(sgk))
      }
    }
  #gsref old - end
  
  pm <- matrix(vnew,2,nspec,byrow=T)
 
  Gs    <- gssMat(gspec,aspec,blitespec,bmoistspec)
#so high D Gs props aren't negative
  Gs[Gs<priorgt.lo] <- priorgt.lo[Gs<priorgt.lo]
  Gs[Gs>priorgt.hi] <- priorgt.hi[Gs>priorgt.hi]
  for(j in 1:nsite){ 
    Gs[M[,j]<moistdef,site.all[,1]==j] <- NA
    Gt[M[,j]<moistdef,site.all[,1]==j] <- NA	
	}	
	  
  
	  
  pnow <- dnorm(aspec,priormat['prior.A',],priormat['prior.sdA',],log=T)         
  if(Dreg) pnow <-  dnorm(aspec,priormat['prior.A',]*pm[1,],priormat['prior.sdA',],log=T)  

	pnow  <- pnow + tapply(apply(as.matrix(dnorm(Gt[-1,],glagMat(Gt[-nt,],Gs[-1,],dt,tau[site.all[,'specindex']]),
                         sqrt(sigma*Evec),log=T)),2,sum,na.rm=T),site.all[,'specindex'],sum,na.rm=T) +
             dnorm(gspec,priormat['prior.G',],priormat['prior.sdG',],log=T)
         
  if(Dreg) mod <- gspec  

  Gp    <- gssMat(pm[1,],pm[2,],blitespec,bmoistspec)
  Gp[Gp<priorgt.lo] <- priorgt.lo[Gp<priorgt.lo]
  Gp[Gp>priorgt.hi] <- priorgt.hi[Gp>priorgt.hi]
  for(j in 1:nsite) Gp[M[,j]<moistdef,site.all[,1]==j] <- NA
  
  pnew <- dnorm(pm[2,],priormat['prior.A',],priormat['prior.sdA',],log=T)         
  if(Dreg) pnew <-  dnorm(pm[2,],priormat['prior.A',]*pm[1,],priormat['prior.sdA',],log=T)  

  pnew  <- pnew + tapply(apply(as.matrix(dnorm(Gt[-1,],glagMat(Gt[-nt,],Gp[-1,],dt,tau[site.all[,'specindex']]),
                         sqrt(sigma*Evec),log=T)),2,sum,na.rm=T),site.all[,'specindex'],sum,na.rm=T) +
             dnorm(pm[1,],priormat['prior.G',],priormat['prior.sdG',],log=T)
	
  if(nspec>1){
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    wa <- which(z < aa,arr.ind=T)
    gspec[wa] <- pm[1,wa]
    aspec[wa] <- pm[2,wa]
    ug[wa] <- 1
    }

  if(nspec==1){
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    if(z<aa) {gspec <- pm[1,1]
             aspec <- pm[2,1] 
             ug <- 1}
    }
      
###########################
#VPD Parameters - end
###########################

###########################
#Light Parameters - begin
###########################

  plite <- tnorm.mvt(as.vector((blitespec)),as.vector((blitespec)),pcovQ,    ##
                    as.vector(priormat[c('loB1','loB2'),]),                  ##
                    as.vector(priormat[c('hiB1','hiB2'),]))                  ##
  pl    <- matrix(plite,2,byrow=F)                                           ##

  Gs    <- gssMat(gspec,aspec,blitespec,bmoistspec)
  Gs[Gs<priorgt.lo] <- priorgt.lo[Gs<priorgt.lo]
  Gs[Gs>priorgt.hi] <- priorgt.hi[Gs>priorgt.hi]

    pnow  <- tapply(apply(as.matrix(dnorm(Gtmat[-1,],glagMat(Gtmat[-nt,],Gs[-1,],dt,tau[site.all[,'specindex']]),
                         sqrt(sigma*Evec),log=T)),2,sum,na.rm=T),site.all[,'specindex'],sum,na.rm=T) +
			 apply(dnorm(blitespec,priormat[c('prior.B1','prior.B2'),],
                   priormat[c('prior.sd1','prior.sd2'),],log=T),2,sum)
        
  Gp    <- gssMat(gspec,aspec,pl,bmoistspec)
  Gp[Gp<priorgt.lo] <- priorgt.lo[Gp<priorgt.lo]
  Gp[Gp>priorgt.hi] <- priorgt.hi[Gp>priorgt.hi]

    pnew  <- tapply(apply(as.matrix(dnorm(Gtmat[-1,],glagMat(Gtmat[-nt,],Gp[-1,],dt,tau[site.all[,'specindex']]),
                         sqrt(sigma*Evec),log=T)),2,sum,na.rm=T),site.all[,'specindex'],sum,na.rm=T) +
             apply(dnorm(pl,priormat[c('prior.B1','prior.B2'),],
                   priormat[c('prior.sd1','prior.sd2'),],log=T),2,sum)

  if(nspec>1){
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    wa <- which(z < aa,arr.ind=T)
    blitespec[,wa]  <- pl[,wa]
    ul[wa] <- 1
  }

  if(nspec==1){
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    if(z<aa) {blitespec <- pl
             ul <- 1}
    }

###########################
#Light Parameters - end
###########################

###########################
#Moist Parameters - begin
###########################

  pmoist <- tnorm.mvt(as.vector((bmoistspec)),as.vector((bmoistspec)),pcovM,    ##
                    as.vector(priormat[c('loB3','loB4'),]),                  ##
                    as.vector(priormat[c('hiB3','hiB4'),]))                  ##
  pm    <- matrix(pmoist,2,byrow=F)                                           ##

  Gs    <- gssMat(gspec,aspec,blitespec,bmoistspec)
  Gs[Gs<priorgt.lo] <- priorgt.lo[Gs<priorgt.lo]
  Gs[Gs>priorgt.hi] <- priorgt.hi[Gs>priorgt.hi]

    pnow  <- tapply(apply(as.matrix(dnorm(Gtmat[-1,],glagMat(Gtmat[-nt,],Gs[-1,],dt,tau[site.all[,'specindex']]),
                         sqrt(sigma*Evec),log=T)),2,sum,na.rm=T),site.all[,'specindex'],sum,na.rm=T) +
             apply(dnorm(bmoistspec,priormat[c('prior.B3','prior.B4'),],
                   priormat[c('prior.sd3','prior.sd4'),],log=T),2,sum)
        
  Gp    <- gssMat(gspec,aspec,blitespec,pm)
  Gp[Gp<priorgt.lo] <- priorgt.lo[Gp<priorgt.lo]
  Gp[Gp>priorgt.hi] <- priorgt.hi[Gp>priorgt.hi]

    pnew  <- tapply(apply(as.matrix(dnorm(Gtmat[-1,],glagMat(Gtmat[-nt,],Gp[-1,],dt,tau[site.all[,'specindex']]),
                         sqrt(sigma*Evec),log=T)),2,sum,na.rm=T),site.all[,'specindex'],sum,na.rm=T) +
             apply(dnorm(pm,priormat[c('prior.B3','prior.B4'),],
                   priormat[c('prior.sd3','prior.sd4'),],log=T),2,sum)

  if(nspec>1){
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    wa <- which(z < aa,arr.ind=T)
    bmoistspec[,wa]  <- pm[,wa]
    um[wa] <- 1
  }

  if(nspec==1){
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    if(z<aa) {bmoistspec <- pm
             um <- 1}
    }
###########################
#Moist Parameters - begin
###########################

  list(gsr = gspec,a = aspec, bl = blitespec, bm = bmoistspec, ug = ug, um = um, ul = ul)

}


#update process error

update_sigma <- function(){

  Evec <- matrix( (1 - exp(-dt%*%t(tau[site.all[,'specindex']]^-1))), (nt-1), dim(site.all)[1])
  ssq  <- (Gtmat[-1,] - glagMat(Gtmat[-nt,],Gsmat[-1,],dt,tau[site.all[,'specindex']]) )^2/Evec

  u1 <- s1 + length(Gtmat[-1,])/2
  u2 <- s2 + .5*sum(ssq,na.rm=T)
  1/rgamma(1,u1,u2)
}

#update stomatal lags

update_tau <- function(){

  ptau <- tnorm(nspec,priormat['lotau',1],priormat['hitau',1],tau,.000005)

  Enow  <- matrix( (1 - exp(-dt%*%t(tau[site.all[,'specindex']]^-1))), (nt-1), dim(site.all)[1])
  Enew  <- matrix( (1 - exp(-dt%*%t(ptau[site.all[,'specindex']]^-1))), (nt-1), dim(site.all)[1])

  gtnow <- glagMat(Gtmat[-nt,],Gsmat[-1,],dt,tau[site.all[,'specindex']])
  gtnew <- glagMat(Gtmat[-nt,],Gsmat[-1,],dt,ptau[site.all[,'specindex']])

  if(dim(site.all)[1]>1){
    pnow  <- apply(dnorm(Gtmat[-1,],gtnow,sqrt(Enow*sigma),log=T),2,sum,na.rm=T) 
	pnow  <- tapply(pnow,site.all[,'specindex'],sum,na.rm=T) + dnorm(tau,priormat['prior.tau',],priormat['prior.sdtau',],log=T)
      
    pnew  <- apply(dnorm(Gtmat[-1,],gtnew,sqrt(Enew*sigma),log=T),2,sum,na.rm=T) 
	pnew  <- tapply(pnew,site.all[,'specindex'],sum,na.rm=T) + dnorm(ptau,priormat['prior.tau',],priormat['prior.sdtau',],log=T)
    
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    wa <- which(z < aa,arr.ind=T)
    tau[wa] <- ptau[wa]
    }
  if(nspec==1){
    pnow  <- sum(dnorm(Gtmat[-1,],gtnow,sqrt(Enow*sigma),log=T),na.rm=T) +
             dnorm(tau,priormat['prior.tau',],priormat['prior.sdtau',],log=T)
      
    pnew  <- sum(dnorm(Gtmat[-1,],gtnew,sqrt(Enew*sigma),log=T),na.rm=T) +
             dnorm(ptau,priormat['prior.tau',],priormat['prior.sdtau',],log=T)
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    if(z<aa) tau <- ptau
    }
        
  tau
}

#update capacitance

update_beta <- function(){ ##FIX
  
  uk   <- rep(0,nspec)

  if(!year %in% seq(1980,2020,by=4))
     DT <- dt*365*24*60*60
  
  if(year %in% seq(1980,2020,by=4))
     DT <- dt*366*24*60*60

	 #draw capacitance paramteers
  if(nspec>1) 
    pbeta <- tnorm.mvt(beta,beta,pcovK,priormat['lobeta',],priormat['hibeta',])
  if(nspec==1) 
    pbeta <- tnorm(1,priormat['lobeta',],priormat['hibeta',],beta,pcovK)

  C <- matrix(pred_jt(aig,bag,bstat),nt,nprobe,byrow=T)  

	gmat <- Gtmat
    jnow <- Jtmat
	wnow <- Wtmat
	jnew <- jnow
	wnew <- wnow
	
	for(t in 2:nt){
	  
	  tmp <- WJcalc(wnew[t-1,],Gtmat[t-1,]*e.qt[t-1,],DT[t-1],pbeta)
	    wnew[t,] <- tmp$W
	    jnew[t,] <- tmp$J
	  }
	jnow<- jnow[,site.sensor]
	jnew<- jnew[,site.sensor]
	
	jnow[is.na(Jdata)] <- NA
	jnew[is.na(Jdata)] <- NA
	  
    #v <- 
	#V <- 
	
	pnow  <- tapply(apply(dnorm(Jdata[-1,],jnow[-1,]*C[-1,],sqrt(verror),log=T),2,sum,na.rm=T),specindex,sum)
    pnew  <- tapply(apply(dnorm(Jdata[-1,],jnew[-1,]*C[-1,],sqrt(verror),log=T),2,sum,na.rm=T),specindex,sum)
	
	pnow  <- pnow + dnorm(beta,priormat['prior.beta',],priormat['prior.beta',],log=T)
    pnew  <- pnew + dnorm(pbeta,priormat['prior.beta',],priormat['prior.beta',],log=T)
    
    aa <- exp(pnew - pnow)
    z  <- runif(nspec,0,1)
    wa <- which(z < aa,arr.ind=T)
    beta[wa] <- pbeta[wa]
    uk[wa]  <- 1    
	
  list(beta = beta, uk = uk)
}

#update data error
update_verror <- function(){

  jmat <- Jtmat[,site.sensor]
  jmat[is.na(Jdata)] <- NA

  mnow <- jmat*matrix(pred_jt(aig,bag,bstat),nt,nprobe,byrow=T)  
  
  u1 <- v1 + length(notmiss)/2
  u2 <- v2 + .5*sum( (Jdata[notmiss] - mnow[notmiss])^2 ,na.rm=T)
  1/rgamma(1,u1,u2)

}

#update GT for capacitance model
update_gt_CAP <- function(){


  Ev <- 1 - exp(-dt%*%t(tau[site.all[,'specindex']]^-1))   #stomatal lag
  smat <- matrix(0,dim(site.all)[1],nsite)

  deficit[1,] <- deficit_0

  C <- matrix(pred_jt(aig,ba,bstat),nt,nprobe,byrow=T)  
  jkeep <- (Jdata*0+1)
        
  
  if(!year %in% seq(1980,2020,by=4))
     DT <- dt*365*24*60*60
  
  if(year %in% seq(1980,2020,by=4))
     DT <- dt*366*24*60*60

  
  for(t in 1:nt){

    smat <- smat*0
    V    <- rep(0,dim(site.all)[1])
    v    <- V
         

	if(t %in% starts){	 
	
	  tmp <- WJcalc(Wtmat[t,]*0,0,DT[t],beta)  #FIX
        Wtmat[t,] <- tmp$W
        Jtmat[t,] <- tmp$J
	  }
	  
    if(!t %in% starts){                                   #transition from t-1 to t
    #sap flux
	  
      tmp <- WJcalc(Wtmat[t-1,],Gtmat[t-1,]*e.qt[t-1,],DT[t-1],beta)
        Wtmat[t,] <- tmp$W
        Jtmat[t,] <- tmp$J
	
      B <- Gtmat[t-1,] + (Gsmat[t,] - Gtmat[t-1,])*Ev[t-1,]
      v <- v + B/sigma/Ev[t-1,]
      V <- V + 1/sigma/Ev[t-1,]
    }
    if(!t %in% stops){                                  #transition from t to t+1
	  JKij <- (Jdata[t+1,] - beta[site.sensor] * C[t,] * Wtmat[t,site.sensor] *( 1 - beta[site.sensor])/DT[t]) * 
	          (beta[site.sensor] * C[t,] * e.qt[t,site.sensor])
	  
	  JK2  <- (beta[site.sensor] * C[t,] * e.qt[t,site.sensor])^2
	
      Vo  <- JK2 * jkeep[t,] /verror
      vo  <- JKij * jkeep[t,] /verror
	  
      V <- V + tapply(Vo,site.sensor,sum,na.rm=T)
      v <- v + tapply(vo,site.sensor,sum,na.rm=T)

      A  <- (Gtmat[t+1,] - Gsmat[t+1,]*Ev[t,])*(1 - Ev[t,])
      V  <- V + (1 - Ev[t,])^2/sigma/Ev[t,]
      v  <- v + A/sigma/Ev[t,]
    }
    
	#draw Gtmat[t,]
    
	Gtmat[t,] <- tnorm(dim(site.all)[1],priorgt.lo[t,],priorgt.hi[t,],v/V,sqrt(1/V))

  }

  list(Gtmat = Gtmat, Jtmat = Jtmat, Wtmat = Wtmat)

}

#update random intercepts
update_rint <- function(){  

  jmat <- matrix(Jtmat[,site.sensor],nt)
  jmat[is.na(Jdata)] <- NA
  rint <- meana[specindex]

    j2ik <- tapply(jmat[notmiss]^2,notmiss[,2],sum,na.rm=T)
    jjik <- tapply(jmat[notmiss]*Jdata[notmiss],notmiss[,2],sum,na.rm=T)

  
    evec <- exp(-(depth - bag[1,specindex])^2/2/(bag[2,specindex]^2))
    evec[depth<bag[1,specindex]] <- 1
    
    
    v <- evec*(bstat[specindex]^Status)*jjik/verror + meana[specindex]/va
    V <- 1/((evec^2)*(bstat[specindex]^(2*Status))*j2ik/verror + 1/va)

    rint <- tnorm(nprobe,lor,hir,v*V,sqrt(V))
  
  ## rescale to intercept
  rint <- rint * meana[specindex] / tapply(rint,specindex,mean)[specindex]  
  #rint
  #temporary

 #  ai <- 1/bag[1]*(1 - bag[2]*sqrt(pi/2))
 # rep(ai,nprobe)
 rint

}

#variance for random effects
update_va <- function(){    

  z1 <- va1 + nprobe/2
  z2 <- va2 + .5*sum( (aig-meana[specindex])^2) 
  1/rgamma(1,z1,z2)
}

###########################################


getTime <- function(input) {
  ly <- seq(1980,2020,by=4)
  xxx<-which(input[,'year'] %in% ly)

  if(max(input[,'Time'])>=30)
    out <- input[,'year'] + (input[,'DOY']-1)/365 + trunc(input[,'Time']/100)/24/365 + 
	       (input[,'Time']-trunc(input[,'Time']/100)*100)/60/24/365
  if(max(input[,'Time'])<30)
    out <- input[,'year'] + (input[,'DOY']-1)/365 + input[,'Time']/24/365
  
	
  if(length(xxx)>0){
    if(max(input[xxx,'Time'])>=30) 
      out[xxx] <- input[xxx,'year'] + (input[xxx,'DOY']-1)/366 + trunc(input[xxx,'Time']/100)/24/366 + 
	              (input[xxx,'Time']-trunc(input[xxx,'Time']/100)*100)/60/24/366
    if(max(input[xxx,'Time'])<30) 
      out[xxx] <- input[xxx,'year'] + (input[xxx,'DOY']-1)/366 + input[xxx,'Time']/24/366
    }
  out
  }

getDay <- function(input) {
ly <- seq(1980,2020,by=4)
  outx<-input[,'year'] + (input[,'DOY']-1)/365
  xxx<-which(input[,'year'] %in% ly)
  outx[xxx]<-input[xxx,'year'] + (input[xxx,'DOY']-1)/366
  outx
  }

###############################
snip <- function(tm,st1, st2){
  ly <- seq(1980,2020,by=4)

  tm <- tm-trunc(tm)
  hhmm <- (0:2)/24
  
  eq <- which(findInterval(tm*365 - trunc(tm*365),hhmm)==1) #times at equilibrium
  if(trunc(timeall[1]) %in% ly)
    eq <- which(findInterval(tm*366 - trunc(tm*366),hhmm)==1) #times at equilibrium
  
  snew1 <- st1
  snew2 <- st2
  for(sn in 1:length(st1)){
    snew1[sn] <- eq[min(which(eq>st1[sn]))]
    snew2[sn] <- eq[max(which(eq<st2[sn]))]
    }
	
  list(snew1 = snew1, snew2 = snew2)    
  }  
  
###########################################  
#LOAD DATA
inData <- function(){

  for(j in 1:length(year)){
    if(j == 1){
	  b2       <- as.matrix(read.csv(SENS.FILE[j],header=T))  #import sensor metadata
      datamat  <- as.matrix(read.csv(FLUX.FILE[j],header=T))  #import sensor flux data

      leafarea <- as.matrix(read.csv(LAI.FILE[j],header=T))   #import leaf area data
      saparea  <- as.matrix(read.csv(SAI.FILE[j],header=T))   #import sapwood area data

      Q.dat<-as.matrix(read.csv(PAR.FILE[j],header=T))        #import light data
      M.dat<-as.matrix(read.csv(SM.FILE[j],header=T))         #import soil moisture data
      D.dat<-as.matrix(read.csv(VPD.FILE[j],header=T))        #import vpd data
      T.dat<-as.matrix(read.csv(TEMP.FILE[j],header=T))       #import temperature data
      }
	if(j>1){
	  datamat  <- rbind(datamat,as.matrix(read.csv(FLUX.FILE[j],header=T)))
	  leafarea <- rbind(leafarea,as.matrix(read.csv(LAI.FILE[j],header=T)))
	  saparea  <- rbind(saparea,as.matrix(read.csv(SAI.FILE[j],header=T)))
	  
	  Q.dat <- rbind(Q.dat,as.matrix(read.csv(PAR.FILE[j],header=T)))
	  M.dat <- rbind(M.dat,as.matrix(read.csv(SM.FILE[j],header=T)))
	  D.dat <- rbind(D.dat,as.matrix(read.csv(VPD.FILE[j],header=T)))
	  T.dat <- rbind(T.dat,as.matrix(read.csv(TEMP.FILE[j],header=T)))
	  }
    }
  leafarea<-leafarea[is.na(leafarea[,1])==0,]
  
  LAIspecies <- leafarea         #leaf area (m^2 m^-2)
  SAPspecies <- saparea/10000    #sapwood area (m^2 m^-2)
  
  
  #indices for covariates
  Q.col   <- colnames(Q.dat)[-(1:4)]
  D.col   <- colnames(D.dat)[-(1:4)]
  M.col   <- colnames(M.dat)[-(1:4)]
  T.col   <- colnames(T.dat)[-(1:4)]

  #
  if(!is.na(SP)) {
    
	b2      <- b2[b2[,'Species']==SP,]
	datamat <- datamat[,c(1:4,grep(SP,colnames(datamat)))]
	
    }
    probe   <- b2[,'colname']      #probe name
    pindex  <- which(colnames(datamat) %in%  probe,arr.ind=T) #columns with flux data
	
  #species names and indices
  species <- as.character(sort(unique(b2[,'Species'])))
  nspec   <- length(species)
  specindex <- match(b2[,'Species'],species)

  #plot names and indices
  plots   <- as.character(sort(unique(b2[,'Ring'])))
  nplot   <- length(plots)
  plotindex <- match(b2[,'Ring'],plots)
  
  #ratio of sapwood at different depths
  SAratio <- numeric(0)

  SITE  <- unique(b2[,'Site'])
  nsite <- length(SITE)
  
  specvars <- matrix(NA,nspec,nsite*2)
  cs <- c('n','SAP')
  colnames(specvars) <- outer(cs,SITE,paste,sep='-')
  rownames(specvars) <- species

  spectab <- table(b2[,'Species'],b2[,'Site'])
  for (i in 1:nsite) 
    specvars[,paste('n-',SITE[i],sep="")] <- spectab[,SITE[i]]
	
  

  Q     <- numeric(0)
  D     <- numeric(0)
  M     <- numeric(0)
  Temp  <- numeric(0)
  Jdata <- numeric(0)

 datacols <- c('year','JDT','DOY','Time')

 dgaps   <- numeric(0)

 
   #get data range
   checkrow <- apply(is.finite(datamat[,pindex]),1,sum,na.rm=T)
   cc      <- min(which(cumsum(checkrow) > 0),na.rm=T)
   cm      <- max(which(checkrow > 0),na.rm=T)

   #get gaps
   crow  <- which(checkrow > 0)
   gaps  <- diff(crow)
   dd    <- which(gaps > gapfill)
   dgaps <- c(dgaps,dd)

   skipit <- numeric(0)
   if(cc>1) skipit <- c(1:(cc-1))

   if(length(dd)>0)
   for(jj in 1:length(dd)){
     skipit <- c(skipit,c((1+crow[dd[jj]]):(crow[dd[jj]] + gaps[dd[jj]])))
   }
   if(cm<nrow(datamat)) skipit <- c(skipit,c((cm+1):nrow(datamat)))
   
   #create keep vectors
   checkrow <- c(1:nrow(datamat))
   if(length(skipit)>0) checkrow <- checkrow[-skipit]

   Q.keep    <- match(getTime(datamat[checkrow,]),getTime(Q.dat))
   D.keep    <- match(getTime(datamat[checkrow,]),getTime(D.dat))
   M.keep    <- match(getTime(datamat[checkrow,]),getTime(M.dat))
   Temp.keep <- match(getTime(datamat[checkrow,]),getTime(T.dat))
   LAI.keep  <- match(getDay(datamat[checkrow,]),getDay(leafarea))
   SAI.keep  <- match(getTime(datamat[checkrow,]),getTime(saparea))
  
   SITE <- matrix(unlist(strsplit(Q.col,"\\.")),ncol=length(Q.col))[2,]
#   if(!MULTCAN) SITE <- can.choose

  if(DATA=='face'){
  	siteindex <- match(b2[,'Ring'],SITE)
  	}else{
  	siteindex <- match(b2[,'Site'],SITE)}
   
   if(DATA=='face') tmp <- match(paste(species[specindex],".R",SITE[siteindex],".",Fert,sep=""),colnames(leafarea))
   if(DATA!='face') tmp <- match(paste(species[specindex],SITE[siteindex],sep="."),colnames(leafarea))

   
     Q      <- matrix(Q.dat[Q.keep,Q.col[matrix(unlist(strsplit(Q.col,"\\.")),2)[2,]==SITE]],ncol=1)           #automate ncol
     D      <- matrix(D.dat[D.keep,D.col[matrix(unlist(strsplit(D.col,"\\.")),2)[2,]==SITE]],ncol=1)           #automate ncol
     M      <- matrix(M.dat[M.keep,M.col[matrix(unlist(strsplit(M.col,"\\.")),2)[2,]==SITE]],ncol=1)           #automate ncol
     if(max(M,na.rm=T) > 1) M <- M/100                     #if in percent change to fraction
     Temp   <- matrix(T.dat[Temp.keep,T.col[matrix(unlist(strsplit(T.col,"\\.")),2)[2,]==SITE]],ncol=1)        #automate ncol

  #exclude time points where no data are available at 
  #  beginning and end of study period 
  Jdata  <- datamat[checkrow,pindex]
  SAPtree <- saparea[SAI.keep,tmp]/10000
  LAItree <- leafarea[LAI.keep,tmp]
  SAPspecies <- saparea[SAI.keep,sort(unique(tmp))]/10000
  LAIspecies <- leafarea[LAI.keep,sort(unique(tmp))]

    SAPspecies <- matrix(0,nrow(LAItree),ncol=nspec)
    LAIspecies <- matrix(0,nrow(LAItree),ncol=nspec)
    for(s in 1:nspec){
	  if(DATA=='face') tmp <- match(paste(species[specindex[specindex==s]],".R",SITE[siteindex[specindex==s]],".",Fert,sep=""),colnames(leafarea))
      if(DATA!='face') tmp <- match(paste(species[specindex[specindex==s]],SITE[siteindex[specindex==s]],sep="."),colnames(leafarea))

	  SAPspecies[,s] <- apply(saparea[SAI.keep,tmp]/10000,1,mean,na.rm=T)
      LAIspecies[,s] <- apply(leafarea[LAI.keep,tmp],1,mean,na.rm=T)	
	  }
	
  datamat <- datamat[checkrow,]
  
  #year  <- unique(datamat[,'year'])
  yrvec <- datamat[,'year']
  timeall <- getTime(datamat)
  nt   <- length(timeall)
  dt   <- diff(timeall)
  Jdata[!is.finite(Jdata)] <- NA

  #vectors of starts and stops to remove data gaps
  gaps <- which(diff(timeall) > diff(timeall)[1]*gapfill) + 1  #where new seq begins
  starts <- c(1,gaps)
  stops  <- c((gaps-1),nt)  
  #dt[stops[-length(stops)]] <- dt[starts[-1]]  #what does this do?

  #plot covariates
  if(length(M)==nt){
    par(mfrow=c(5,1),mar=c(2,4,2,2))
    plot(timeall,Temp,type='l')
    plot(timeall,D,type='l')
    plot(timeall,Q,type='l')
    plot(timeall,M,type='l',ylim=range(M,na.rm=T))
  }
  
  if(length(M)>nt){
    par(mfrow=c(5,1),mar=c(2,4,2,2))
    plot(timeall,Temp[,1],type='l')
    plot(timeall,D[,1],type='l')
    plot(timeall,Q[,1],type='l')
    plot(timeall,M[,1],type='l',ylim=range(M,na.rm=T))
  }
  
#  ntime <- length(keeptime)
 
  ntree <- nrow(b2) #number of sensors

  Status <- as.integer(b2[,'Status'])  #0 - suppressed, 1 - codominant

  #Jdata[wdD,] <- NA

  nprobe  <- ncol(Jdata) #number of sensors, why repeat this???

  #plot sap flux data
  plot(timeall,Jdata[,1],col=1,ylim=c(0,100),type='l')
  for(j in 1:nprobe){
   lines(timeall,Jdata[,j],col=j)
  }

  # draw depth values
    rad <- matrix(NA,dim(b2)[1],2)
    
	for(j in 1:length(species))
	  rad[which(b2[,'Species']==species[j]),2] <- 
	    dib(species[j],as.numeric(b2[b2[,'Species']==species[j],paste('DBH',min(year),sep="")]))/2 	#is it okay to use the minimum year?
  
    for(j in 1:dim(rad)[1])
      rad[j,1] <- hd(b2[j,'Species'],rad[j,2]*2)/2
      colnames(rad) <- c('rin','rout')

	for(j in 1:nspec)
	  for(i in 1:nsite)
	    specvars[species[j],paste('SAP',SITE[i],sep="-")] <- 
		  sum(pi*rad[specindex==j & siteindex==i,2]^2 - pi*rad[specindex==j & siteindex==i,1]^2)/10000
	  
	
	depth <- as.numeric(b2[,'Depth'])/apply(rad,1,diff)
	
	
  
  #if(RELM & MULTCAN) M <- (M - matrix(apply(M,2,min),nt,ncol(M),byrow=T))/
  #              (matrix(apply(M,2,max),nt,ncol(M),byrow=T) - matrix(apply(M,2,min),nt,ncol(M),byrow=T))

  #if(RELM & !MULTCAN) M <- (M - min(M))/(max(M) - min(M))
             	
  list(Q = Q, M = M, D = D, Temp = Temp, ntree = ntree, nprobe = nprobe, rad = rad, 
       depth = depth, Status = Status, Jdata = Jdata, starts = starts, 
       stops = stops, gaps = gaps, timeall = timeall, nt = nt, plots = plots, nplot = nplot, 
       LAIspecies = LAIspecies, LAItree = LAItree, specindex = specindex, plotindex = plotindex, 
       siteindex = siteindex, SAPspecies = SAPspecies, SAPtree = SAPtree, probe = probe, SITE = SITE, 
       species = species, nspec = nspec, specvars = specvars, SAratio = SAratio)

}

spltname <- function(x,on) return(matrix(unlist(strsplit(x,on)),length(x),byrow=T))


getSection <- function(){

  wt      <- which(trunc(timeall)%in%year & (timeall-trunc(timeall)) >= intval[1] & (timeall-trunc(timeall)) <= intval[2])
  nt      <- length(wt)

  starts <- starts[starts %in% wt]
  stops <- stops[stops %in% wt]

  st1     <- unique(wt[findInterval(starts,wt)])
  st2     <- unique(wt[findInterval(stops,wt)])
 
  if(!wt[1] %in% st1)  st1 <- c(wt[1],st1)
  if(!wt[nt] %in% st2) st2 <- c(st2,wt[nt])

  ##begin -- start and stop timeseries at 12am
  tmp <- snip(timeall,st1,st2)
    st1 <- tmp$snew1
    st2 <- tmp$snew2
	
  stkeep <- which(st2-st1>0)
    st1 <- st1[stkeep]
    st2 <- st2[stkeep]
  
  tmp <- as.integer()
  for(j in 1:length(st1))
    tmp <- c(tmp,st1[j]:st2[j])
  wt <- tmp
  
	
	
  timenew <- timeall[wt]

  st1     <- st1 - wt[1] + 1
  st2     <- st2 - wt[1] + 1

  timeall <- timenew
  nt      <- length(timeall)

  abline(v=c(timeall[1],timeall[nt]),lwd=5,col='red')

  gaps <- which(diff(timeall) > diff(timeall)[1]*gapfill) + 1  #where new seq begins
  starts <- c(1,gaps)
  stops  <- c((gaps-1),nt)  

  Jdata <- Jdata[wt,]

  np <- nsite

    Temp  <- matrix(Temp[wt,],ncol=np)
    Q     <- matrix(Q[wt,],ncol=np)
    M     <- matrix(M[wt,],ncol=np)
    D     <- matrix(D[wt,],ncol=np)

	
  LAItree    <- matrix(LAItree[wt,],length(wt))
  LAIspecies <- matrix(LAIspecies[wt,],length(wt))
  SAPtree    <- matrix(SAPtree[wt,],length(wt))
  SAPspecies <- matrix(SAPspecies[wt,],length(wt))
  
  deficit <- matrix(0,nt,np)
  dt    <- diff(timeall)
  dt[stops[-length(stops)]] <- dt[starts[-1]]
  dt[stops[-length(stops)]] <- dt[starts[-1]]
  
  dt[dt > dt[1]*2] <- dt[1]*2  #remove large dt's
  
  #gap fill soil moisture --temporary
    for(ss in 1:np) {
      y  <- diff(log(M[,ss]))
      wy <- which(y < .02,arr.ind=T)
      y  <- y[wy]
      x <- cbind(rep(1,nt),M[,ss],Temp[,ss])[-nt,]
      x  <- x[wy,]
      b  <- solve(crossprod(x))%*%crossprod(x,y)
      for(t in 1:(nt-1)){
  	    if(!is.finite(M[t+1,ss])){
  	 	  M[t+1,ss] <- M[t,ss]*exp(as.vector(c(1,M[t,ss],Temp[t,ss]))%*%b)
  	}}}
	
if(np==1){  
  par(mfrow=c(5,1),mar=c(2,4,2,2))
  plot(timeall,Temp,type='l')
  plot(timeall,D,type='l')
  plot(timeall,Q,type='l')
  plot(timeall,M,type='l')
  plot(timeall,Jdata[,1],ylim=c(0,100),type='l')
  for(j in 1:nprobe){
   lines(timeall,Jdata[,j],col=j)
   }

  for(j in 1:length(starts)){
   xx <- timeall[starts[j]:stops[j]]
   lines(xx,rep(100,length(xx)),lwd=5)
   }
  }
if(np>1) {
  par(mfrow=c(5,1),mar=c(2,4,2,2))
  plot(timeall,Temp[,1],type='l')
  plot(timeall,D[,1],type='l')
  plot(timeall,Q[,1],type='l')
  plot(timeall,M[,1],type='l')
  plot(timeall,Jdata[,1],ylim=c(0,100))
  for(j in 1:nprobe){
   lines(timeall,Jdata[,j],col=j)
   }

  for(j in 1:length(starts)){
   xx <- timeall[starts[j]:stops[j]]
   lines(xx,rep(100,length(xx)),lwd=5)
   }
  }  

  list(starts = starts, stops = stops, Jdata = Jdata, Temp = Temp, 
       Q = Q, M = M, D = D, deficit = deficit, timeall = timeall, dt = dt,
	   LAItree = LAItree, LAIspecies = LAIspecies, SAPtree = SAPtree, 
	   SAPspecies = SAPspecies)
}


########################
## Plotting Functions ##
########################

## plot midday Js by parameter
midJs <- function(Js,covar,label,i) {
  TT <- which((365*(timeall-2007)-trunc(365*(timeall-2007)))>.45 & 
              (365*(timeall-2007)-trunc(365*(timeall-2007)))<.55)
  Jscol <- which(specindex==i)
  covar <- covar[TT]
  for (m in 1:length(Jscol)) {
    if(m==1) plot(covar,Js[TT,Jscol[m]],xlab=label,ylab='Js',main=species[i],
                  xlim=range(covar,na.rm=T),ylim=range(Js[TT,],na.rm=T))
    if(m>1) points(covar,Js[TT,Jscol[m]],col=m)
  }
}

plot.lags <-function() {
 meanlag  <- rep(0,nprobe+nspec)
hours    <- 24*(365*(timeall-trunc(timeall))-trunc(365*(timeall-trunc(timeall))))   
hours <- which(hours>=0 & hours<=24) 

par(mfrow=c(3,4))
  QT <- qt(LAIspecies[j], mean(SAPtree[specindex==j]), D[hours], Temp[hours])
  for (j in 1:nprobe) {
    ccx <- ccf(D[hours],Jdata[hours,j],type='correlation',lag.max=6,
    na.action=na.contiguous,plot=F)
    print(paste(probe[j],species[specindex[j]],Status[j],sep='-'))
    print(ccx)
    plot(ccx,type='l',main=paste(probe[j],species[specindex[j]],Status[j],sep='-'))
    meanlag[j] <- ccx[[4]][which(ccx[[1]]==max(ccx[[1]]))]
    abline(v=meanlag[j])
  }
  
  ccx <- ccf(D[hours],apply(Jdata[hours,which(specindex==1)],1,mean,na.rm=T),type='correlation',lag.max=6,
  na.action=na.contiguous,plot=F)
  print(species[1])
  print(ccx)
  plot(ccx,type='l',main='average')
  meanlag[nprobe+1] <- ccx[[4]][which(ccx[[1]]==max(ccx[[1]]))]
  abline(v=meanlag[nprobe+1])
  for(j in 2:nspec) {
    ccx <- ccf(D[hours],apply(Jdata[hours,which(specindex==j)],1,mean,na.rm=T),type='correlation',lag.max=6,
    na.action=na.exclude,plot=F)
    print(species[j])
    print(ccx)
    lines(ccx[['lag']],ccx[['acf']],col=j)
    meanlag[nprobe+j] <- ccx[[4]][which(ccx[[1]]==max(ccx[[1]]))]
    abline(v=meanlag[nprobe+j],col=j)  
  }
  legend('topleft',legend=species,text.col=1:4)
 

pdf('lags.pdf')
par(mfrow=c(1,1))
for(j in 1:nprobe){
  plot(timeall,Jdata[,j],ylim=c(0,100),type='l',lwd=3,xlim=c(2007.68,2007.69))
#  lines(timeall,Jpred[,j],col='red')
  title(paste(species[specindex[j]],Status[j],probe[j],sep='-'))
  lines(timeall[-((nt-2):nt)],Jdata[-(1:3),j],col='blue',lwd=2,lty=2)
  lines(timeall,D*10,col='blue',lwd=2)
  lines(timeall,Q/50,col='orange',lwd=2)
#  lines(timeall,gmean[,specindex[j]],col='purple')
  legend('topleft',legend=c('Jdata','D*10','Q/50'),text.col=c('black','blue','orange'))
  abline(v=seq(2007,2008,by=1/365),lty=2)
}

#dev.print(device=postscript,file='conductanceBySpec.ps',width=7, height=10, horizontal=FALSE)
dev.off()

}

loop.conv <- function() {
    par(mfrow=c(5,max(2,nspec)),mai=c(.2,.4,.4,.1))
    for(jj in 1:nspec) {
      plot(max(1,g-1000):g,ggibbs[max(1,g-1000):g,jj],type='l',xlab='gibbs step',
           ylab=colnames(ggibbs)[jj],main=colnames(ggibbs)[jj])
	  lines((g-39):g,ggibbs[(g-39):g,jj],col='red')
      }
    for(jj in 1:nspec) {
      plot(max(1,g-1000):g,agibbs[max(1,g-1000):g,jj],type='l',xlab='gibbs step',
           ylab=colnames(agibbs)[jj],main=colnames(agibbs)[jj])
	  lines((g-39):g,agibbs[(g-39):g,jj],col='red')
	  }
    for(jj in 1:nspec) {
      plot(lgibbs[max(1,g-1000):g,jj*2-1],lgibbs[max(1,g-1000):g,jj*2],type='l',xlab='bl1',
           ylab='bl2',main=species[jj])
	  lines(lgibbs[(g-39):g,jj*2-1],lgibbs[(g-39):g,jj*2],col='red')
      }
    for(jj in 1:nspec) {
      plot(mgibbs[max(1,g-1000):g,jj*2-1],mgibbs[max(1,g-1000):g,jj*2],type='l',xlab='bm1',
           ylab='bm2',main=species[jj])
	  lines(mgibbs[(g-39):g,jj*2-1],mgibbs[(g-39):g,jj*2],col='red')
      }
   plot(max(1,g-1000):g,vgibbs[max(1,g-1000):g,1],type='l',xlab='gibbs step',
           ylab=colnames(vgibbs)[1],main=colnames(vgibbs)[1],log='y')
      plot(max(1,g-1000):g,vgibbs[max(1,g-1000):g,2],type='l',xlab='gibbs step',
           ylab=colnames(vgibbs)[2],main=colnames(vgibbs)[2],log='y')   
      plot(max(1,g-1000):g,vgibbs[max(1,g-1000):g,3],type='l',xlab='gibbs step',
           ylab=colnames(vgibbs)[3],main=colnames(vgibbs)[3],log='y')   
  
  if(mean(agibbs[(g-39):g,],na.rm=T)>5)
    boxplot(as.data.frame(agibbs[(g-39):g,]/ggibbs[(g-39):g,]),labels=species,
	  main='lambda:Gref')
  if(mean(agibbs[(g-39):g,],na.rm=T)<5)
    boxplot(as.data.frame(agibbs[(g-39):g,]),labels=species,
	  main='lambda')

    if(nspec>1)
	  boxplot(as.data.frame(vgibbs[(g-39):g,grep('beta',colnames(vgibbs))]),labels=species,
	    main='beta')
    if(nspec==1)
	  plot(max(1,g-1000):g, vgibbs[max(1,g-1000):g,grep('beta',colnames(vgibbs))],type='l',ylab='beta',
	    main='beta')	
      lines((g-39):g,vgibbs[(g-39):g,grep('beta',colnames(vgibbs))],col='red')
      		
  }
	  
	  
cov.plot <- function(){
  for(j in 1:length(year)){
    ly  <- seq(1980,2020,by=4)
    ykeep <- which(trunc(timeall)==year[j])
    tmp <- which(trunc(((timeall[ykeep]-trunc(timeall[ykeep]))*365-
           trunc((timeall[ykeep]-trunc(timeall[ykeep]))*365))*24)==12)
    if(year[j] %in% ly)
    tmp <- which(trunc(((timeall[ykeep]-trunc(timeall[ykeep]))*366-
           trunc((timeall[ykeep]-trunc(timeall[ykeep]))*366))*24)==12)
    if(j==1) mid <- tmp
    if(j>1)  mid <- c(mid,tmp)
    }
  
  par(mfrow=c(1,1),family='serif')
  pairs(cbind(D[mid,1],Q[mid,1],M[mid,1]),labels=c(expression(italic(D[t])),expression(italic(Q[t])),expression(italic(M[jt]))))
  cov.cor <- cor(cbind(D[mid,1],Q[mid,1],M[mid,1]))
  
  cov.cor

  }
  
#boundary line analysis
bla <- function(minQ, minM){

  #partitions for light
  Qpart <- matrix(c(minQ,
                    2000),1,2)

  #partitions for moisture
  Mpart <- matrix(c(  minM,
                      1.00),1,2)

  #partitions for vapor pressure deficit
  tmpD  <- seq(0.6,4.8,by=0.2)
  Dpart <- matrix(c(tmpD,tmpD+0.2),length(tmpD),2)

  setindex <- expand.grid(1:dim(Qpart)[1],1:dim(Mpart)[1],1:length(tmpD))

  #creates vector keep with indices for subset b and individual i
  setfunc <- function(b,tind){
    keep <- which(Q >= Qpart[setindex[b,1],1] & Q <= Qpart[setindex[b,1],2] &
                  M >= Mpart[setindex[b,2],1] & M <= Mpart[setindex[b,2],2] &
                  D >= Dpart[setindex[b,3],1] & D <= Dpart[setindex[b,3],2])
    keep <- keep[is.finite(Jdata[keep,tind])]
  
    keep
    }

  library(outliers)
  
  Gs <- Jdata/e.qt[,specindex]  
  
  #means and standard deviations for each bin
   meanGs <- matrix(NA, dim(setindex)[1], nprobe)
     sdGs <- meanGs  
      nGs <- meanGs

    BL <- matrix(NA,0,3)
   
    for(tree in 1:nprobe){
      for(j in 1:dim(setindex)[1]){
        keep <- setfunc(j,tree)
	    if(length(keep)==0) next
	    meanGs[j,tree] <- mean(Gs[keep,tree])
	      sdGs[j,tree] <- sd(Gs[keep,tree])
	       nGs[j,tree] <- length(keep)
	    if(length(keep)>3 & length(keep)<30){
	    
		
		  nk <- 1
		  while(nk==1){
		    out <- dixon.test(Gs[keep,tree])
		  
		    out.rem <- as.integer(strsplit(out$alternative," ")[[1]][1]=='highest')
	      
		    nokeep <- max(order(Gs[keep,tree]))
		    if(out.rem==0) nokeep <- min(order(Gs[keep,tree]))
		
		    #remove worst outlier
		    if(out$p.value<=.05) {
		      keep <- keep[-nokeep]
		      }
		    #if on outliers, move on
		    if(out$p.value>.05) nk <- 0
		    if(length(keep)<3)  nk <- 0
		  
		    }  
		  }
	  
	    if(length(keep)<3) next
	    BLkeep <- which(Gs[keep,tree]>(meanGs[j,tree]+sdGs[j,tree]))
      
	    #calculate mean boundary only if n>=5
	    if(length(BLkeep)<5) next
	    BL <- rbind(BL,cbind(Dpart[setindex[j,3],1]+.1,tree,Gs[keep[BLkeep],tree]))
	    }
      }

  cf <- matrix(NA,nprobe,3)
    colnames(cf) <- c('gref','lam','ratio')

  for(j in 1:nprobe) {
    if(length(BL[BL[,2]==j,3])==0) next
	out <- lm(BL[BL[,2]==j,3] ~ log(BL[BL[,2]==j,1]))$coef
    cf[j,1:2] <- out
    cf[j,3] <- -cf[j,2]/cf[j,1]
    }  

  write.csv(cbind(probe,cf),file=paste('BLAoutput',year,'csv',sep='.'),
    quote=F,row.names=F)

}  

source("Code/SF.Priors.r")
